#!/usr/bin/env node
/**
 * 检查出图的分辨率是否达到投稿底线。
 *
 *     node tools/check_fig_dpi.mjs results/pbmc3k/figures
 *     node tools/check_fig_dpi.mjs results/lymph_node/figures
 *
 * **为什么需要这个工具：**
 *
 * 分辨率是"静默失败"的典型 —— 图照样生成、文件大小正常、`check_figures.mjs`
 * 照样报"有墨"，只有放大或印刷时才发现字形和细线发虚。而发虚的原因
 * 往往不是画图代码，是**配置把默认值覆盖掉了**：
 *
 *     common.py:  ana.setdefault("figure_dpi", 300)     # 只在键缺失时生效
 *     config.yml: figure_dpi: 150                       # <- 这个赢了
 *
 * 实测就踩过：只改了 `common.py` 的默认值，CI 跑绿，日志里没有任何异常，
 * 但从 artifact 里读 PNG 头发现宽度仍是 907 px（= 6.05 in × 150）。
 * **只有把 artifact 拿下来量像素才能发现。**
 *
 * **怎么算 dpi：**
 *
 * PNG 的像素数 = figsize(英寸) × dpi，而 PDF 的 `/MediaBox` 就是
 * figsize(英寸) × 72（点）。两者相除，figsize 约掉，得到真实 dpi：
 *
 *     dpi = png_px / (MediaBox_pt / 72)
 *
 * 所以**不需要**在配置里写死"这张图应该是多少毫米"—— 物理尺寸从 PDF 读，
 * 像素数从 PNG 读，两边独立，算出来的 dpi 是实测值。
 *
 * 判据是**两条边都要达标**。只查宽度会漏掉"高度被压扁"的情况
 * （`figsize` 的高度是另一个表达式，可以各自出错）。
 *
 * 退出码：0 = 全部达标；1 = 有图不达标或读不出来。
 * 纯 Node、零依赖。
 */

import { readFileSync, readdirSync, statSync } from "node:fs";
import { join, basename } from "node:path";

/** 投稿底线。低于它就要么重出图，要么承认这张图不能印刷。 */
const MIN_DPI = 300;
/** 允许的舍入余量：matplotlib 写 MediaBox 时可能少几个 ulp。 */
const TOL = 1.0;
const PT_PER_IN = 72;

function readPngSize(buf) {
  if (buf.readUInt32BE(0) !== 0x89504e47) throw new Error("不是 PNG");
  // IHDR 必须是第一个 chunk：8 字节签名 + 4 长度 + 4 类型 = 偏移 16
  const w = buf.readUInt32BE(16);
  const h = buf.readUInt32BE(20);
  if (!w || !h) throw new Error("IHDR 尺寸为 0");
  return { w, h };
}

/**
 * 从 PDF 里读第一个 /MediaBox。
 *
 * 用 latin1 读：PDF 的交叉引用表里可能有二进制字节，按 utf8 解码会
 * 替换成 U+FFFD 并**改变字节偏移**，而 MediaBox 是 ASCII 文本、
 * 用 latin1 逐字节保留原样，正则照样能匹配。
 */
function readMediaBox(buf) {
  const txt = buf.toString("latin1");
  const m = /\/MediaBox\s*\[\s*([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)\s*\]/.exec(txt);
  if (!m) throw new Error("读不到 /MediaBox");
  const [, x0, y0, x1, y1] = m.map(Number);
  const w = Math.abs(x1 - x0);
  const h = Math.abs(y1 - y0);
  if (!w || !h) throw new Error("MediaBox 尺寸为 0");
  return { w, h };
}

function main() {
  const dir = process.argv[2];
  if (!dir) {
    console.error("用法: node tools/check_fig_dpi.mjs <图目录>");
    console.error("  例: node tools/check_fig_dpi.mjs results/pbmc3k/figures");
    process.exit(1);
  }
  let entries;
  try {
    entries = readdirSync(dir);
  } catch (e) {
    console.error(`读不到目录: ${dir}`);
    process.exit(1);
  }

  // 跳过 `__` 开头的临时/中间文件 —— 它们不是交付图。
  const pngs = entries.filter((f) => f.endsWith(".png") && !f.startsWith("__")).sort();
  if (pngs.length === 0) {
    console.error(`[FAIL] ${dir} 里一张 PNG 都没有 —— 这个检查等于没跑`);
    process.exit(1);
  }

  const rows = [];
  const problems = [];

  for (const f of pngs) {
    const pngPath = join(dir, f);
    const pdfPath = join(dir, f.replace(/\.png$/, ".pdf"));
    const name = basename(f, ".png");

    let png;
    try {
      png = readPngSize(readFileSync(pngPath));
    } catch (e) {
      problems.push(`${name}: PNG 读不出来（${e.message}）`);
      continue;
    }

    let box;
    try {
      statSync(pdfPath);
      box = readMediaBox(readFileSync(pdfPath));
    } catch (e) {
      // **没有 PDF 兄弟文件 = 无法验证物理尺寸 = 不能算通过。**
      // 静默跳过会让"少出一份 PDF"这种问题从检查里溜走。
      problems.push(`${name}: 没有可读的 PDF 兄弟文件（${e.message}），无法验证 dpi`);
      continue;
    }

    const dpiX = png.w / (box.w / PT_PER_IN);
    const dpiY = png.h / (box.h / PT_PER_IN);
    const mmW = (box.w / PT_PER_IN) * 25.4;
    const mmH = (box.h / PT_PER_IN) * 25.4;
    rows.push({ name, png, mmW, mmH, dpiX, dpiY });
    if (dpiX < MIN_DPI - TOL || dpiY < MIN_DPI - TOL) {
      problems.push(
        `${name}: dpi ${dpiX.toFixed(0)}x${dpiY.toFixed(0)} < ${MIN_DPI}` +
          `（${mmW.toFixed(1)}x${mmH.toFixed(1)} mm -> ${png.w}x${png.h} px）`
      );
    }
  }

  console.log(`图目录: ${dir}`);
  console.log(
    "  " +
      "图".padEnd(34) +
      "物理尺寸(mm)".padEnd(20) +
      "像素".padEnd(14) +
      "dpi"
  );
  for (const r of rows) {
    console.log(
      "  " +
        r.name.padEnd(34) +
        `${r.mmW.toFixed(1)}x${r.mmH.toFixed(1)}`.padEnd(20) +
        `${r.png.w}x${r.png.h}`.padEnd(14) +
        `${r.dpiX.toFixed(0)}x${r.dpiY.toFixed(0)}`
    );
  }

  const buckets = new Map();
  for (const r of rows) {
    const k = `${Math.round(r.dpiX / 50) * 50}`;
    buckets.set(k, (buckets.get(k) || 0) + 1);
  }
  console.log("  dpi 分布:");
  for (const k of [...buckets.keys()].sort((a, b) => a - b)) {
    console.log(`    ~${k.padStart(4)}      ${buckets.get(k)}`);
  }

  if (problems.length) {
    console.error("");
    for (const p of problems) console.error(`  [FAIL] ${p}`);
    console.error(`\n[FAIL] ${problems.length} 张图未达到 ${MIN_DPI} dpi（或无法验证）`);
    process.exit(1);
  }
  console.log(`\n全部 ${rows.length} 张图都达到 ${MIN_DPI} dpi`);
}

main();
