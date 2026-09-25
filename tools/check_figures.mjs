#!/usr/bin/env node
/**
 * tools/check_figures.mjs — 像素级检查图不是空白
 *
 * **为什么需要它。** 一个"成功生成"的 PNG 完全可能是空白的：坐标轴画了
 * 但数据点没画、图例占了整个画布、或者数据全被过滤掉导致散点为空。
 * 这类问题在 CI 里表现为绿色 job + 一个 0 字节以外的正常文件，
 * 只有打开看才发现是白图。
 *
 * 做法：解 PNG（zlib inflate + 逐行反滤波），统计非背景像素占比。
 * 纯 JS，不依赖任何图形库。
 *
 * 判定：非背景像素占比 < MIN_INK 判为空白。
 * 阈值不能设太高 —— 一个只有坐标轴的图也有少量墨迹，那正是要抓的。
 *
 * ── 双向检查（S2-3，2026-09-24 新增）────────────────────────────
 *
 * 原来只有下限（挡"白图"），挡不住另一个方向：**墨迹铺满、糊成一片**
 * 同样读不出信息（色标压死、背景填满、热图全黑），而它在下限眼里完美。
 *
 * 两个方向都用**实测标定**，不是拍脑袋：
 *   对图库 161 张实测（geo 53 / scrna 33 / spatial 75）：
 *     min 2.39% / p05 4.26% / 中位 24.12% / p95 64.79% / **max 87.71%**
 *   最密的那张是 `03-03-03-unit3-he-reference`（H&E 组织学照片），
 *   本来就该铺满 —— 所以上限必须留出余量，不能按"中位数"卡。
 *
 * 分档：
 *   FAIL  < MIN_INK(0.2%)   或  > INK_FAIL_MAX(96%)  —— 空白 / 糊死，结构性缺陷
 *   WARN  < INK_WARN_MIN(2%) 或 > INK_WARN_MAX(92%)  —— 提示，不判红
 *
 * 实测标定下当前 161 张**零误判**（最低 2.39% / 最高 87.71%，
 * 距 FAIL 两侧各留 ≥2.3pp 与 ≥8pp 余量）。
 *
 * ── 内容贴边检查（Q-26 / E-49，2026-09-25 新增）──────────────────
 *
 * 墨迹检查有一个盲区：**被裁掉的图照样有墨**。`02-08-01` 的标题两侧
 * 各被切掉一段（左端只剩 `dKnk:`、右端断在 `the 60`），墨迹占比完全正常。
 *
 * 但裁切有一个物理后果：**墨迹一直延伸到画布边缘**。所以量非背景像素的
 * 外接框到四边的距离，左右任一 < EDGE_MIN_PX 就判红。
 *
 * 阈值同样是实测标定，不是拍脑袋：
 *   scrna 34 张：min(L,R) 分布 {0:1, 10:9, 11:2, 12:6, 13:13, 14:2, 41:1}
 *   spatial 76 张：min(L,R) 分布 {10:4, 12:13, 13:38, 14:13, 15:1, 16:3, 20:3, 41:1}
 *   → 110 张里唯一 < 10 px 的就是被裁的 `02-08-01`（L0/R0），
 *     阈值取 3 px **零误伤**。
 *
 * **只看左右。** 实测 `02-05-05-unit1-pseudotime-by-cluster` 的上边距是 0
 * —— 那是布局取舍（子图顶到边），不是内容装不下；横向贴边才是
 * "内容比画布宽、被 savefig 切掉"的信号。
 * 暗底图（照片类整幅都是"墨"）外接框必是满幅，量不出边距，跳过。
 *
 * 用法: node tools/check_figures.mjs <目录> [更多目录...]
 */

import { readdirSync, readFileSync, existsSync, statSync, writeFileSync } from "node:fs";
import { join } from "node:path";
import { inflateSync } from "node:zlib";

// 非背景像素占比低于此值判为空白
const MIN_INK = 0.002;
// 双向检查的提示/判红档（见文件头标定说明）
const INK_WARN_MIN = 0.02;
const INK_WARN_MAX = 0.92;
const INK_FAIL_MAX = 0.96;
// 内容贴边判红阈值（px）：非背景像素外接框到画布左/右缘的距离下限。
// 实测 110 张图（scrna 34 + spatial 76）里，除被裁的 02-08-01（0 px）外
// 最小的也有 10 px —— 3 px 零误伤。
const EDGE_MIN_PX = 3;

function decodePng(buf) {
  if (buf.readUInt32BE(0) !== 0x89504e47) throw new Error("不是 PNG");
  let pos = 8;
  let width = 0, height = 0, bitDepth = 0, colorType = 0;
  const idat = [];
  while (pos < buf.length) {
    const len = buf.readUInt32BE(pos);
    const type = buf.toString("ascii", pos + 4, pos + 8);
    const data = buf.subarray(pos + 8, pos + 8 + len);
    if (type === "IHDR") {
      width = data.readUInt32BE(0);
      height = data.readUInt32BE(4);
      bitDepth = data[8];
      colorType = data[9];
    } else if (type === "IDAT") {
      idat.push(data);
    } else if (type === "IEND") {
      break;
    }
    pos += 12 + len;
  }
  if (bitDepth !== 8) throw new Error(`只支持 8 位深度，实际 ${bitDepth}`);
  const channels = { 0: 1, 2: 3, 3: 1, 4: 2, 6: 4 }[colorType];
  if (!channels) throw new Error(`不支持的 colorType ${colorType}`);

  const raw = inflateSync(Buffer.concat(idat));
  const stride = width * channels;
  const out = Buffer.alloc(height * stride);

  // 逐行反滤波（PNG 的 5 种滤波器）
  for (let y = 0; y < height; y++) {
    const ft = raw[y * (stride + 1)];
    const src = raw.subarray(y * (stride + 1) + 1, y * (stride + 1) + 1 + stride);
    const dst = out.subarray(y * stride, (y + 1) * stride);
    const prev = y > 0 ? out.subarray((y - 1) * stride, y * stride) : null;
    for (let x = 0; x < stride; x++) {
      const a = x >= channels ? dst[x - channels] : 0;
      const b = prev ? prev[x] : 0;
      const c = prev && x >= channels ? prev[x - channels] : 0;
      let v = src[x];
      if (ft === 1) v += a;
      else if (ft === 2) v += b;
      else if (ft === 3) v += (a + b) >> 1;
      else if (ft === 4) {
        const p = a + b - c;
        const pa = Math.abs(p - a), pb = Math.abs(p - b), pc = Math.abs(p - c);
        v += pa <= pb && pa <= pc ? a : pb <= pc ? b : c;
      }
      dst[x] = v & 0xff;
    }
  }
  return { width, height, channels, pixels: out };
}

function inkFraction(png) {
  const { width, height, channels, pixels } = png;
  // 取四角的多数色作为背景色（图通常白底，但不要假定）
  const corner = (x, y) => {
    const o = (y * width + x) * channels;
    return [pixels[o], pixels[o + 1], pixels[o + 2]];
  };
  const cs = [corner(0, 0), corner(width - 1, 0), corner(0, height - 1),
              corner(width - 1, height - 1)];
  const bg = [0, 1, 2].map((i) => {
    const v = cs.map((c) => c[i]).sort((a, b) => a - b);
    return v[Math.floor(v.length / 2)];
  });

  let ink = 0;
  const total = width * height;
  // 非背景像素的外接框（Q-26 / E-49：内容贴边 = 被裁的信号）
  let minX = width, maxX = -1, minY = height, maxY = -1;
  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      const o = (y * width + x) * channels;
      const d = Math.abs(pixels[o] - bg[0]) + Math.abs(pixels[o + 1] - bg[1]) +
                Math.abs(pixels[o + 2] - bg[2]);
      if (d > 30) {
        ink++;
        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
      }
    }
  }
  // 背景亮度一起返回：全黑图会让"四角取背景"把黑当成背景、算出 0% 墨迹，
  // 于是**糊死被误报成空白**。判据是相对量，必须再报一个绝对量才分得清
  // （与 AGENTS 里"检查器要挡住自己的盲区"同一条）。
  const bgLuma = (bg[0] + bg[1] + bg[2]) / 3;
  const margins = maxX < 0 ? null : {
    left: minX, right: width - 1 - maxX,
    top: minY, bottom: height - 1 - maxY,
  };
  return { frac: ink / total, bgLuma, margins };
}

const dirs = process.argv.slice(2);
if (dirs.length === 0) {
  console.error("用法: node tools/check_figures.mjs <目录> [更多目录...]");
  process.exit(1);
}

let n = 0, blank = 0;
const warnings = [];
const edgeHits = [];
for (const dir of dirs) {
  if (!existsSync(dir)) continue;
  for (const e of readdirSync(dir, { withFileTypes: true })) {
    if (!e.isFile() || !e.name.endsWith(".png")) continue;
    const p = join(dir, e.name);
    n++;
    const size = statSync(p).size;
    try {
      const png = decodePng(readFileSync(p));
      const { frac, bgLuma, margins } = inkFraction(png);
      // 暗背景（亮度 < 128）说明图不是白底 —— 此时"四角取背景"的墨迹占比
      // 不可解释，直接判红让人去看，而不是给出一个可能反过来的结论。
      const darkBg = bgLuma < 128;
      if (frac < MIN_INK) {
        console.error(`  [${darkBg ? "糊死" : "空白"}] ${p}  ${png.width}x${png.height}  ` +
                      `墨迹占比 ${(frac * 100).toFixed(4)}%  背景亮度 ${bgLuma.toFixed(0)}  (${size} B)` +
                      (darkBg ? "  ← 整幅暗底，等同糊死" : ""));
        blank++;
      } else if (frac > INK_FAIL_MAX) {
        console.error(`  [糊死] ${p}  ${png.width}x${png.height}  ` +
                      `墨迹占比 ${(frac * 100).toFixed(2)}% > ${(INK_FAIL_MAX * 100).toFixed(0)}%  ` +
                      `(色标压死 / 背景填满，图读不出信息)`);
        blank++;
      } else {
        // 内容贴边（Q-26 / E-49）：墨迹顶到左右缘 = 内容比画布宽、被切了。
        // 暗底图整幅都是"墨"，外接框必满幅，量不出边距 —— 跳过。
        const edgeBad = !darkBg && margins &&
                        (margins.left < EDGE_MIN_PX || margins.right < EDGE_MIN_PX);
        if (edgeBad) {
          console.error(`  [贴边] ${e.name}  ${png.width}x${png.height}  ` +
                        `左右边距 ${margins.left}/${margins.right} px ` +
                        `< ${EDGE_MIN_PX} px  ← 内容比画布宽，两侧被 savefig 切掉了`);
          blank++;
          edgeHits.push({ figure: e.name, left: margins.left, right: margins.right });
        } else {
          console.log(`  [OK]   ${e.name}  ${png.width}x${png.height}  ` +
                      `墨迹 ${(frac * 100).toFixed(2)}%` +
                      (margins ? `  边距 L/R ${margins.left}/${margins.right}` : ""));
          // 双向提示：超出常见范围但不判红
          if (frac < INK_WARN_MIN || frac > INK_WARN_MAX) {
            warnings.push(`${e.name}  ${(frac * 100).toFixed(2)}%`);
          }
        }
      }
    } catch (err) {
      console.error(`  [错误] ${p}: ${err.message}`);
      blank++;
    }
  }
}

if (warnings.length > 0) {
  console.log(`\n[WARN] ${warnings.length} 张图墨迹占比超出常见范围 ` +
              `[${(INK_WARN_MIN * 100).toFixed(0)}%, ${(INK_WARN_MAX * 100).toFixed(0)}%] —— 仅提示：`);
  for (const w of warnings) console.log(`    ${w}`);
  console.log("  偏疏多为散点/森林图（本来就少），偏密多为热图/组织学照片（本来就满）。");
  console.log("  只要不是「空白」或「糊死」，具体疏密由人工终审判断。");
}

// ---- WARN 落盘（2026-09-24 审计 P1-9：WARN 必须有稳定消费入口，否则等于噪声）----
// 与 geo 仓 check_fig_sizes.mjs 同款约定：写 warn_report.json，
// ① 人工亲读图**之前**先看（selfcheck 汇总提示）；② 随 artifact 上传，跨轮对比
// "WARN 集合是否稳定"——稳定 = 已知审美取舍，新出现 = 回归信号。
// **注意本文件与 spatial 仓的 check_figures.mjs 必须逐字相同**（既有约定）。
if (warnings.length > 0) {
  const report = {
    generatedAt: new Date().toISOString(),
    tool: "check_figures.mjs",
    gates: { inkFailBlank: INK_FAIL_MIN, inkFailSaturated: INK_FAIL_MAX,
            inkWarnBand: [INK_WARN_MIN, INK_WARN_MAX] },
    summary: { total: n, blank: blank, warn: warnings.length },
    items: warnings.map((w) => {
      const m = w.match(/^(\S+)\s+([\d.]+)%$/);
      const frac = m ? Number(m[2]) / 100 : null;
      return {
        figure: m ? m[1] : w,
        kind: "ink-fraction-warn",
        value: frac,
        band: [INK_WARN_MIN, INK_WARN_MAX],
        hint: frac !== null && frac < INK_WARN_MIN
          ? "墨迹偏疏——多为散点/森林图（本来就少），人工确认是否真空白"
          : "墨迹偏密——多为热图/组织学照片（本来就满），人工确认是否糊死",
      };
    }),
  };
  try {
    writeFileSync(join(dir, "warn_report.json"), JSON.stringify(report, null, 2), "utf8");
    console.log(`\n[WARN 报告] 已写出 ${join(dir, "warn_report.json")}（${warnings.length} 条）`);
  } catch (e) {
    console.log(`\n[WARN 报告] 写出失败（不阻断）：${e.message}`);
  }
}

if (n === 0) {
  console.error("没有找到任何 PNG —— 图目录不对或图没生成");
  process.exit(1);
}
if (blank > 0) {
  if (edgeHits.length > 0) {
    console.error(`\n其中 ${edgeHits.length} 张内容贴到画布左右缘（边距 < ${EDGE_MIN_PX} px）` +
                  `—— 内容比画布宽，两侧已被 savefig 静默切掉：`);
    for (const h of edgeHits) {
      console.error(`    ${h.figure}  L${h.left}/R${h.right} px`);
    }
    console.error("  修法：调大 figsize 宽度、或把标题/图例改短（折行）。" +
                  "注意宽度仍要落在 89/136/183 mm 三档之内。");
  }
  console.error(`\n${blank}/${n} 张图有问题`);
  process.exit(1);
}
console.log(`\n${n} 张图全部非空白、非糊死、内容不贴边`);
