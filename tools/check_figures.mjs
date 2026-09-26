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
 *   对图库实测（geo 53 / scrna 33 / spatial 75，S2-3 当时的快照）：
 *     min 2.39% / p05 4.26% / 中位 24.12% / p95 64.79% / max 87.71%
 *   最密的那张是 `03-03-03-unit3-he-reference`（H&E 组织学照片），
 *   本来就该铺满 —— 所以上限必须留出余量，不能按"中位数"卡。
 *
 *   **R-04c 之后那张 H&E 图升到 94.24%**（修好等比例后组织铺满画布，
 *   见 `governance/15_ERROR_LEDGER.md` E-63）。它落在 WARN 带里、不判红，
 *   但**正是它第一次触发了下面那个从未执行过的 WARN 落盘分支** ——
 *   把两个潜伏已久的缺陷炸了出来。
 *
 * 分档：
 *   FAIL  < MIN_INK(0.2%)   或  > INK_FAIL_MAX(96%)  —— 空白 / 糊死，结构性缺陷
 *   WARN  < INK_WARN_MIN(2%) 或 > INK_WARN_MAX(92%)  —— 提示，不判红
 *
 * 实测标定下当时 161 张**零误判**（最低 2.39% / 最高 87.71%，
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
 * ── WARN 落盘（P1-9，2026-09-24 新增；E-63 修好，2026-09-26）──────
 *
 * WARN 必须有稳定消费入口，否则等于噪声：写 `warn_report.json`
 * （路径与 geo 仓 `check_fig_sizes.mjs` 同款约定 = `dirname(argv[2])/warn_report.json`），
 * ① 人工亲读图**之前**先看；② 随 artifact 上传，跨轮对比"WARN 集合是否稳定"
 * ——稳定 = 已知审美取舍，新出现 = 回归信号。
 *
 * **这段代码从落地到 E-63 之间一次都没执行过**（WARN 带一直没被触发），
 * 于是两个缺陷在里面潜伏了很久：`INK_FAIL_MIN` 根本没声明（真常量叫
 * `MIN_INK`），以及 `dir` 在循环外已离开作用域、报告**永远写不出来**。
 * 第一个让门禁直接 `ReferenceError` 崩掉（假红），第二个被 `catch` 吞成
 * "写出失败（不阻断）"（假绿）。教训：**只在 WARN 分支里的代码需要
 * 一个能强制走到它的自检** —— 见下面的 `--selftest`。
 *
 * 用法:
 *   node tools/check_figures.mjs <目录> [更多目录...]
 *   node tools/check_figures.mjs --selftest
 */

import { readdirSync, readFileSync, existsSync, statSync, writeFileSync, mkdtempSync, rmSync } from "node:fs";
import { join, dirname } from "node:path";
import { inflateSync, deflateSync } from "node:zlib";
import { tmpdir } from "node:os";

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

// ---- 核心：扫描目录集合，返回结构化结果（不 print、不 exit）------------------
// 抽成函数是为了 `--selftest` 能直接断言。**这是 E-63 的教训**：
// 只活在 WARN 分支里的代码需要一个能强制走到它的入口。
function checkDirs(dirs, { log = console.log, err = console.error } = {}) {
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
          err(`  [${darkBg ? "糊死" : "空白"}] ${p}  ${png.width}x${png.height}  ` +
              `墨迹占比 ${(frac * 100).toFixed(4)}%  背景亮度 ${bgLuma.toFixed(0)}  (${size} B)` +
              (darkBg ? "  ← 整幅暗底，等同糊死" : ""));
          blank++;
        } else if (frac > INK_FAIL_MAX) {
          err(`  [糊死] ${p}  ${png.width}x${png.height}  ` +
              `墨迹占比 ${(frac * 100).toFixed(2)}% > ${(INK_FAIL_MAX * 100).toFixed(0)}%  ` +
              `(色标压死 / 背景填满，图读不出信息)`);
          blank++;
        } else {
          // 内容贴边（Q-26 / E-49）：墨迹顶到左右缘 = 内容比画布宽、被切了。
          // 暗底图整幅都是"墨"，外接框必满幅，量不出边距 —— 跳过。
          const edgeBad = !darkBg && margins &&
                          (margins.left < EDGE_MIN_PX || margins.right < EDGE_MIN_PX);
          if (edgeBad) {
            err(`  [贴边] ${e.name}  ${png.width}x${png.height}  ` +
                `左右边距 ${margins.left}/${margins.right} px ` +
                `< ${EDGE_MIN_PX} px  ← 内容比画布宽，两侧被 savefig 切掉了`);
            blank++;
            edgeHits.push({ figure: e.name, left: margins.left, right: margins.right });
          } else {
            log(`  [OK]   ${e.name}  ${png.width}x${png.height}  ` +
                `墨迹 ${(frac * 100).toFixed(2)}%` +
                (margins ? `  边距 L/R ${margins.left}/${margins.right}` : ""));
            // 双向提示：超出常见范围但不判红
            if (frac < INK_WARN_MIN || frac > INK_WARN_MAX) {
              warnings.push(`${e.name}  ${(frac * 100).toFixed(2)}%`);
            }
          }
        }
      } catch (err2) {
        err(`  [错误] ${p}: ${err2.message}`);
        blank++;
      }
    }
  }
  return { n, blank, warnings, edgeHits };
}

// WARN 集合 → 报告对象。**gates 里引用的每个常量都必须真实存在** ——
// 第一版这里写的是从未声明的 `INK_FAIL_MIN`，而它只在 WARN 分支里求值，
// 于是从落地起潜伏到 E-63 才炸（那之前 WARN 带一次都没被触发过）。
function buildWarnReport({ n, blank, warnings }) {
  return {
    generatedAt: new Date().toISOString(),
    tool: "check_figures.mjs",
    gates: { inkFailBlank: MIN_INK, inkFailSaturated: INK_FAIL_MAX,
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
}

// WARN 落盘（P1-9）：路径约定与 geo 仓 check_fig_sizes.mjs 一致 ——
// `dirname(argv[2])/warn_report.json`，即与 `figures/` 同级、随 artifact 上传。
//
// **抽成函数是为了自检能走到真代码。** E-63 的第二半就是这段的旧写法：
// 它引用了上面 for 循环的循环变量 `dir`，出了循环已离开作用域 →
// `dir is not defined`，被 catch 吞成"写出失败（不阻断）"→ 报告永远不存在
// 而门禁照样绿。
//
// 自检用例必须调**这个函数**，不能在自检里另写一遍路径拼接 —— 那样
// 把 bug 注回去自检也不会红（我第一版就是这么写的，反向标定当场抓到）。
function writeWarnReport(dirs, payload) {
  // `dirname("a/b/figures/")` 会给出 `a/b/figures`（尾斜杠把最后一段当文件名）——
  // 先剥掉尾分隔符，否则报告会落进 figures/ 里面、与图混在一起。
  const first = String(dirs[0]).replace(/[\\/]+$/, "");
  const outPath = join(dirname(first), "warn_report.json");
  writeFileSync(outPath, JSON.stringify(buildWarnReport(payload), null, 2), "utf8");
  return outPath;
}

// ---- 自检：用合成 PNG 强制走到每个分支 ---------------------------------------
// 合成图不依赖 matplotlib，也不依赖任何真实产物 —— 所以 CI 和本机行为一致。
const CRC_TABLE = (() => {
  const t = new Int32Array(256);
  for (let i = 0; i < 256; i++) {
    let c = i;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    t[i] = c;
  }
  return t;
})();
function crc32(buf) {
  let c = -1;
  for (let i = 0; i < buf.length; i++) c = CRC_TABLE[(c ^ buf[i]) & 0xff] ^ (c >>> 8);
  return (c ^ -1) >>> 0;
}
function encodePng(width, height, rgbAt) {
  const raw = Buffer.alloc(height * (width * 3 + 1));
  let o = 0;
  for (let y = 0; y < height; y++) {
    raw[o++] = 0; // filter: none
    for (let x = 0; x < width; x++) {
      const [r, g, b] = rgbAt(x, y);
      raw[o++] = r; raw[o++] = g; raw[o++] = b;
    }
  }
  const chunk = (type, data) => {
    const len = Buffer.alloc(4);
    len.writeUInt32BE(data.length, 0);
    const td = Buffer.concat([Buffer.from(type, "ascii"), data]);
    const crc = Buffer.alloc(4);
    crc.writeUInt32BE(crc32(td), 0);
    return Buffer.concat([len, td, crc]);
  };
  const ihdr = Buffer.alloc(13);
  ihdr.writeUInt32BE(width, 0);
  ihdr.writeUInt32BE(height, 4);
  ihdr[8] = 8; ihdr[9] = 2; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;
  return Buffer.concat([
    Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]),
    chunk("IHDR", ihdr),
    chunk("IDAT", deflateSync(raw)),
    chunk("IEND", Buffer.alloc(0)),
  ]);
}

// 白底 + 指定比例的深色像素（按行填充），可控地落进任一档
function synth(width, height, inkFrac, { margin = 6 } = {}) {
  const inner = width - 2 * margin;
  const need = Math.round(inkFrac * width * height);
  return encodePng(width, height, (x, y) => {
    if (x < margin || x >= width - margin || y < margin || y >= height - margin) {
      return [255, 255, 255];
    }
    const idx = (y - margin) * inner + (x - margin);
    return idx < need ? [40, 40, 40] : [255, 255, 255];
  });
}
// 内容顶到左右缘（模拟被 savefig 裁切）：几条横贯整幅的细线，
// 墨迹占比很低（不会先落进"糊死"档），但外接框左右边距为 0。
function synthEdge(width, height) {
  return encodePng(width, height, (_x, y) => (y >= 50 && y < 55 ? [40, 40, 40] : [255, 255, 255]));
}
// 整幅暗底（模拟"糊死"，四角取背景会把黑当背景 → 必须靠 bgLuma 分出来）
function synthDark(width, height) {
  return encodePng(width, height, () => [10, 10, 10]);
}

function selftest() {
  const dir = mkdtempSync(join(tmpdir(), "figselftest-"));
  const quiet = { log: () => {}, err: () => {} };
  const cases = [];
  const put = (name, buf) => writeFileSync(join(dir, name), buf);
  const clear = () => {
    for (const e of readdirSync(dir)) rmSync(join(dir, e), { force: true });
  };
  const scan = () => checkDirs([dir], quiet);

  try {
    // 用例 1：正常图（墨迹约 24%，落在 WARN 带内）→ 零 WARN 零红
    put("a-normal.png", synth(200, 200, 0.24));
    let r = scan();
    cases.push(["正常图零 WARN", r.warnings.length === 0 && r.blank === 0,
                `warn=${r.warnings.length} blank=${r.blank}`]);

    // 用例 2：偏密（94%，就是 R-04c 后那张 H&E 的档位）→ 1 条 WARN、不判红
    clear();
    put("b-dense.png", synth(400, 400, 0.94));
    r = scan();
    cases.push(["94% 落 WARN 带且不判红", r.warnings.length === 1 && r.blank === 0,
                `warn=${r.warnings.length} blank=${r.blank}`]);

    // 用例 3：偏疏（0.4%）→ 1 条 WARN、不判红
    clear();
    put("c-sparse.png", synth(400, 400, 0.004));
    r = scan();
    cases.push(["0.4% 落 WARN 带且不判红", r.warnings.length === 1 && r.blank === 0,
                `warn=${r.warnings.length} blank=${r.blank}`]);

    // 用例 4：全白 → 判红（空白）
    clear();
    put("d-blank.png", synth(120, 120, 0.0));
    r = scan();
    cases.push(["全白判红", r.blank >= 1, `blank=${r.blank}`]);

    // 用例 5：整幅暗底 → 判红（糊死，靠 bgLuma 分出来）
    clear();
    put("e-dark.png", synthDark(120, 120));
    r = scan();
    cases.push(["整幅暗底判红", r.blank >= 1, `blank=${r.blank}`]);

    // 用例 6：内容贴到左右缘 → 判红（贴边）
    clear();
    put("f-edge.png", synthEdge(200, 200));
    r = scan();
    cases.push(["内容贴边判红", r.edgeHits.length >= 1, `edge=${r.edgeHits.length}`]);

    // 用例 7（E-63 回归 1）：WARN 报告能构造出来 —— 引用的常量全部存在。
    // 旧代码在这里抛 `ReferenceError: INK_FAIL_MIN is not defined`。
    let reportOk = false, reportDetail = "";
    try {
      const rep = buildWarnReport({ n: 1, blank: 0, warnings: ["x.png  94.00%"] });
      reportOk = rep.gates.inkFailBlank === MIN_INK && rep.gates.inkFailSaturated === INK_FAIL_MAX &&
                 rep.items.length === 1 && rep.items[0].value === 0.94;
      reportDetail = JSON.stringify(rep.gates);
    } catch (e) {
      reportDetail = `${e.constructor.name}: ${e.message}`;
    }
    cases.push(["WARN 报告可构造（E-63 回归 1）", reportOk, reportDetail]);

    // 用例 8（E-63 回归 2）：报告真的能写到盘上 —— 走 main() 用的**同一个函数**。
    // 旧代码用的是循环外的 `dir`，已离开作用域 → 永远写不出来还被 catch 吞掉。
    let writeOk = false, writeDetail = "";
    try {
      const outPath = writeWarnReport([dir], { n: 1, blank: 0, warnings: ["x.png  94.00%"] });
      const back = JSON.parse(readFileSync(outPath, "utf8"));
      writeOk = back.tool === "check_figures.mjs" && back.summary.warn === 1;
      writeDetail = outPath;
    } catch (e) {
      writeDetail = `${e.constructor.name}: ${e.message}`;
    }
    cases.push(["WARN 报告可写盘（E-63 回归 2）", writeOk, writeDetail]);

    // 用例 9：目录参数带尾斜杠时，报告不能落进 figures/ 里面。
    // `dirname("a/figures/")` = `a/figures` —— 尾斜杠把最后一段当成了文件名。
    let tailOk = false, tailDetail = "";
    try {
      const p = writeWarnReport([dir + "/"], { n: 1, blank: 0, warnings: ["x.png  94.00%"] });
      tailOk = dirname(p) === dirname(dir) && !p.startsWith(join(dir, "warn_report"));
      tailDetail = p;
    } catch (e) {
      tailDetail = `${e.constructor.name}: ${e.message}`;
    }
    cases.push(["尾斜杠目录报告不落进 figures/（E-63 回归 3）", tailOk, tailDetail]);

    let failed = 0;
    for (const [name, pass, detail] of cases) {
      if (!pass) failed++;
      console.log(`  ${pass ? "[OK]  " : "[FAIL]"} ${name}${pass ? "" : `  ← ${detail}`}`);
    }
    console.log(`\n自检${failed === 0 ? "通过" : "失败"}（${cases.length} 个用例，含 3 条 E-63 回归）`);
    return failed === 0 ? 0 : 1;
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
}

function main() {
  const dirs = process.argv.slice(2);
  if (dirs.length === 0) {
    console.error("用法: node tools/check_figures.mjs <目录> [更多目录...]");
    console.error("      node tools/check_figures.mjs --selftest");
    process.exit(1);
  }

  const { n, blank, warnings, edgeHits } = checkDirs(dirs);

  if (warnings.length > 0) {
    console.log(`\n[WARN] ${warnings.length} 张图墨迹占比超出常见范围 ` +
                `[${(INK_WARN_MIN * 100).toFixed(0)}%, ${(INK_WARN_MAX * 100).toFixed(0)}%] —— 仅提示：`);
    for (const w of warnings) console.log(`    ${w}`);
    console.log("  偏疏多为散点/森林图（本来就少），偏密多为热图/组织学照片（本来就满）。");
    console.log("  只要不是「空白」或「糊死」，具体疏密由人工终审判断。");
  }

  // WARN 落盘（P1-9）：见 writeWarnReport 的注释（E-63 就出在这里）。
  if (warnings.length > 0) {
    try {
      const outPath = writeWarnReport(dirs, { n, blank, warnings });
      console.log(`\n[WARN 报告] 已写出 ${outPath}（${warnings.length} 条）—— 亲读图前先看这个`);
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
}

if (process.argv.includes("--selftest")) {
  process.exit(selftest());
} else {
  main();
}
