#!/usr/bin/env node
/**
 * 调色板门禁：从 scripts/lib/common.py 里解析出实际色值，按判据验证。
 *
 * 移植自 geo-normal-pipeline-skill/tools/check_palette.mjs（色彩数学逐行相同），
 * 但判据按本仓库的情况重写了 —— Python 侧的色板结构与 R 侧不同。
 *
 * **为什么本仓库需要这道门禁：** 实测 `umap_clusters` 有 10 个簇，
 * 而 `PAL_CYCLE` 只有 5 色，于是 0≡5 / 1≡6 / 2≡7 / 3≡8 / 4≡9 ——
 * 图上簇 5 与簇 0 是**上下相邻的两个同色团块**，读者看不出分界在哪。
 * 这类问题**凭眼睛看不出是"配色不够"还是"聚类错了"**，只有算才看得见。
 *
 * 判据：
 *   1. 色相相差 < 15° 视为同一颜色（近灰色不承载色相语义，跳过）
 *   2. 三种色盲下 OKLab 距离 >= 0.05
 *   3. 白底对比度 >= 2.0（`muted` 例外 —— 它专门用于随机基线这类次要参照）
 *   4. `PAL_CYCLE` 至少 10 色（实测 resolution=1.0 给 10 个簇；
 *      少于簇数就必然有簇撞色，而撞色**不会报错**，只会让两张图看起来一样）
 *   5. **`PAL_CYCLE` 与 `assets/publication.mplstyle` 的 `axes.prop_cycle`
 *      必须逐色相同。** 这两个地方**各定义了一份**，改一处就会出现
 *      "rcParams 说 A、样式文件说 B"，而 `apply_style()` 里
 *      `rcParams["axes.prop_cycle"] = cycler(color=PAL_CYCLE)` 会让
 *      common.py 那一份悄悄赢 —— 样式文件于是变成一句谎话。
 *
 * 用法：node tools/check_palette.mjs
 */
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const here = dirname(fileURLToPath(import.meta.url));
const repo = join(here, "..");
const commonPath = join(repo, "scripts", "lib", "common.py");
const stylePath = join(repo, "assets", "publication.mplstyle");

// ---- 色彩空间（与 geo 侧逐行相同）-----------------------------------------
const srgbToLinear = (c) => (c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);
const linearToSrgb = (c) => {
  const v = Math.min(1, Math.max(0, c));
  return v <= 0.0031308 ? v * 12.92 : 1.055 * v ** (1 / 2.4) - 0.055;
};

function hexToRgb(h) {
  const s = h.replace("#", "");
  return [0, 2, 4].map((i) => parseInt(s.slice(i, i + 2), 16) / 255);
}

const M1 = [
  [0.4122214708, 0.5363325363, 0.0514459929],
  [0.2119034982, 0.6806995451, 0.1073969566],
  [0.0883024619, 0.2817188376, 0.6299787005],
];
const M2 = [
  [0.2104542553, 0.793617785, -0.0040720468],
  [1.9779984951, -2.428592205, 0.4505937099],
  [0.0259040371, 0.7827717662, -0.808675766],
];
const mul = (M, v) => M.map((row) => row.reduce((s, x, i) => s + x * v[i], 0));

function oklab(rgb) {
  const lin = rgb.map(srgbToLinear);
  return mul(M2, mul(M1, lin).map(Math.cbrt));
}

function oklch(rgb) {
  const [L, a, b] = oklab(rgb);
  return { L, C: Math.hypot(a, b), h: ((Math.atan2(b, a) * 180) / Math.PI + 360) % 360 };
}

const CVD = {
  protanopia: [
    [0.152286, 1.052583, -0.204868],
    [0.114503, 0.786281, 0.099216],
    [-0.003882, -0.048116, 1.051998],
  ],
  deuteranopia: [
    [0.367322, 0.860646, -0.227968],
    [0.280085, 0.672501, 0.047413],
    [-0.01182, 0.04294, 0.968881],
  ],
  tritanopia: [
    [1.255528, -0.076749, -0.178779],
    [-0.078411, 0.930809, 0.147602],
    [0.004733, 0.691367, 0.3039],
  ],
};

const simulate = (rgb, kind) => mul(CVD[kind], rgb.map(srgbToLinear)).map(linearToSrgb);
const dist = (a, b) => {
  const [x, y, z] = oklab(a);
  const [p, q, r] = oklab(b);
  return Math.hypot(x - p, y - q, z - r);
};
const hueGap = (a, b) => {
  const d = Math.abs(a - b) % 360;
  return Math.min(d, 360 - d);
};

function relLum(rgb) {
  const [r, g, b] = rgb.map(srgbToLinear);
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}
function contrast(fg, bg) {
  const l1 = relLum(fg);
  const l2 = relLum(bg);
  return (Math.max(l1, l2) + 0.05) / (Math.min(l1, l2) + 0.05);
}

// ---- 从 common.py 解析色值 -------------------------------------------------
const src = readFileSync(commonPath, "utf8");

// PAL 是个普通 dict：`"name": "#RRGGBB",`
const palBlock = src.match(/^PAL\s*=\s*\{([\s\S]*?)^\}/m);
if (!palBlock) throw new Error("common.py 里找不到 PAL = { ... }");
const PAL = {};
for (const m of palBlock[1].matchAll(/"(\w+)"\s*:\s*"(#[0-9A-Fa-f]{6})"/g)) {
  PAL[m[1]] = m[2].toUpperCase();
}
if (Object.keys(PAL).length === 0) throw new Error("PAL 解析出 0 个色值");

// PAL_CYCLE 是列表，元素写成 PAL["name"] 或字面量 —— 两种都要认
//
// **不能用 `\[([\s\S]*?)\]`**：列表元素本身带方括号（`PAL["blue"]`），
// 非贪婪匹配会停在**里面**那个 `]` 上，于是解析出 0 个色值。
// 锚到行尾的 `]` 才是列表的收尾。
const cycBlock = src.match(/^PAL_CYCLE\s*=\s*\[([\s\S]*?)\]\s*(?:#.*)?$/m);
if (!cycBlock) throw new Error("common.py 里找不到 PAL_CYCLE = [ ... ]");
const CYCLE = [];
for (const tok of cycBlock[1].matchAll(/PAL\[\s*"(\w+)"\s*\]|"(#[0-9A-Fa-f]{6})"/g)) {
  if (tok[1] !== undefined) {
    if (!(tok[1] in PAL)) throw new Error(`PAL_CYCLE 引用了不存在的 PAL["${tok[1]}"]`);
    CYCLE.push(PAL[tok[1]]);
  } else {
    CYCLE.push(tok[2].toUpperCase());
  }
}
if (CYCLE.length === 0) throw new Error("PAL_CYCLE 解析出 0 个色值");

// ---- 从 .mplstyle 解析分类循环色 -------------------------------------------
const styleSrc = readFileSync(stylePath, "utf8");
const cycLine = styleSrc.match(/^axes\.prop_cycle:.*$/m);
if (!cycLine) throw new Error("publication.mplstyle 里找不到 axes.prop_cycle");
const STYLE_CYCLE = [...cycLine[0].matchAll(/'#?([0-9A-Fa-f]{6})'/g)].map((m) =>
  ("#" + m[1]).toUpperCase()
);
if (STYLE_CYCLE.length === 0) throw new Error("axes.prop_cycle 解析出 0 个色值");

const WHITE = hexToRgb("#FFFFFF");
const MIN_CONTRAST = 2.0;
const MIN_CVD = 0.05;
const MIN_HUE = 15;
const MIN_CYCLE = 10;
/**
 * 对比度只查**真正会被画成实心点/线的颜色**。
 *
 * 一个定义了但没人用的色值不会让任何一张图看不清 —— 拿它判红会让门禁
 * 变成"色值表整洁度检查"。实测 `PAL` 里有 11 个色值，而全仓库脚本里
 * 只直接引用了 `PAL["primary"]` 与 `PAL["highlight"]` 两个；
 * 其余的分类色是经 `PAL_CYCLE` -> rcParams 隐式消费的。
 *
 * 所以判据范围 = `PAL_CYCLE` + 语义别名。未使用的条目另行**列出但不判失败**
 * （可见即可，避免"删了才发现别处在用"）。
 */
const CHECKED = new Set([...CYCLE, PAL["primary"], PAL["highlight"]].filter(Boolean));

const problems = [];

// ---- 判据 5：两处定义必须一致 ----------------------------------------------
if (CYCLE.join(",") !== STYLE_CYCLE.join(",")) {
  problems.push(
    `PAL_CYCLE 与 .mplstyle 的 axes.prop_cycle 不一致 ——\n` +
      `      common.py : ${CYCLE.join(" ")}\n` +
      `      .mplstyle : ${STYLE_CYCLE.join(" ")}\n` +
      `      两处各定义了一份，改一处就会让另一处变成谎话`
  );
}

// ---- 判据 4：色数够不够 ----------------------------------------------------
if (CYCLE.length < MIN_CYCLE) {
  problems.push(
    `PAL_CYCLE 只有 ${CYCLE.length} 色 < ${MIN_CYCLE}：簇数超过它就会有簇撞色` +
      `（实测 resolution=1.0 给 10 个簇，0≡5 / 1≡6 / …，图上完全分不出来）`
  );
}

// ---- 判据 1 + 2：分类循环色两两检查 ----------------------------------------
const cycNames = CYCLE.map((hx, i) => `cycle ${i + 1} ${hx}`);
for (let i = 0; i < CYCLE.length; i++) {
  for (let j = i + 1; j < CYCLE.length; j++) {
    const a = CYCLE[i];
    const b = CYCLE[j];
    const ca = oklch(hexToRgb(a));
    const cb = oklch(hexToRgb(b));
    if (ca.C >= 0.04 && cb.C >= 0.04) {
      const g = hueGap(ca.h, cb.h);
      if (g < MIN_HUE) {
        problems.push(`色相 ${g.toFixed(1)}° < ${MIN_HUE}°：${cycNames[i]} 与 ${cycNames[j]} 是同一颜色`);
      }
    }
    for (const kind of Object.keys(CVD)) {
      const d = dist(simulate(hexToRgb(a), kind), simulate(hexToRgb(b), kind));
      if (d < MIN_CVD) {
        problems.push(`${kind} 下 Δ=${d.toFixed(3)} < ${MIN_CVD}：${cycNames[i]} 与 ${cycNames[j]} 无法分辨`);
      }
    }
  }
}

// ---- 判据 3：白底对比度（只查会被画出来的颜色）-----------------------------
const rows = [];
const unused = [];
for (const [n, hx] of Object.entries(PAL)) {
  const c = contrast(hexToRgb(hx), WHITE);
  const checked = CHECKED.has(hx);
  rows.push(`  ${n.padEnd(12)} ${hx}  ${c.toFixed(2)}:1${checked ? "" : "   （未使用，不判）"}`);
  if (!checked) unused.push(`${n} ${hx}`);
  if (checked && c < MIN_CONTRAST) {
    problems.push(`白底对比度 ${c.toFixed(2)}:1 < ${MIN_CONTRAST}：${n} ${hx} 作为实心圆点看不清`);
  }
}

// ---- 报告 ------------------------------------------------------------------
console.log(`检查 ${commonPath}`);
console.log(`PAL ${Object.keys(PAL).length} 个色值 / PAL_CYCLE ${CYCLE.length} 色 / 样式文件 ${STYLE_CYCLE.length} 色`);
console.log(`\n白底对比度：\n${rows.join("\n")}`);
console.log(`\n分类循环色（按 draw 顺序）：`);
CYCLE.forEach((hx, i) => {
  const { L, C, h } = oklch(hexToRgb(hx));
  console.log(`  ${String(i + 1).padStart(2)}. ${hx}  L=${L.toFixed(3)} C=${C.toFixed(3)} h=${h.toFixed(1)}°`);
});
console.log(
  `\n两两检查：${CYCLE.length} 色，${(CYCLE.length * (CYCLE.length - 1)) / 2} 对，` +
    `判据 色相>=${MIN_HUE}° / 色盲 Δ>=${MIN_CVD} / 对比度>=${MIN_CONTRAST}`
);
if (unused.length) {
  console.log(
    `\n定义了但没被画出来的色值 ${unused.length} 个（不判失败，仅供清理时参考）：\n  ` +
      unused.join("\n  ")
  );
}

if (problems.length) {
  console.log(`\n未通过 ${problems.length} 项：`);
  for (const p of problems) console.log("  ✗ " + p);
  process.exit(1);
}
console.log("\n调色板检查通过。");
