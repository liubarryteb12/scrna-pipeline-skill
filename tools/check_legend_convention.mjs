#!/usr/bin/env node
/**
 * tools/check_legend_convention.mjs — 图例规范门禁（用户约定 v2）
 *
 * 规范（三仓 AGENTS 同一条：geo 31 / scrna 26 / spatial 26）：
 *
 *   **图例一律画在图框外、右侧、纵向排列。**
 *
 *     ① 框内图例压住数据点 —— 读者分不清哪块是数据、哪块是说明；
 *     ② 顶部图例会**把主图压扁变形**（实测校准图 / KM / UMAP 被压扁）。
 *
 * 三仓的 `theme_paper()` / `.mplstyle` 已经设好默认值，所以**源码里
 * 大部分调用是对的** —— 本检查器抓的是"有人就地覆盖了默认值"：
 *
 *   R (ggplot2)         不合格：legend.position = "top" / "bottom" / c(x, y)
 *                                （c(x,y) 是框内坐标；"left" 也不合右侧约定）
 *   Python (matplotlib) 不合格：ax.legend(loc="outside ...")   ← 直接 ValueError
 *                                fig.legend(...) 缺 ncol=1 或缺 outside
 *
 * ── 为什么需要它 ─────────────────────────────────────────────
 *
 * `check_r_syntax.mjs` / `check_py_syntax.mjs` 都看不见图例位置：
 * 前者只查括号配平与名字拼写，后者只查未定义名字。
 * 一张图例压在数据上的图，在**所有**既有门禁眼里都是合格的
 * （文件存在、有墨迹、图名合规、配色合规、图幅合规）。
 *
 * 这正是 `15_ERROR_LEDGER.md` E-06「门禁本身有盲区」的又一例：
 * **"检查通过"只说明被检查的那件事通过了。**
 *
 * ── 判据的严格度（刻意保守）─────────────────────────────────
 *
 * 只判**确定错**的两类：
 *   - 显式写成 top/bottom/c(x,y)（约定明确禁止）
 *   - ax.legend 用了 outside（运行期必然 ValueError，实测崩过整个 job）
 *
 * **不判**"没写 legend.position" —— 那是继承 theme 默认值，是正确写法。
 * 也不判 `"none"`（无图例）与 `"left"`：前者合法，后者是风格取舍，
 * 交给人工终审。**静默豁免会让检查退化成没有检查，过度判死会让
 * 合法写法过不去** —— 两者都要避免，所以边界写在这里。
 *
 * 用法:
 *   node tools/check_legend_convention.mjs            # 扫本仓库
 *   node tools/check_legend_convention.mjs --verbose  # 列出全部调用
 */

import { readFileSync, readdirSync, existsSync } from "node:fs";
import { join, dirname, basename } from "node:path";
import { fileURLToPath } from "node:url";

const VERBOSE = process.argv.includes("--verbose");
const repoRoot = dirname(dirname(fileURLToPath(import.meta.url)));
const scriptsDir = join(repoRoot, "scripts");
const libDir = join(scriptsDir, "lib");

/**
 * 收集 src 里所有 `name(` 调用（括号配平、跨行）。
 * @returns {{line:number,text:string}[]}
 */
function calls(src, name) {
  const out = [];
  let i = 0;
  while ((i = src.indexOf(name, i)) !== -1) {
    const open = i + name.length;
    if (src[open] !== "(") { i = open; continue; }
    let depth = 0;
    let k = open;
    for (; k < src.length; k++) {
      if (src[k] === "(") depth++;
      else if (src[k] === ")") { depth--; if (depth === 0) { k++; break; } }
    }
    out.push({ line: src.slice(0, i).split("\n").length, text: src.slice(i, k) });
    i = k;
  }
  return out;
}

/** 去掉行注释（`#` 到行尾），保留字符串内的 `#`。 */
function stripComments(src) {
  return src.split("\n").map((l) => {
    let inS = null;
    for (let i = 0; i < l.length; i++) {
      const ch = l[i];
      if (inS) { if (ch === "\\") i++; else if (ch === inS) inS = null; }
      else if (ch === '"' || ch === "'") inS = ch;
      else if (ch === "#") return l.slice(0, i);
    }
    return l;
  }).join("\n");
}

const problems = [];
const notes = [];

function scanR(file) {
  const src = stripComments(readFileSync(file, "utf8"));
  src.split("\n").forEach((line, idx) => {
    const m = /legend\.position\s*=\s*([^,)\n]+)/.exec(line);
    if (!m) return;
    const val = m[1].trim();
    if (VERBOSE) notes.push(`${basename(file)}:${idx + 1}  legend.position = ${val}`);
    // 约定明确禁止：顶部 / 底部 / 框内坐标
    if (/^["'](top|bottom)["']$/.test(val) || /^c\s*\(/.test(val)) {
      problems.push(`${basename(file)}:${idx + 1}  legend.position = ${val}` +
                    `  —— 约定 v2 要求框外右侧（"right"）`);
    }
  });
}

function scanPy(file) {
  const src = readFileSync(file, "utf8");
  for (const c of calls(src, "fig.legend")) {
    const hasNcol = /ncol\s*=/.test(c.text);
    const outside = /loc\s*=\s*["']outside/.test(c.text);
    if (VERBOSE) notes.push(`${basename(file)}:${c.line}  fig.legend  ncol=${hasNcol} outside=${outside}`);
    if (!outside) {
      problems.push(`${basename(file)}:${c.line}  fig.legend 缺 loc="outside right center"` +
                    `  —— 框内图例会压住数据点`);
    } else if (!hasNcol) {
      problems.push(`${basename(file)}:${c.line}  fig.legend 缺 ncol=1` +
                    `  —— 多图例默认会横排，必须显式纵向单列`);
    }
  }
  for (const c of calls(src, "ax.legend")) {
    if (VERBOSE) notes.push(`${basename(file)}:${c.line}  ax.legend`);
    if (/loc\s*=\s*["']outside/.test(c.text)) {
      problems.push(`${basename(file)}:${c.line}  ax.legend 用了 loc="outside ..."` +
                    `  —— 只对 fig.legend 有效，运行期必然 ValueError`);
    }
  }
}

const rFiles = [];
const pyFiles = [];
if (existsSync(scriptsDir)) {
  for (const e of readdirSync(scriptsDir)) {
    if (e.endsWith(".R")) rFiles.push(join(scriptsDir, e));
    if (e.endsWith(".py")) pyFiles.push(join(scriptsDir, e));
  }
}
if (existsSync(libDir)) {
  for (const e of readdirSync(libDir)) if (e.endsWith(".R")) rFiles.push(join(libDir, e));
}

if (rFiles.length === 0 && pyFiles.length === 0) {
  console.error(`[FAIL] ${scriptsDir} 下没找到 .R / .py —— 检查没有实际执行，不能算通过`);
  process.exit(1);
}

for (const f of rFiles) scanR(f);
for (const f of pyFiles) scanPy(f);

const scanned = rFiles.length + pyFiles.length;
console.log(`图例规范检查（约定 v2：框外右侧、纵向单列）`);
console.log(`  扫描 ${scanned} 个脚本（R ${rFiles.length} / Python ${pyFiles.length}）`);

if (VERBOSE && notes.length) {
  console.log("");
  console.log("  全部调用：");
  for (const n of notes) console.log(`    ${n}`);
}

if (problems.length > 0) {
  console.error("");
  console.error(`[FAIL] ${problems.length} 处图例违反约定 v2：`);
  for (const p of problems) console.error(`    ${p}`);
  console.error("");
  console.error("  修法：");
  console.error('    R      : theme(legend.position = "right") —— 不要 top/bottom/c(x, y)');
  console.error('    Python : fig.legend(loc="outside right center", ncol=1)');
  console.error("             loc=\"outside ...\" 只对 fig.legend 有效，ax.legend 会报错。");
  process.exit(1);
}

console.log("");
console.log("  未发现违反约定 v2 的图例写法");
