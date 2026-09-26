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
 * ── E-62：Python 侧曾把注释当代码（已修，附 `--selftest`）────────
 *
 * 修复前 `scanPy()` 直接把**原始源码**喂给 `calls()`，而 `calls()` 用
 * `src.indexOf(name)` 找调用点、靠**括号配平**取实参文本。注释里出现的
 * `fig.legend(` 因此有两个方向的错：
 *
 *   ① **假阳性**：注释里的 `fig.legend(...)` 被当成真调用去判 `loc`；
 *   ② **假阴性**：注释里那个括号**不配平**时，`calls()` 一路吞到文件尾，
 *      它后面**所有真调用永不被检查** —— 门禁照样打绿。
 *
 * 实测最小复现（`--selftest` 用例 F）：注释里同时出现 `fig.legend(` 与
 * 字样 `loc="outside right center"`，紧随其后的真调用缺 `loc` —— **放行**。
 * 这不是假想：本仓 `03_spatial_domains.py` 的注释在说明这个检查器时
 * 恰好写出了 `fig.legend(`，于是它自己的真调用（L134）从未被检查过。
 *
 * 修法：`blankNonCode()` 先把**注释与三引号字符串**逐字符换成空格
 * （保留换行，故行号与偏移都不变），再交给 `calls()`。单引号/双引号
 * 字符串**保留内容**，因为 `loc="outside ..."` 正是靠它判的。
 *
 * 教训（与 E-61 同一条）：**检查器自己也要被检查** —— 所以本文件带
 * `--selftest`，把"故意塞一个错"的用例写进代码里，不靠临时脚本。
 *
 * ── 已知盲区（不修，写在这里）───────────────────────────────
 *
 * `fig.legend(**legend_kw)` 这类**经 kwargs 转发**的调用，纯文本检查器
 * 看不见 `loc`/`ncol`。跨函数数据流分析成本高且会引入新的假阴性，
 * 所以约定 v2 的两个关键参数应当在 helper 里**焊死**（`fig.legend(
 * loc="outside right center", ncol=1, **legend_kw)`），把可变部分留给
 * 调用方 —— 这样检查器看得见，调用方也无法静默覆盖约定。
 *
 * 用法:
 *   node tools/check_legend_convention.mjs            # 扫本仓库
 *   node tools/check_legend_convention.mjs --verbose  # 列出全部调用
 *   node tools/check_legend_convention.mjs --selftest # 只跑内建标定用例
 */

import { readFileSync, readdirSync, existsSync } from "node:fs";
import { join, dirname, basename } from "node:path";
import { fileURLToPath } from "node:url";

const VERBOSE = process.argv.includes("--verbose");
const SELFTEST = process.argv.includes("--selftest");
const repoRoot = dirname(dirname(fileURLToPath(import.meta.url)));
const scriptsDir = join(repoRoot, "scripts");
const libDir = join(scriptsDir, "lib");

/**
 * 把**注释**与**三引号字符串**的字符换成空格（保留 `\n`，故偏移/行号不变）。
 *
 * 单引号/双引号字符串**原样保留** —— `loc="outside right center"` 的内容
 * 是判据本身，抹掉就没法判了。
 *
 * 为什么要抹注释：`calls()` 靠 `indexOf` + 括号配平取实参，注释里的
 * `fig.legend(` 会被当调用；括号不配平时还会一路吞掉后面的真调用（E-62）。
 * 为什么要抹三引号字符串：文档字符串里成片的代码示例不是被执行的调用。
 *
 * @param {string} src
 * @returns {string} 与 src **等长**的字符串
 */
function blankNonCode(src) {
  const a = src.split("");
  const n = src.length;
  const blank = (from, to) => {
    for (let j = from; j < to; j++) if (a[j] !== "\n") a[j] = " ";
  };
  let i = 0;
  while (i < n) {
    const ch = src[i];
    if (ch === "#") {
      let j = i;
      while (j < n && src[j] !== "\n") j++;
      blank(i, j);
      i = j;
      continue;
    }
    if (ch === '"' || ch === "'") {
      const t = src.substr(i, 3);
      if (t === '"""' || t === "'''") {
        let j = i + 3;
        while (j < n && src.substr(j, 3) !== t) {
          if (src[j] === "\\") j++;
          j++;
        }
        j = Math.min(n, j + 3);
        blank(i, j);
        i = j;
        continue;
      }
      // 单/双引号字符串：跳过但**不抹**，内容要留给判据。
      let j = i + 1;
      while (j < n && src[j] !== ch && src[j] !== "\n") {
        if (src[j] === "\\") j++;
        j++;
      }
      i = Math.min(n, j + 1);
      continue;
    }
    i++;
  }
  return a.join("");
}

/**
 * 收集 src 里所有 `name(` 调用（括号配平、跨行）。
 *
 * **必须传 `blankNonCode()` 处理过的源码** —— 否则注释里的调用点会污染
 * 结果（E-62）。行号与实参文本仍按原样返回，因为抹白是等长的。
 *
 * @param {string} src 已抹掉注释与三引号字符串的源码
 * @param {string} name 如 `"fig.legend"`
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

/**
 * 扫一段 Python 源码，返回它的问题与备注（不写全局，方便 `--selftest`）。
 * @param {string} src
 * @param {string} label 报错前缀（通常是文件名）
 */
function scanPySource(src, label) {
  const code = blankNonCode(src);
  const probs = [];
  const nts = [];
  for (const c of calls(code, "fig.legend")) {
    const hasNcol = /ncol\s*=/.test(c.text);
    const outside = /loc\s*=\s*["']outside/.test(c.text);
    if (VERBOSE) nts.push(`${label}:${c.line}  fig.legend  ncol=${hasNcol} outside=${outside}`);
    if (!outside) {
      probs.push(`${label}:${c.line}  fig.legend 缺 loc="outside right center"` +
                 `  —— 框内图例会压住数据点`);
    } else if (!hasNcol) {
      probs.push(`${label}:${c.line}  fig.legend 缺 ncol=1` +
                 `  —— 多图例默认会横排，必须显式纵向单列`);
    }
  }
  for (const c of calls(code, "ax.legend")) {
    if (VERBOSE) nts.push(`${label}:${c.line}  ax.legend`);
    if (/loc\s*=\s*["']outside/.test(c.text)) {
      probs.push(`${label}:${c.line}  ax.legend 用了 loc="outside ..."` +
                 `  —— 只对 fig.legend 有效，运行期必然 ValueError`);
    }
  }
  return { probs, nts };
}

function scanPy(file) {
  const r = scanPySource(readFileSync(file, "utf8"), basename(file));
  problems.push(...r.probs);
  notes.push(...r.nts);
}

// ── 内建标定（E-62 教训：检查器自己也要被检查）────────────────────
//
// 每个用例 = 一段源码 + 期望的问题条数。**必须包含反向用例**：
// 只跑"干净源码不报错"看不出假阴性，只跑"脏源码报错"看不出假阳性。
if (SELFTEST) {
  const OK = 'fig.legend(loc="outside right center", ncol=1)';
  const CASES = [
    ["基线：正确写法", [OK], 0],
    ["真缺陷：缺 loc", ["fig.legend(ncol=1)"], 1],
    ["真缺陷：缺 ncol", ['fig.legend(loc="outside right center")'], 1],
    ["真缺陷：ax.legend 用 outside", ['ax.legend(loc="outside right center")'], 1],
    ["假阳性回归：注释里有违规写法", [`# ${'fig.legend(loc="top")'}`, OK], 0],
    ["假阳性回归：注释里有 fig.legend( 字样", ['# 见 `fig.legend(` 的说明', OK], 0],
    ["假阳性回归：docstring 里有违规写法",
     ['def f():\n    """\n    ' + 'fig.legend(loc="top")\n    """\n', OK], 0],
    ["假阴性回归：注释里的 `fig.legend(` 不得吞掉后面的真缺陷",
     ['# 见 `fig.legend(` 的说明', "fig.legend(ncol=1)"], 1],
    ["假阴性回归：注释里带 loc 字样也不得吞掉后面的真缺陷",
     ['# 判据是 `fig.legend(` 里的 loc="outside right center" 与 ncol=1',
      "fig.legend(ncol=1)"], 1],
    ["假阴性回归：注释里的括号不配平",
     ["# 早期写法 fig.legend( 已废弃", OK], 0],
  ];
  let bad = 0;
  console.log("check_legend_convention.mjs 自检（内建标定用例）");
  for (const [name, lines, want] of CASES) {
    const { probs } = scanPySource(lines.join("\n") + "\n", "case.py");
    const ok = probs.length === want;
    if (!ok) bad++;
    console.log(`  ${ok ? "[ok]" : "[FAIL]"} ${name}  期望 ${want} 条，实际 ${probs.length} 条`);
    if (!ok) probs.forEach((p) => console.log(`         ${p}`));
  }
  console.log("");
  if (bad > 0) {
    console.error(`[FAIL] 自检 ${bad}/${CASES.length} 个用例不通过 —— 检查器本身有问题。`);
    process.exit(1);
  }
  console.log(`  自检通过（${CASES.length} 个用例，含 4 条假阳性/假阴性回归）`);
  process.exit(0);
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
