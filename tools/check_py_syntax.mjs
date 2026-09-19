#!/usr/bin/env node
/**
 * tools/check_py_syntax.mjs — 所有 Python 脚本的语法检查
 *
 * **为什么需要它。** CI 里跑真实流程要装 scanpy 全家桶、下载数据、
 * 跑几分钟。语法错误在第 1 秒就能发现，却要等 5 分钟才暴露。
 * 这一步在依赖安装之前跑，把"手滑打错字"这类问题挡在最前面。
 *
 * 用 `python -m py_compile`（真正的编译，不是正则匹配）。它会写出
 * __pycache__，检查完清掉。
 *
 * 另外做三条本仓库特有的检查（见 AGENTS.md 规则）：
 *   1. 不允许在 scripts/ 里用 print() 直接输出日志 —— 必须走 common 的
 *      log_info/log_warn/log_error（否则 CI 日志没有时间戳和级别）
 *   2. 不允许硬编码 results/ 或 data/ 路径 —— 必须从 cfg 派生
 *   3. **用到 common 里的符号就必须 import 它。**
 *
 * 第 3 条是补的：`py_compile` 只做编译，**编译期看不出未定义名字** ——
 * 漏 import 一个 `W_SINGLE` 时它照样报"语法通过"，要等运行时才炸。
 * 实测就是这么漏的，白跑了一轮流水线。
 *
 * 用法: node tools/check_py_syntax.mjs
 */

import { execFileSync } from "node:child_process";
import { readdirSync, readFileSync, rmSync, existsSync, statSync } from "node:fs";
import { join, relative } from "node:path";

const REPO = new URL("..", import.meta.url).pathname.replace(/^\/([A-Za-z]:)/, "$1");
const SCAN_DIRS = ["scripts", "tools"];

function walk(dir, out = []) {
  if (!existsSync(dir)) return out;
  for (const e of readdirSync(dir, { withFileTypes: true })) {
    const p = join(dir, e.name);
    if (e.isDirectory()) {
      if (e.name === "__pycache__" || e.name.startsWith(".")) continue;
      walk(p, out);
    } else if (e.name.endsWith(".py")) {
      out.push(p);
    }
  }
  return out;
}

const files = SCAN_DIRS.flatMap((d) => walk(join(REPO, d)));
if (files.length === 0) {
  console.error("没有找到任何 .py 文件 —— 目录结构不对？");
  process.exit(1);
}

let failed = 0;

// ---- 1. 真正的语法编译 -----------------------------------------------------
console.log(`编译检查 ${files.length} 个 Python 文件`);
try {
  execFileSync("python", ["-m", "py_compile", ...files], { stdio: "pipe" });
} catch (e) {
  const out = `${e.stdout ?? ""}${e.stderr ?? ""}`;
  console.error("语法检查失败:\n" + out);
  failed++;
}

// 清掉 py_compile 产生的 __pycache__
for (const d of SCAN_DIRS) {
  const pc = join(REPO, d, "__pycache__");
  if (existsSync(pc)) rmSync(pc, { recursive: true, force: true });
}

// ---- 2. 本仓库特有的规则 ---------------------------------------------------
const RULES = [
  {
    name: "scripts/ 里不能用裸 print() 输出日志",
    // 允许 print 出现在 __main__ 的极少数地方？不允许 —— 一律走 common。
    test: (line) => /^\s*print\s*\(/.test(line),
    dirs: ["scripts"],
    hint: "改用 common 的 log_info / log_warn / log_error",
    // common.py 自己实现日志，内部用 print
    exempt: ["scripts/lib/common.py"],
  },
  {
    name: "不能硬编码 results/ 或 data/ 路径",
    test: (line) => /["'`](results|data)\//.test(line),
    dirs: ["scripts"],
    hint: "从 cfg['output']['results_dir'|'data_dir'|'figures_dir'] 派生",
    exempt: ["scripts/lib/common.py"],
  },
];

for (const f of files) {
  const rel = relative(REPO, f).replace(/\\/g, "/");
  if (!RULES.some((r) => r.dirs.some((d) => rel.startsWith(d + "/")))) continue;
  const lines = readFileSync(f, "utf8").split(/\r?\n/);
  for (const rule of RULES) {
    if (!rule.dirs.some((d) => rel.startsWith(d + "/"))) continue;
    if (rule.exempt.includes(rel)) continue;
    lines.forEach((line, i) => {
      if (rule.test(line)) {
        console.error(`  ${rel}:${i + 1}  ${rule.name}`);
        console.error(`      ${line.trim()}`);
        console.error(`      -> ${rule.hint}`);
        failed++;
      }
    });
  }
}

// ---- 3. 用到这几个 common 符号就必须 import -------------------------------
//
// 为什么是**一份写死的名单**而不是"common 导出的所有名字"：
// 后者会把函数参数名当成用法 —— `alignment.py` 里
// `def verify_alignment(adata, log_info=None, ...)` 的 log_info 是参数，
// 不是漏 import，第一版就误报了。要正确处理得做作用域分析，
// 那是重写一个 linter。
//
// 这份名单是**实际会漏的那一组**：新加的样式/几何符号。漏了它们在
// 运行时才 NameError（`py_compile` 看不出未定义名字），而 CI 要跑几分钟
// 才会撞上 —— 实测已经因此白跑过一轮。
const WATCH_SYMBOLS = [
  "W_SINGLE", "W_ONE_HALF", "W_DOUBLE", "mm",
  "PAL", "PAL_CYCLE", "apply_style",
];

/** 去掉注释与字符串字面量，避免"名字只出现在注释里"的误报。 */
function stripComments(src) {
  return src
    .replace(/"""[\s\S]*?"""/g, '""')
    .replace(/'''[\s\S]*?'''/g, "''")
    .replace(/#[^\n]*/g, "")
    .replace(/"(?:[^"\\\n]|\\.)*"/g, '""')
    .replace(/'(?:[^'\\\n]|\\.)*'/g, "''");
}

for (const f of files) {
  const rel = relative(REPO, f).replace(/\\/g, "/");
  if (rel === "scripts/lib/common.py") continue;
  const src = readFileSync(f, "utf8");

  // 该文件从 common import 了哪些名字（支持多行括号块）
  // 用 /#[^\n]*/ 而不是 /#.*$/ —— `.` 不匹配换行，用 $ 会把注释后面
  // 所有名字一起吃掉，产生"明明 import 了却报缺失"的误报。
  const imported = new Set();
  for (const m of src.matchAll(/^from\s+common\s+import\s*\(([\s\S]*?)\)/gm)) {
    for (const n of m[1].split(",")) {
      const name = n.replace(/#[^\n]*/, "").trim().split(/\s+as\s+/).pop().trim();
      if (name) imported.add(name);
    }
  }
  for (const m of src.matchAll(/^from\s+common\s+import\s+([^(\n]+)$/gm)) {
    for (const n of m[1].split(",")) {
      const name = n.replace(/#[^\n]*/, "").trim().split(/\s+as\s+/).pop().trim();
      if (name) imported.add(name);
    }
  }

  // 剥掉注释和字符串再找用法，并把 import 块本身也去掉
  const body = stripComments(src).replace(/^from\s+common\s+import[\s\S]*?\)\s*$/gm, "");
  const missing = WATCH_SYMBOLS.filter(
    (name) => !imported.has(name) && new RegExp(`(?<![\\w.])${name}\\b`).test(body)
  );
  if (missing.length) {
    console.error(`  ${rel}  用到了 common 的符号但没 import:`);
    console.error(`      ${missing.join(", ")}`);
    console.error(`      -> 加进 'from common import (...)'，否则运行时才 NameError`);
    failed++;
  }
}

if (failed > 0) {
  console.error(`\n${failed} 项检查失败`);
  process.exit(1);
}
console.log("语法与规则检查全部通过");
