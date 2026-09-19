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
 * 另外做两条本仓库特有的检查（见 AGENTS.md 规则）：
 *   1. 不允许在 scripts/ 里用 print() 直接输出日志 —— 必须走 common 的
 *      log_info/log_warn/log_error（否则 CI 日志没有时间戳和级别）
 *   2. 不允许硬编码 results/ 或 data/ 路径 —— 必须从 cfg 派生
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

if (failed > 0) {
  console.error(`\n${failed} 项检查失败`);
  process.exit(1);
}
console.log("语法与规则检查全部通过");
