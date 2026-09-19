#!/usr/bin/env node
/**
 * tools/check_artifact_paths.mjs — 交叉检查"脚本写出的 data/ 文件"与
 * "workflow 上传的 artifact 清单"
 *
 * **这个工具的存在理由（Part 1 踩过的真坑）。**
 * 验收清单检查的是 **runner 工作目录**里的文件，artifact 的 `path:`
 * 是一个**完全独立的列表**。两者不一致时：
 *   - 验收通过（文件确实生成了）
 *   - job 绿（没报错）
 *   - 但用户下载 artifact 后发现核心产物不在里面
 * 实测漏了 5 个文件，其中包括 44 MB 的核心表达矩阵。
 *
 * 做法：
 *   1. 从 workflow 里解析 artifact 的 `path:` 块
 *   2. 从 scripts/**.py 里抓 `data_dir / "X"` 与 `Path(...data_dir...) / "X"` 字面量
 *   3. 差集 = 写了但没上传的文件
 *
 * 只检查 data/（大文件、二进制、真正会被漏的那一类）。
 * results/ 整个目录上传，不需要逐文件核对。
 *
 * 用法: node tools/check_artifact_paths.mjs
 * 退出码: 0 = 一致；1 = 有文件写了但不在 artifact 清单里
 */

import { readdirSync, readFileSync, existsSync } from "node:fs";
import { join, relative } from "node:path";

const REPO = new URL("..", import.meta.url).pathname.replace(/^\/([A-Za-z]:)/, "$1");
const WF = join(REPO, ".github", "workflows", "scrna_analysis.yml");

if (!existsSync(WF)) {
  console.error(`找不到 workflow: ${WF}`);
  process.exit(1);
}

// ---- 1. 解析 artifact 的 path: 块 -------------------------------------------
const wfLines = readFileSync(WF, "utf8").split(/\r?\n/);
const pathIdx = wfLines.findIndex((l) => /^\s*path:\s*\|/.test(l));
if (pathIdx < 0) {
  console.error("未能从 workflow 里解析出 artifact path 清单（找不到 `path: |`）");
  process.exit(1);
}
const pathIndent = wfLines[pathIdx].match(/^(\s*)/)[1].length;
const artifactPaths = [];
for (let i = pathIdx + 1; i < wfLines.length; i++) {
  const line = wfLines[i];
  if (line.trim() === "") continue;
  const m = line.match(/^(\s*)(\S.*)$/);
  if (!m) continue;
  // 缩进回到 path: 那一层或更浅 -> 块结束
  if (m[1].length <= pathIndent) break;
  artifactPaths.push(m[2].trim());
}
if (artifactPaths.length === 0) {
  console.error("artifact path 清单解析为空 —— workflow 格式变了？");
  process.exit(1);
}

// ---- 2. 从 scripts 里抓 data_dir 下的文件名字面量 ----------------------------
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

const pyFiles = walk(join(REPO, "scripts"));
const written = new Map();   // 文件名 -> 首次出现的 文件:行
// 匹配:  data_dir / "X"
const RE = /data_dir\s*\/\s*"([^"]+)"/g;

for (const f of pyFiles) {
  const rel = relative(REPO, f).replace(/\\/g, "/");
  readFileSync(f, "utf8").split(/\r?\n/).forEach((line, i) => {
    for (const m of line.matchAll(RE)) {
      const name = m[1];
      if (name.includes("/")) continue;          // 只关心 data/ 顶层文件
      // **跳过目录。** `data_dir / "cache"` 是解压/下载缓存目录，
      // 不是产物。用"有没有扩展名"判断，够用且不会误伤。
      if (!name.includes(".")) continue;
      if (!written.has(name)) written.set(name, `${rel}:${i + 1}`);
    }
  });
}

if (written.size === 0) {
  console.warn("警告：没有从 scripts 里抓到任何 data_dir 下的文件名 —— 正则可能过时了");
}

// ---- 3. 求差集 --------------------------------------------------------------
// artifact 路径形如 `data/<dataset_id>/raw.h5ad`，脚本里是
// `data_dir / "raw.h5ad"` —— 两边都归一到 **basename** 再比。
// （第一版只剥了 `data/` 前缀，没处理数据集 id 那一段，全是误报。）
const inArtifact = new Set(
  artifactPaths.filter((p) => p.startsWith("data/"))
    .map((p) => p.split("/").pop()));

const missing = [...written.entries()]
  .filter(([name]) => !inArtifact.has(name))
  .sort();

console.log(`artifact 清单: ${artifactPaths.length} 项，其中 data/ 下 ${inArtifact.size} 项`);
console.log(`脚本写出的 data/ 文件: ${written.size} 个`);

// 反向：artifact 里列了但脚本从不写（拼错文件名 / 删了代码没删清单）
const neverWritten = [...inArtifact].filter((n) => !written.has(n)).sort();

let bad = 0;
if (missing.length > 0) {
  console.error(`\n以下 ${missing.length} 个文件会被脚本写出，但**不在 artifact 清单里**`);
  console.error("（验收会通过、job 会绿，但用户下载 artifact 时拿不到）:");
  for (const [name, where] of missing) {
    console.error(`  - data/${name}   (写在 ${where})`);
  }
  bad++;
}
if (neverWritten.length > 0) {
  console.error(`\nartifact 清单里有 ${neverWritten.length} 项脚本从不写出（可能拼错或已废弃）:`);
  for (const n of neverWritten) console.error(`  - data/${n}`);
  bad++;
}

if (bad > 0) {
  console.error("\n请让两边一致后再提交");
  process.exit(1);
}
console.log("data/ 产物与 artifact 清单一致");
