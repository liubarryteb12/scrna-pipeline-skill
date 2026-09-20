#!/usr/bin/env node
/**
 * 出图命名门禁：图名必须与它所在的代码阶段对应。
 *
 *     node tools/check_fig_names.mjs
 *
 * 命名格式（五个字段，用 `-` 连起来）：
 *
 *     <阶段>-<模块>-<图>-unit<单元>-<名称>
 *      02      03     01   unit1     umap-clusters
 *
 *   * **阶段** = 三大部分之一：`01` geo / `02` scrna / `03` spatial
 *   * **模块** = 脚本序号，取自脚本文件名（`03_cluster_annotate.py` -> `03`）
 *   * **图**   = 该脚本内第几张图，两位数字，从 `01` 起
 *   * **单元** = 同一张图里的功能单元。**拆成单图后共用一个图号** ——
 *     这样"这几个文件原本是一张图"这个来源信息还在。
 *   * **名称** = 小写连字符的 slug
 *
 * **为什么要写成门禁而不是靠人记：** 图名和脚本序号是**两处**，
 * 而"图名里的模块号写错了"这件事**没有任何东西能发现** ——
 * 图照样生成、CI 照样绿、验收照样过，只是读者按图名去 `scripts/` 里找
 * 对应代码时会找错文件。这正是"两处定义必然分叉"那一类。
 *
 * 检查四条：
 *   1. 每个出图调用的名字都符合上面的格式
 *   2. **模块号与所在脚本的文件名一致**
 *   3. 同一脚本内图号从 `01` 起连续，同一图号下的 `unit` 从 `1` 起连续
 *   4. 全仓库没有重名（重名会让后写的图覆盖先写的）
 *
 * 退出码：0 = 全部合规；1 = 有不合规。
 * 纯 Node、零依赖。
 */

import { readFileSync, readdirSync } from "node:fs";
import { join, basename } from "node:path";
import { fileURLToPath } from "node:url";
import { dirname } from "node:path";

const here = dirname(fileURLToPath(import.meta.url));
const repo = join(here, "..");

/** 三大部分 -> 阶段号。仓库目录名 -> 阶段号。 */
const PART_BY_REPO = {
  "geo-normal-pipeline-skill": "01",
  "scrna-pipeline-skill": "02",
  "spatial-pipeline-skill": "03",
};
const repoName = basename(repo);
const PART = PART_BY_REPO[repoName];
if (!PART) {
  console.error(`认不出仓库 ${repoName} —— 阶段号表里没有它`);
  process.exit(1);
}

/**
 * 剥掉注释，**保持行数不变**（行号要用来报错）。
 *
 * 不剥会**误报**：注释里写了 `save_fig()` 时计数会被当成一个调用点，
 * 而它当然认不出字面量 —— 于是报"名字是拼出来的"。
 * Python 的 `#` 与 R 的 `#` 规则相同，逐行扫、跟踪引号状态：
 * `#` 只有在**引号外**才是注释开始。
 */
function stripComments(src) {
  return src
    .split(/\r?\n/)
    .map((line) => {
      let q = null;
      for (let i = 0; i < line.length; i++) {
        const c = line[i];
        if (q) {
          if (c === q) q = null;
        } else if (c === '"' || c === "'") {
          q = c;
        } else if (c === "#") {
          return line.slice(0, i);
        }
      }
      return line;
    })
    .join("\n");
}

/**
 * 收集脚本里所有**符合命名格式的字符串字面量**。
 *
 * **为什么不要求名字必须写在 `save_fig` 的实参位置上：**
 * 经辅助函数传字面量是合法写法，而字面量照样可 grep。
 * 要求"实参位置必须是字面量"会逼出把函数体复制两遍这种更差的代码。
 *
 * 代价要说清楚：**账目对齐只能保证"字面量数 >= 出图调用数"**，
 * 挡不住"一个调用用了拼出来的名字、同时另有一个没用的字面量"。
 * 所以拼出来的调用会**逐条列出**（可见），只是不判失败 ——
 * 静默豁免会让检查退化成没有检查，完全判死又会让合法写法过不去。
 */
function collectLiterals(src) {
  const out = [];
  const lineOf = (idx) => src.slice(0, idx).split(/\r?\n/).length;
  // 只认"长得像图名"的字面量：以 <阶段>- 开头
  const re = new RegExp(`"(${PART}-[^"]*)"`, "g");
  let m;
  while ((m = re.exec(src)) !== null) {
    out.push({ name: m[1].replace(/\.(pdf|png)$/i, ""), line: lineOf(m.index) });
  }
  return out;
}

/** 数出图调用点，并取出每个调用写的名字（实参不是字面量时是 null）。 */
function collectCalls(src) {
  const calls = [];
  const lineOf = (idx) => src.slice(0, idx).split(/\r?\n/).length;
  // Python 侧：`save_fig(cfg, "name"` 是唯一写法。
  // R 侧（geo）写 `save_pdf(file.path(res, "name.pdf")` —— 两个模式都留着，
  // 这样三个仓库共用同一份门禁，不会各自漂移。
  const patterns = [
    { re: /save_fig\s*\(/g, lit: /save_fig\(\s*cfg\s*,\s*"([^"]+)"/g },
    { re: /save_pdf\s*\(/g, lit: /save_pdf\(\s*file\.path\([^,]+,\s*"([^"]+)"\s*\)/g },
  ];
  const byLine = new Map();
  for (const { lit } of patterns) {
    let m;
    while ((m = lit.exec(src)) !== null) {
      byLine.set(lineOf(m.index), m[1].replace(/\.(pdf|png)$/i, ""));
    }
  }
  for (const { re } of patterns) {
    let m;
    while ((m = re.exec(src)) !== null) {
      const ln = lineOf(m.index);
      calls.push({ line: ln, name: byLine.get(ln) || null });
    }
  }
  return calls;
}

const NAME_RE = new RegExp(`^${PART}-(\\d{2})-(\\d{2})-unit(\\d+)-[a-z0-9]+(?:-[a-z0-9]+)*$`);

const problems = [];
const notes = [];
const seen = new Map();
let nCallsTotal = 0;

const scriptsDir = join(repo, "scripts");
const scripts = readdirSync(scriptsDir)
  .filter((f) => /^\d+_.*\.(R|py)$/.test(f))
  .sort();

for (const f of scripts) {
  const moduleNo = f.slice(0, 2);
  const src = stripComments(readFileSync(join(scriptsDir, f), "utf8"));
  const lits = collectLiterals(src);
  const calls = collectCalls(src);
  nCallsTotal += calls.length;
  if (calls.length === 0) continue;

  // **同一个名字在脚本里出现多次是正常的** —— `figs_written` /
  // `figs_missing` 这类记账列表会再引一遍。所以按名字去重。
  const figs = [];
  const byName = new Map();
  for (const item of lits) {
    const m = NAME_RE.exec(item.name);
    if (!m) {
      problems.push(
        `${f}:${item.line} 名字 "${item.name}" 不符合 <阶段>-<模块>-<图>-unit<单元>-<名称>，` +
          `应为 ${PART}-${moduleNo}-NN-unitN-<slug>`
      );
      continue;
    }
    if (m[1] !== moduleNo) {
      problems.push(
        `${f}:${item.line} 名字里的模块号是 ${m[1]}，但它在脚本 ${f}（模块 ${moduleNo}）里 —— ` +
          `读者会去 scripts/${m[1]}_* 找代码`
      );
    }
    if (m[2] === "00") problems.push(`${f}:${item.line} 图号不能是 00（从 01 起）`);
    if (byName.has(item.name)) continue;
    byName.set(item.name, item.line);
    figs.push({ fig: Number(m[2]), unit: Number(m[3]), name: item.name, line: item.line });

    // **跨脚本重名才算问题**（同一脚本内重复是记账引用）
    if (seen.has(item.name)) {
      problems.push(
        `${f}:${item.line} 名字 "${item.name}" 与 ${seen.get(item.name)} 重复 —— 两个脚本写同一个文件`
      );
    } else {
      seen.set(item.name, `${f}:${item.line}`);
    }
  }

  // **同一个脚本里两处 save 调用写同一个名字** = 后写的覆盖先写的。
  // 这个只能从**调用点**看：字面量去重之后它已经看不出来了。
  const callNames = calls.map((c) => c.name).filter(Boolean);
  const dupInScript = callNames.filter((n, i) => callNames.indexOf(n) !== i);
  for (const n of new Set(dupInScript)) {
    problems.push(`${f}: 有两处 save 调用都写 "${n}" —— 后写的会覆盖先写的`);
  }

  // **账目对齐**：字面量数不能少于出图调用数。
  // 少了就说明有调用用了拼出来的名字，而那个名字没被任何检查看过。
  if (lits.length < calls.length) {
    problems.push(
      `${f}: ${calls.length} 处出图调用，但只有 ${lits.length} 个合规图名字面量 —— ` +
        `有调用用的是拼出来的名字，静态检查看不见它`
    );
  }
  for (const c of calls) {
    if (!c.name) {
      notes.push(`${f}:${c.line} 出图调用的实参不是字面量（经辅助函数传名，字面量在别处）`);
    }
  }

  // 图号从 01 起连续
  const figNos = [...new Set(figs.map((x) => x.fig))].sort((a, b) => a - b);
  figNos.forEach((n, i) => {
    if (n !== i + 1) {
      problems.push(`${f}: 图号 ${String(n).padStart(2, "0")} 不连续（应为 ${String(i + 1).padStart(2, "0")}）`);
    }
  });
  // 同一图号下 unit 从 1 起连续
  for (const n of figNos) {
    const units = figs.filter((x) => x.fig === n).map((x) => x.unit).sort((a, b) => a - b);
    units.forEach((u, i) => {
      if (u !== i + 1) {
        problems.push(`${f}: 图 ${String(n).padStart(2, "0")} 的 unit 编号 ${u} 不连续（应为 ${i + 1}）`);
      }
    });
  }
}

console.log(`检查 ${repoName}（阶段 ${PART}）`);
console.log(`  ${scripts.length} 个脚本，${nCallsTotal} 处出图调用，${seen.size} 个唯一图名`);
const multi = [...seen.keys()].filter((n) => /-unit[2-9]/.test(n));
console.log(`  多单元图：${multi.length ? multi.join(", ") : "（无）"}`);
if (notes.length) {
  console.log(`\n以下调用的实参不是字面量（可见，不判失败）：`);
  for (const n of notes) console.log("  · " + n);
}

if (problems.length) {
  console.error(`\n未通过 ${problems.length} 项：`);
  for (const p of problems) console.error("  ✗ " + p);
  process.exit(1);
}
console.log("\n出图命名检查通过。");
