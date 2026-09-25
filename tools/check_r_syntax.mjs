#!/usr/bin/env node
/**
 * check_r_syntax.mjs — 无 R 运行时的 R 静态检查（scrna 版）
 *
 * 移植自 `geo-normal-pipeline-skill/tools/check_r_syntax.mjs`，去掉了 geo 专属的
 * 检查（run_XX 编排函数、cfg$ 字段、ggplot 副标题、common.R 尺寸符号 —— 本仓库
 * 的 R 侧只有一个独立脚本，没有那些结构），保留了通用的四条：
 *
 *   1. 括号 / 引号配平
 *   2. 相邻字符串字面量没有逗号（R 没有隐式拼接，是语法错误）
 *   3. sprintf 格式串里的裸 %
 *   4. stats:: 误用基础绘图函数（本仓库 R 侧不画图，这条实际不会触发，
 *      但留着成本为零，且将来加 R 出图时立刻有用）
 *
 * **本仓库新增第五条，也是最有价值的一条**：R 脚本的 CLI 契约与 Python 调用点
 * 必须一致 —— R 里 `DEFAULTS` 的参数名、Python 里 `--xxx` 的拼写、
 * R 写出的文件名、Python 读回的文件名。这四样东西分布在两个语言里，
 * **没有任何编译器能发现它们不一致**：
 *   - Python 传了 R 不认识的参数 → R 侧 `stop("未知参数")`，但那是运行时；
 *   - R 写了 `a.csv`、Python 读 `b.csv` → Python 报"退出码 0 但没写出文件"，
 *     看起来像 R 崩了，其实只是名字对不上；
 *   - R 参数名拼错（`--n_nets` vs `--n_net`）→ 同样到运行时才炸。
 * 这正是"两处定义必然分叉"那一类，必须在推送前静态挡掉。
 *
 * 用法: node tools/check_r_syntax.mjs [目录，默认 scripts]
 */

import { readFileSync, readdirSync, statSync, existsSync } from 'node:fs'
import { join, relative, resolve, dirname } from 'node:path'
import { fileURLToPath } from 'node:url'

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..')
const targetDir = resolve(root, process.argv[2] ?? 'scripts')

const problems = []
const files = []

function walk(dir) {
  for (const entry of readdirSync(dir)) {
    const full = join(dir, entry)
    if (statSync(full).isDirectory()) walk(full)
    else if (entry.endsWith('.R')) files.push(full)
  }
}
walk(targetDir)
files.sort()

if (files.length === 0) {
  // **一个 R 文件都没找到时必须报错，不能静默通过。**
  // 目录改名 / 脚本被挪走时，整段检查会静默消失 —— 门禁看着还在、
  // 其实已经不检查任何东西了。这类"门禁自己失效"比门禁报错危险得多。
  console.log('发现 1 个问题:')
  console.log('  - scripts/ 下没有找到任何 .R 文件 —— R 静态检查无法执行')
  process.exit(1)
}

/**
 * 去掉注释与字符串字面量，只留下结构性字符。
 * 同时报告未闭合的字符串（R 里这会导致后续整段被吞掉，报错位置离真因很远）。
 */
function stripLiterals(source, file) {
  let out = ''
  let i = 0
  let line = 1
  while (i < source.length) {
    const ch = source[i]
    if (ch === '\n') { line++; out += '\n'; i++; continue }

    if (ch === '#') {
      while (i < source.length && source[i] !== '\n') i++
      continue
    }
    if (ch === '"' || ch === "'") {
      const quote = ch
      const startLine = line
      i++
      let closed = false
      while (i < source.length) {
        if (source[i] === '\\') { i += 2; continue }
        if (source[i] === '\n') { line++; i++; continue }
        if (source[i] === quote) { closed = true; i++; break }
        i++
      }
      if (!closed) problems.push(`${file}:${startLine} 字符串未闭合 (${quote})`)
      out += '""'
      continue
    }
    if (ch === '`') {  // R 的反引号标识符，例如 `%||%`
      i++
      while (i < source.length && source[i] !== '`') i++
      i++
      out += 'X'
      continue
    }
    out += ch
    i++
  }
  return out
}

function checkBalance(code, file) {
  const pairs = { ')': '(', ']': '[', '}': '{' }
  const stack = []
  let line = 1
  for (let i = 0; i < code.length; i++) {
    const ch = code[i]
    if (ch === '\n') { line++; continue }
    if (ch === '(' || ch === '[' || ch === '{') stack.push({ ch, line })
    else if (ch in pairs) {
      const top = stack.pop()
      if (!top) problems.push(`${file}:${line} 多余的 '${ch}'`)
      else if (top.ch !== pairs[ch]) {
        problems.push(`${file}:${line} '${ch}' 与第 ${top.line} 行的 '${top.ch}' 不匹配`)
      }
    }
  }
  for (const item of stack) problems.push(`${file}:${item.line} '${item.ch}' 未闭合`)
}

/**
 * 相邻字符串字面量 —— R 没有隐式字符串拼接。
 *
 *   x <- ("第一段"
 *         "第二段")        # <- unexpected string constant
 *
 * Python / C 会把相邻字面量接起来，R **不会**，这是语法错误。
 * 而 `c("a", "b")` / `paste0("a", "b")` 是逗号分隔的多参数，合法。
 * 所以判据是：**两个字符串字面量之间除了空白/注释什么都没有**。
 *
 * **为什么必须在 stripLiterals 之前查：** 那个函数把每个字符串换成 `""`，
 * 于是 `("a" "b")` 变成 `("" "")`，括号配平完全正常 —— 这个错实测在
 * geo 侧从 check_r_syntax.mjs 眼皮底下溜过去，到云端 `source()` 才炸
 * （run 35482359507，09_export_targets.R:258）。
 */
function checkImplicitConcat(source, file) {
  let i = 0
  let line = 1
  while (i < source.length) {
    const ch = source[i]
    if (ch === '\n') { line++; i++; continue }
    if (ch === '#') { while (i < source.length && source[i] !== '\n') i++; continue }
    if (ch !== '"' && ch !== "'") { i++; continue }

    const quote = ch
    const startLine = line
    i++
    while (i < source.length) {
      if (source[i] === '\\') { i += 2; continue }
      if (source[i] === '\n') { line++; i++; continue }
      if (source[i] === quote) { i++; break }
      i++
    }

    let j = i
    let jLine = line
    for (;;) {
      const c = source[j]
      if (c === undefined) break
      if (c === '\n') { jLine++; j++; continue }
      if (c === ' ' || c === '\t' || c === '\r') { j++; continue }
      if (c === '#') { while (j < source.length && source[j] !== '\n') j++; continue }
      break
    }
    const next = source[j]
    if (next === '"' || next === "'") {
      problems.push(
        `${file}:${startLine} 相邻字符串字面量没有逗号 —— R 没有隐式字符串拼接，` +
        `用 paste0(...) 或 c(...) 显式连接（下一段在第 ${jLine} 行）`
      )
    }
  }
}

// ---- 检查：基础绘图函数不得写成 stats::（错误台账 E-03）--------------------
// 实测踩过两次：`stats::abline` / `stats::axis` —— 它们都在 **graphics** 命名空间。
// 本仓库 R 侧禁止 library()，每个函数都要写全名，所以写错命名空间会直接报
// "'X' is not an exported object from 'namespace:stats'"，而那个报错
// 不会提示"应该换成 graphics"。
const GRAPHICS_FNS = new Set([
  'plot', 'image', 'axis', 'abline', 'mtext', 'par', 'layout', 'legend',
  'text', 'title', 'lines', 'points', 'box', 'grid', 'rect', 'polygon',
  'hist', 'barplot', 'pie', 'contour', 'persp', 'pairs', 'matplot',
])
const STATS_WRONGLY = /stats::([A-Za-z_.][A-Za-z0-9_.]*)/g

// ---- 检查：sprintf 格式串里的裸 %（错误台账 E-02）--------------------------
// 实测：`sprintf("... top 2% of both ...", rho)` 报 `too few arguments` ——
// 文本里的 `% o` 被当成八进制转换符 `%o`，多吃一个参数。
//
// **判据要窄，否则误报淹没真信号。** 真正危险的只有一种模式：`%` 后面紧跟
// **空格**、再跟字母。合法的 `%%`、`%d`、`%.2f`、`%s`、`%-52s` 一律不报。
const BAD_PERCENT = /(?<!%)% +[a-zA-Z]/g
function checkSprintfPercent(source, file) {
  const lines = source.split(/\r?\n/)
  lines.forEach((line, i) => {
    if (!/sprintf\s*\(/.test(line)) return
    if (/^\s*#/.test(line)) return
    for (const m of line.matchAll(/"((?:[^"\\]|\\.)*)"/g)) {
      const lit = m[1]
      if (!lit.includes('%')) continue
      BAD_PERCENT.lastIndex = 0
      if (BAD_PERCENT.test(lit)) {
        problems.push(
          `${file}:${i + 1} sprintf 格式串里有裸 %（"${lit.slice(0, 40)}..."）—— ` +
          `文本里的 % 必须写 %%，否则 R 当成转换符多吃参数，报 too few arguments`)
      }
    }
  })
}

for (const file of files) {
  const rel = relative(root, file).split('\\').join('/')
  const source = readFileSync(file, 'utf8')
  const code = stripLiterals(source, rel)
  checkBalance(code, rel)
  // **必须在原始 source 上查**，不能用 stripLiterals 的结果（见函数注释）
  checkImplicitConcat(source, rel)
  checkSprintfPercent(source, rel)

  for (const m of code.matchAll(STATS_WRONGLY)) {
    if (GRAPHICS_FNS.has(m[1])) {
      problems.push(
        `${rel}: stats::${m[1]} 不存在 —— 基础绘图函数在 graphics 命名空间，` +
        `写 stats:: 会报 "not an exported object from 'namespace:stats'"`)
    }
  }
}

// ============================================================================
// 检查：R ↔ Python 的 CLI 契约
// ============================================================================
//
// R 脚本的 `DEFAULTS` 参数名、Python 调用点的 `--xxx` 拼写、R 写出的文件名、
// Python 读回的文件名 —— 这四样分布在两个语言里，**没有任何编译器能发现
// 它们不一致**。每一条的失败方式都很难从报错反推：
//
//   * Python 传了 R 不认识的参数 → R `stop("未知参数")`，运行时才炸；
//   * R 参数名拼错（`--n_nets` vs `--n_net`）→ 同上；
//   * R 写了 a.csv、Python 读 b.csv → Python 报"退出码 0 但没写出文件"，
//     看起来像 R 崩了，其实只是名字对不上。
//
// 这就是"两处定义必然分叉"那一类，必须在推送前静态挡掉。

// R 侧：DEFAULTS 里的参数名
const rParamNames = new Set()
// R 侧：**真正写进 out_dir 的**文件名。
// 只收 `file.path(out_dir, "xxx")` 里的字面量 —— 不能收全部字符串，
// 否则 R 自测里造临时文件的 `matrix.txt` / `meta.json`（在 tmpdir 下，
// 跟 Python 一点关系都没有）会被当成漏接的产物报出来。
const rOutputNames = new Set()
// Python 侧：传给 Rscript 的 --xxx
const pyFlags = new Set()
// Python 侧：读回的文件名
const pyReadNames = new Set()

for (const file of files) {
  const source = readFileSync(file, 'utf8')
  // 注释要先剥掉：文件头的用法示例里有 `--input in_matrix.txt`，
  // 那不是产物，但正则看得见。
  const code = source.replace(/^\s*#.*$/gm, '')

  // DEFAULTS <- list(input = NULL, meta = NULL, ..., selftest = FALSE)
  const d = /DEFAULTS\s*<-\s*list\(([\s\S]*?)\n\s*\)/.exec(code)
  if (d) {
    for (const m of d[1].matchAll(/(?:^|,)\s*([a-z_][a-z0-9_]*)\s*=/gm)) {
      rParamNames.add(m[1])
    }
  }

  for (const m of code.matchAll(/file\.path\(\s*out_dir\s*,\s*"([^"]+)"/g)) {
    rOutputNames.add(m[1])
  }
}

// Python 侧：找调用 Rscript 的那个文件
const pyFiles = []
function walkPy(dir) {
  for (const entry of readdirSync(dir)) {
    const full = join(dir, entry)
    if (statSync(full).isDirectory()) walkPy(full)
    else if (entry.endsWith('.py')) pyFiles.push(full)
  }
}
walkPy(targetDir)

/**
 * 取出 `cmd = [ ... ]` 这个列表的**完整**文本。
 *
 * 不能写成 /cmd\s*=\s*\[([\s\S]*?)\]/ —— 列表元素里就有 `]`：
 * `str(tparams["n_net"])` 的下标方括号会让非贪婪匹配在那里截断，
 * 于是**后面的参数全部收不到**，报出 5 条"Python 没传 --n_cells"的假阳性。
 * 所以必须按括号配平扫描。
 */
function extractCmdList(source) {
  const out = []
  for (const m of source.matchAll(/cmd\s*=\s*\[/g)) {
    let i = m.index + m[0].length
    let depth = 1
    const start = i
    while (i < source.length && depth > 0) {
      const ch = source[i]
      if (ch === '"' || ch === "'") {         // 跳过字符串字面量
        const q = ch
        i++
        while (i < source.length && source[i] !== q) {
          if (source[i] === '\\') i++
          i++
        }
        i++
        continue
      }
      if (ch === '[') depth++
      else if (ch === ']') depth--
      if (depth === 0) break
      i++
    }
    out.push(source.slice(start, i))
  }
  return out
}

for (const file of pyFiles) {
  const source = readFileSync(file, 'utf8')
  // **只在真正拼 Rscript 命令的那段里收参数。**
  // 收全文件的 `"--xxx"` 会把 main_analysis.py 的 argparse 选项
  // （`--config` / `--steps`）也算进来，然后报"R 不认识 --config"。
  // 判据：列表里出现 TENIFOLD_R 的才算调用 R 脚本的那一段。
  for (const block of extractCmdList(source)) {
    if (!/TENIFOLD_R/.test(block)) continue
    for (const f of block.matchAll(/"(--[a-z_][a-z0-9_]*)"/g)) pyFlags.add(f[1])
  }
  for (const m of source.matchAll(/["']([A-Za-z0-9_.-]+\.(?:csv|json|txt|tsv))["']/g)) {
    pyReadNames.add(m[1])
  }
}

if (rParamNames.size > 0 && pyFlags.size > 0) {
  // 1) Python 传的每个 --flag 都必须是 R 认识的参数
  for (const f of [...pyFlags].sort()) {
    const name = f.slice(2)
    if (!rParamNames.has(name)) {
      problems.push(
        `Python 传了 ${f}，但 R 脚本的 DEFAULTS 里没有 ${name} —— ` +
        `R 侧会 stop("未知参数")。R 认识: ${[...rParamNames].sort().join(', ')}`)
    }
  }
  // 2) R 的每个参数（除 selftest）都应当被 Python 传到
  for (const name of [...rParamNames].sort()) {
    if (name === 'selftest') continue
    if (!pyFlags.has(`--${name}`)) {
      problems.push(
        `R 脚本的 DEFAULTS 里有参数 ${name}，但 Python 调用点没有传 --${name} —— ` +
        `它会静默使用 R 侧的默认值，而那个默认值不在配置里，改不动`)
    }
  }
}

// 3) 文件名契约：R 写出的每个文件，Python 必须读回；反之亦然。
//    只比 basename（两侧一个拼 tmp 目录、一个拼 out_dir）。
if (rOutputNames.size > 0 && pyReadNames.size > 0) {
  const rBase = new Set([...rOutputNames].map(n => n.split('/').pop()))
  const pyBase = new Set([...pyReadNames].map(n => n.split('/').pop()))
  // Python 侧也写文件（落盘给用户），那些不是"从 R 读回"的，所以只查
  // **R 写了但 Python 完全没提到**的方向 —— 那一定是漏接了。
  for (const n of [...rBase].sort()) {
    if (!pyBase.has(n)) {
      problems.push(
        `R 脚本会写出 ${n}，但 Python 侧没有任何地方提到它 —— ` +
        `这份产物要么被丢掉了，要么文件名拼错了。R 产出: ${[...rBase].sort().join(', ')}`)
    }
  }
}

console.log(`检查了 ${files.length} 个 R 文件`)
console.log(`R 参数: ${[...rParamNames].sort().join(', ') || '（无）'}`)
console.log(`Python 传参: ${[...pyFlags].sort().join(', ') || '（无）'}`)

if (problems.length > 0) {
  console.log(`\n发现 ${problems.length} 个问题:`)
  for (const p of problems) console.log(`  - ${p}`)
  process.exit(1)
}
console.log('\n静态检查通过（注意：这不等于 R 能跑通，仍需真实执行验证）')
