#!/usr/bin/env node
// ============================================================================
// tools/check_doc_refs.mjs — 文档里引用的仓库文件必须真的存在（不需要 R / Python）
// ============================================================================
// 为什么需要它：**这条规则是被一个真实的手工发现逼出来的。**
//
// 实测 `AGENTS.md` 规则 29 写着「`_smoke_survival_diag.py` 把 pick_horizons /
// KM-at-horizon / cut 分组逐行转写，20 条断言」—— 而那个文件**早就不在仓库里**
// （临时脚本，验证完就删了）。`_smoke_export_targets.py` 同样。
//
// 危害不是"读者白找一次"：**下一个 agent 会照着文档去执行一个不存在的脚本**，
// 或者更糟 —— 以为那一步的验证没做过而重做一遍。
//
// 而这一类错误**当时没有任何门禁能挡住**，因为三个仓库的 workflow 都有
// `paths:` 过滤（只跑 scripts/ tools/ assets/ requirements.txt），
// **`*.md` 的改动根本不触发 CI**。所以文档里指向已删文件的死链接
// 可以一直躺着，CI 每次都是绿的。
//
// 现在配一个独立的 `docs_check.yml`（只有这一步，约 20 秒），
// 它才是"文档改动也要过 CI"的落点。
//
// ---------------------------------------------------------------------------
// 判据（两条，都刻意选高精度、低误报）
// ---------------------------------------------------------------------------
// A. **带仓库目录前缀的路径**：`scripts/…` / `tools/…` / `assets/…` /
//    `references/…` / `.github/…` → 必须存在。
//    这些前缀下的东西**都是提交进仓库的**，不存在就是死链接。
//    （`data/` 与 `results/` 故意不在名单里 —— 它们是运行时目录、被 gitignore，
//    文档里提到它们是正确的，不该报错。）
//
// B. **裸文件名 + 源码类扩展名**：`foo.py` / `foo.mjs` / `foo.R` / `foo.sh` /
//    `foo.md` / `foo.yml` → 必须在仓库里某个位置存在（按 basename 匹配）。
//    **只收源码/配置类扩展名**：`.csv` / `.json` / `.h5ad` / `.png` / `.pdf`
//    是**运行时产物**，文档里提到它们是预期行为（它们本来就不在仓库里）。
//    **必须以字词字符开头** —— 否则 `.mplstyle` 这种"光一个扩展名"的写法
//    会被当成文件名（实测误报过）。
//
// ---------------------------------------------------------------------------
// 逃生舱：怎么合法地提到一个"已经删掉的"或"别人家的"文件
// ---------------------------------------------------------------------------
// 在提到该文件的**同一行、或前后两行之内**写出下面任一标记即可
// （**是 ±2 行不是 ±1** —— 中文 Markdown 一句话常折三行，
//  解释落在第 2 行外会被误判成死链接，实测踩过）：
//
//   1. 文件已删：  已删 / 已删除 / 已移除 / 删了 / 删掉 / 不再存在 /
//                  曾经 / 当时的 / 旧名 / 原名 / deleted / removed
//   2. 第三方源码：第三方 / 上游 / 包内 / 源码 / site-packages / 库里 /
//                  该包 / 安装的
//
// 例：
//
//     > 两个当时的转写脚本 `_smoke_survival_diag.py` / `_smoke_export_targets.py`
//     > 已经删掉了（临时的，验证完就清）。
//
//     scFates 包内的 `pseudotime.py` 第 254 行才写入 milestones。
//
// **两个逃生舱都故意做成"要写一句话"的** ——
// 指向已删文件时必须说明它已删，否则读者无从判断是笔误还是历史；
// 指向第三方源码时必须说明是哪个包，否则读者不知道去 site-packages 里翻什么
// （**实测这一条本来就该补**：原文只写 `models.py`，没说那是 SpaGCN 的）。
// 静默豁免会让检查退化成没有检查。
//
// ---------------------------------------------------------------------------
// 已知的自指陷阱：**写这个工具的文档时会触发它自己**
// ---------------------------------------------------------------------------
// 上面这些"实测踩到过 X"的举例里，X 往往就是一个已经不存在的文件名 ——
// 于是文档本身成了死链接。实测两轮都栽在这里（`_smoke_survival_diag.py` /
// `assets/config.yml`）。
//
// **这不是 bug，是这条判据在正常工作** —— 读者确实会去找那个文件。
// 解法就是在举例时把话说全（"当时的错名" / "早就删了"），
// 而那本来就该说 —— 举例时不加限定，读者分不清你在说历史还是现状。
// ============================================================================

import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs'
import { join, dirname, basename, resolve } from 'node:path'

const ROOT = resolve(dirname(new URL(import.meta.url).pathname.replace(/^\/([A-Za-z]:)/, '$1')), '..')

// 扫描哪些文档
const DOCS = []
const addDoc = p => { if (existsSync(join(ROOT, p))) DOCS.push(p) }
for (const f of ['AGENTS.md', 'README.md', 'SKILL.md', 'EXPERIMENTAL_DESIGN.md',
                 'REFERENCES_VERIFICATION.md', 'CONTRIBUTING.md']) addDoc(f)
const refDir = join(ROOT, 'references')
if (existsSync(refDir)) {
  for (const f of readdirSync(refDir).sort()) {
    if (f.endsWith('.md')) DOCS.push(`references/${f}`)
  }
}

if (DOCS.length === 0) {
  console.error('没有找到任何文档（AGENTS.md / README.md / references/*.md）—— 检查器没东西可查')
  process.exit(2)
}

// ---- 仓库内所有文件的 basename（给判据 B 用）--------------------------------
// 跳过运行时目录与 VCS 元数据：那里的文件**不是仓库内容**，
// 拿它们当"存在"的证据会让判据 B 形同虚设。
const SKIP_DIRS = new Set(['.git', 'node_modules', 'data', 'results', '__pycache__',
                           '.venv', 'venv', 'dist', 'build', '.pytest_cache', '.Rproj.user'])
const basenames = new Set()
const relPaths = new Set()
;(function walk(dir, rel) {
  for (const name of readdirSync(dir)) {
    if (SKIP_DIRS.has(name)) continue
    const abs = join(dir, name)
    const r = rel ? `${rel}/${name}` : name
    let st
    try { st = statSync(abs) } catch { continue }
    if (st.isDirectory()) walk(abs, r)
    else { basenames.add(name); relPaths.add(r) }
  }
})(ROOT, '')

// ---- 判据 A 的前缀（都是提交进仓库的目录）----------------------------------
const PREFIXES = ['scripts/', 'tools/', 'assets/', 'references/', '.github/']

// ---- 判据 B 的形态（只收源码/配置，不收运行时产物；必须字词字符开头）--------
const BARE_SRC = /^[A-Za-z0-9_][A-Za-z0-9_.+-]*\.(py|mjs|js|ts|R|r|sh|md|yml|yaml|mplstyle|cff)$/

// ---- 逃生舱标记 -------------------------------------------------------------
// 1. 指向**已删**的文件：读者需要知道这是历史，不是笔误。
//    含"删了/删掉"这类口语说法 —— 它们和"已删"一样明确，
//    而写文档的人更可能顺手写成"早就删了"（实测就吃了这个亏）。
const DELETED_MARK = /已删|已删除|已移除|删了|删掉|不再存在|曾经|当时的|旧名|原名|deleted|removed/i
// 2. 指向**第三方包**的源码：读者需要知道去哪个包里找。
const THIRD_PARTY_MARK = /第三方|上游|包内|源码|site-packages|库里|该包|安装的/i
const EXEMPT_MARK = new RegExp(`${DELETED_MARK.source}|${THIRD_PARTY_MARK.source}`, 'i')

// 反引号里的内容。文档里的路径一律写成 `...`，这是本仓库的既定风格。
const TICK = /`([^`\n]+)`/g

// 去掉 token 尾部的标点（中文全角也要处理）
const trimTail = s => s.replace(/[，。、；：）】》,;:)\]}>.]+$/u, '')

const problems = []
let nCheckedA = 0
let nCheckedB = 0
let nExempt = 0

for (const doc of DOCS) {
  const lines = readFileSync(join(ROOT, doc), 'utf8').split(/\r?\n/)
  lines.forEach((line, i) => {
    for (const m of line.matchAll(TICK)) {
      let tok = trimTail(m[1].trim())
      if (!tok) continue
      // 占位符 / glob / 明显不是路径的，一律跳过
      if (/[<>{}*?$|]/.test(tok)) continue
      if (/\s/.test(tok)) continue
      // 绝对路径 / URL 跳过
      if (/^(https?:|\/|[A-Za-z]:\\)/.test(tok)) continue
      // 带单引号或方括号的是代码片段（layers['counts'] / obs["celltype"]），不是路径
      if (/['"[\]]/.test(tok)) continue

      let isA = false
      let isB = false
      if (PREFIXES.some(p => tok.startsWith(p))) isA = true
      else if (!tok.includes('/') && BARE_SRC.test(tok)) isB = true
      if (!isA && !isB) continue

      // 判据 B 只按 basename 找；判据 A 按完整相对路径找
      const found = isA
        ? (existsSync(join(ROOT, tok)) || relPaths.has(tok))
        : basenames.has(tok)
      if (found) { isA ? nCheckedA++ : nCheckedB++; continue }

      // 没找到 —— 看逃生舱。窗口是 **±2 行**，不是 ±1。
      //
      // **为什么是 2：** 中文 Markdown 按 ~40 字硬折行，一句话经常占三行。
      // 实测 spatial 的 AGENTS.md 里 `models.py` 在第 487 行，而解释它的
      // "（它们是 SpaGCN 与 STAGATE_pyG）**包内**的文件"落在第 489 行 ——
      // ±1 看不见它，一条**本来就写清楚了**的引用被判成死链接。
      // 窗口太窄的后果不是"更严格"，是"逼人把话说得更碎"。
      //
      // **为什么不用更大：** 实测 ±2 相对 ±1 只多豁免**一条**（就是上面那条），
      // 而 ±3 之后"源码"这类常见词会开始误伤 —— 窗口越宽，
      // 静默豁免越容易发生，而静默豁免会让检查退化成没有检查。
      const window = [
        lines[i - 2] || '', lines[i - 1] || '', line, lines[i + 1] || '', lines[i + 2] || '',
      ].join('\n')
      if (EXEMPT_MARK.test(window)) { nExempt++; continue }

      problems.push({ doc, line: i + 1, tok, rule: isA ? 'A' : 'B' })
    }
  })
}

// ---- 报告 -------------------------------------------------------------------
const total = nCheckedA + nCheckedB
console.log(`文档引用检查：${DOCS.length} 份文档`)
console.log(`  判据 A（带仓库目录前缀的路径）  ${nCheckedA} 条通过`)
console.log(`  判据 B（裸源码/配置文件名）     ${nCheckedB} 条通过`)
if (nExempt) console.log(`  逃生舱（前后两行内标了"已删"或"第三方包"）  ${nExempt} 条豁免`)

if (problems.length === 0) {
  console.log(`\n全部 ${total} 条引用都指向真实存在的文件。`)
  process.exit(0)
}

console.error(`\n${problems.length} 条引用指向**不存在的文件**：\n`)
for (const p of problems) {
  console.error(`  ${p.doc}:${p.line}  [判据 ${p.rule}]  ${p.tok}`)
}
console.error(`
死链接比"少写一句话"更糟：读者会去找一个不存在的文件，
或者以为那一步没做过而重做一遍。

三种改法，三选一：
  1. 引用改成真实存在的文件（多半是路径写错或文件改名了）；
  2. 那个文件**确实已经删掉** → 在它的前后两行之内写明
     "已删 / 已删除 / 删了 / 已移除 / 不再存在 / 曾经 / 当时的 / deleted"。
     **只留结论，不留死链接。**
  3. 那是**第三方包的源码**（不是本仓库的） → 写明是哪个包，
     用"包内 / 上游 / 源码 / site-packages / 该包"之类的话。
     只写 \`models.py\` 读者不知道去哪个包里翻。
`)
process.exit(1)
