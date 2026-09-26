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

import { readFileSync, existsSync, readdirSync, statSync,
         writeFileSync, mkdirSync, mkdtempSync, rmSync } from 'node:fs'
import { join, dirname, basename, resolve } from 'node:path'
import { tmpdir } from 'node:os'

const ROOT = resolve(dirname(new URL(import.meta.url).pathname.replace(/^\/([A-Za-z]:)/, '$1')), '..')

// 扫描哪些文档
function listDocs(root) {
  const DOCS = []
  const addDoc = p => { if (existsSync(join(root, p))) DOCS.push(p) }
  for (const f of ['AGENTS.md', 'README.md', 'SKILL.md', 'EXPERIMENTAL_DESIGN.md',
                   'REFERENCES_VERIFICATION.md', 'CONTRIBUTING.md']) addDoc(f)
  const refDir = join(root, 'references')
  if (existsSync(refDir)) {
    for (const f of readdirSync(refDir).sort()) {
      if (f.endsWith('.md')) DOCS.push(`references/${f}`)
    }
  }
  return DOCS
}

// ---- 仓库内所有文件的 basename（给判据 B 用）--------------------------------
// 跳过运行时目录与 VCS 元数据：那里的文件**不是仓库内容**，
// 拿它们当"存在"的证据会让判据 B 形同虚设。
const SKIP_DIRS = new Set(['.git', 'node_modules', 'data', 'results', '__pycache__',
                           '.venv', 'venv', 'dist', 'build', '.pytest_cache', '.Rproj.user'])
function indexRepo(dir) {
  const basenames = new Set()
  const relPaths = new Set()
  ;(function walk(d, rel) {
    for (const name of readdirSync(d)) {
      if (SKIP_DIRS.has(name)) continue
      const abs = join(d, name)
      const r = rel ? `${rel}/${name}` : name
      let st
      try { st = statSync(abs) } catch { continue }
      if (st.isDirectory()) walk(abs, r)
      else { basenames.add(name); relPaths.add(r) }
    }
  })(dir, '')
  return { basenames, relPaths }
}
const { basenames, relPaths } = indexRepo(ROOT)

// ---- 判据 C：**跨仓**引用（workspace 里的姊妹仓库）---------------------------
// 三个仓库的 AGENTS.md 互相点名是常态，写法是
// `scrna-pipeline-skill/scripts/03_spatial_domains.py`。这种 token
// **以前两条判据都不进**（不以 scripts/ 开头、又含 `/` 所以不是判据 B）
// —— 于是**完全不被检查**，死链接可以一直躺着而 CI 是绿的。
// 与 E-62 / E-63 同一类：检查器"没报错"不等于"检查了"。
//
// **但 CI 只 checkout 本仓库**，姊妹目录不存在 —— 那时必须判"**不适用**"
// 而不是"失败"（同 `check_py_names.py` 对纯 R 仓库的处理：判红会让每次 CI
// 都红，而那条告警与被检查的改动毫无关系）。所以：**姊妹目录存在才查**，
// 不存在就跳过，并在报告里写明跳过了多少条 —— 让"没查"和"查过没问题"
// 长得不一样（AGENTS 规则 4 的同一条理由）。
//
// > **姊妹仓库索引必须建在判据 A/B 之前，且判据 C 的分支必须写在
// > `if (!isA && !isB) continue` 之前** —— 第一版把它写在后面，
// > 于是"带仓名前缀"那一条永远走不到（它在 L219 就被 `continue` 掉了）。
// > 见下方判据 C 主体处的详细记录。
function buildSiblings(root) {
  const WS = resolve(root, '..')
  const siblings = new Map() // 目录名 -> { basenames, relPaths }
  try {
    for (const name of readdirSync(WS)) {
      if (name === basename(root) || SKIP_DIRS.has(name)) continue
      const abs = join(WS, name)
      let st
      try { st = statSync(abs) } catch { continue }
      if (!st.isDirectory()) continue
      // **只认真正的仓库**（有 `.git`）。workspace 根下还有 releases/ 等
      // 非仓库目录，它们的结构和仓库很像（`releases/…-v0.1.1/.github/workflows/`），
      // 按 basename 会把一条"缺仓名前缀"的引用**指到一个打包副本上** ——
      // 实测就这样把 `geo_analysis.yml` 提示成了
      // `releases/workflows/geo_analysis.yml`（连路径都是错的）。
      // **指错地方比不指地方更糟**（同 E-61 防复发④）。
      if (!existsSync(join(abs, '.git'))) continue
      siblings.set(name, indexRepo(abs))
    }
  } catch { /* workspace 根不可读（CI 只 checkout 本仓库时是正常的） */ }
  return siblings
}

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

// ---- 扫描主体：**抽成纯函数**，`main()` 与 `--selftest` 都调它 --------------
// E-63 的教训（AGENTS 规则 30.1）：**自检里重实现一遍被测逻辑，等于没测**。
// 自检第一版如果自己另写一遍扫描循环，那么把缺陷注回真代码时自检照样通过。
// 所以这里必须是同一个函数，`root` / `DOCS` / `siblings` 全部走参数。
function scanDocs(root, DOCS, siblings) {
  const { basenames, relPaths } = indexRepo(root)
  const problems = []
  let nCheckedA = 0
  let nCheckedB = 0
  let nCheckedC = 0
  let nExempt = 0
  let nStripped = 0
  let nXrepoSkipped = 0

  for (const doc of DOCS) {
    const lines = readFileSync(join(root, doc), 'utf8').split(/\r?\n/)
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

      // ---- 先剥掉"行号后缀"再判存在性（E-57）-----------------------------
      // 「真实文件 + 行号」是本仓库三仓文档的**既定引用风格**
      // （`scripts/05_trajectory.py:275-279` 让读者直接去源码定位），
      // 但整串拿去 `existsSync` 必然为假。
      //
      // 而且它**在两条判据下错得相反**（实测）：
      //   - 带仓库前缀的 `scripts/…py:275-279` → 走判据 A → **假阳性**；
      //   - 裸文件名的 `01_qc.py:226` / `models.py:19` → `BARE_SRC` 要求以
      //     `.py` 结尾，`:19` 结尾不匹配 → **两条判据都不进、完全不被检查**。
      // 假阳性会训练人忽略告警；盲区则让"指向第三方包源码却没写包名"的引用
      // 一直躺着 —— 而那正是这条判据存在的理由（文件头已论证过）。
      //
      // 所以**两条判据统一**先剥后缀再判存在性，报告里仍打印**原文**。
      // 只剥末尾的行号：`:275` / `:275-279` / `:L275` / `:275,278`。
      // 中间带冒号的（`references/methods.md`）不受影响 —— 正则锚在 `$`。
      const stripped = tok.replace(/:(?:L)?\d+(?:\s*[-,]\s*\d+)*$/, '')
      if (stripped !== tok) nStripped++
      const probe = stripped

      // ---- 逃生舱的窗口先算出来，**判据 A/B/C 共用** ----------------------
      // 窗口是 **±2 行**，不是 ±1。
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

      // ---- 判据 C：跨仓引用 ------------------------------------------------
      // 三个仓库的 AGENTS.md 互相点名是常态，两种写法：
      //
      //   `scrna-pipeline-skill/tools/check_figures.mjs`  ← 带仓名前缀（正确）
      //   `.github/workflows/geo_analysis.yml`            ← 没带前缀（读者不知道去哪找）
      //
      // **第一种两条判据都不进**：不以 `scripts/` 等开头（不是判据 A），
      // 又含 `/` 所以不是判据 B。**它必须在这里就被显式接住** ——
      // 否则会在下面那句 `if (!isA && !isB) continue` 上直接跳过，
      // 判据 C 的"带仓名前缀"分支成了**死代码**。
      //
      // > 这个坑实测踩过：第一版判据 C 写在 `continue` 之后，
      // > 三仓跑下来 `nCheckedC` **恒为 0**，而报告里那行是
      // > `if (nCheckedC)` 守卫的 —— 于是它**连打印都不打印**，
      // > 看起来跟"没有跨仓引用"一模一样。与 E-62 / E-63 同一类：
      // > 检查器"没报错"不等于"检查了"，而且**假阴性还会顺手把自己藏起来**。
      //
      // 姊妹目录**不存在时不查**（CI 只 checkout 本仓库）：那是"不适用"，
      // 不是"失败" —— 判红会让每次 CI 都红，而那条告警与本次改动无关。
      // 报告里分别计数，让"没查"和"查过没问题"长得不一样。
      const slash = probe.indexOf('/')
      const head = slash === -1 ? null : probe.slice(0, slash)
      const base = basename(probe)
      const isXrepo = head !== null && siblings.has(head)

      if (PREFIXES.some(p => probe.startsWith(p))) isA = true
      else if (!probe.includes('/') && BARE_SRC.test(probe)) isB = true

      if (!isA && !isB && !isXrepo) {
        // 含斜杠、又不是判据 A/B 的，绝大多数**本来就不该查**：
        //   - 运行时目录：`data/` / `results/` / `figures/`
        //   - 第三方包内部路径：`SpaGCN/SpaGCN.py` / `STAGATE_pyG/gat_conv.py`
        //   - 无 scheme 的 URL：`r-lib/actions/setup-r-dependencies@v2`
        //   - 代码片段：`if/else` / `72/fig.dpi` / `OMP/MKL/OPENBLAS_NUM_THREADS`
        // 这些在判据 C 之前就一直是"不进任何判据"的，静默跳过即可。
        //
        // **但姊妹目录不存在时性质就不同了**：那时"跨仓引用"这一类
        // **整体无法判断**，必须如实计数并在报告里写出来 ——
        // 否则"没查"和"查过没问题"长得一模一样（AGENTS 规则 4）。
        // 两种情况**不能共用一个计数器**：混在一起会打印出
        // "姊妹仓库目录不存在"这句**假话**（实测踩到，见报告段注释）。
        if (slash !== -1 && siblings.size === 0) nXrepoSkipped++
        continue
      }

      // 带仓名前缀：直接去那个仓里找（**必须在 `found` 之前**，因为
      // `scrna-pipeline-skill/...` 在本仓库里必然找不到）。
      if (isXrepo) {
        const sib = siblings.get(head)
        const rest = probe.slice(slash + 1)
        if (sib.relPaths.has(rest) || sib.basenames.has(base)) { nCheckedC++; continue }
        if (EXEMPT_MARK.test(window)) { nExempt++; continue }
        problems.push({ doc, line: i + 1, tok, rule: 'C' })
        continue
      }

      // 判据 B 只按 basename 找；判据 A 按完整相对路径找
      const found = isA
        ? (existsSync(join(ROOT, probe)) || relPaths.has(probe))
        : basenames.has(probe)
      if (found) { isA ? nCheckedA++ : nCheckedB++; continue }

      // 本仓库找不到 —— 也许它躺在某个姊妹仓库里，只是**没写仓名前缀**。
      // 这时报"缺仓名前缀"而不是"文件不存在"：真实原因是读者不知道去
      // 哪个仓找，而三个仓的目录结构很像（都有 `.github/workflows/`）。
      //
      // 提示里给**真实相对路径**，不是拿 token 拼出来的猜测路径 ——
      // 拼出来的那条实测是错的（`releases/workflows/geo_analysis.yml`）。
      if (siblings.size > 0) {
        const homes = [...siblings.entries()]
          .map(([n, s]) => {
            if (s.relPaths.has(probe)) return `${n}/${probe}`
            const hit = [...s.relPaths].find(r => basename(r) === base)
            return hit ? `${n}/${hit}` : null
          })
          .filter(Boolean)
        if (homes.length > 0) {
          if (EXEMPT_MARK.test(window)) { nExempt++; continue }
          problems.push({ doc, line: i + 1, tok, rule: 'C', homes })
          continue
        }
      }

      // 没找到 —— 看逃生舱。
      if (EXEMPT_MARK.test(window)) { nExempt++; continue }

      problems.push({ doc, line: i + 1, tok, rule: isA ? 'A' : 'B' })
      }
    })
  }
  return { problems, nCheckedA, nCheckedB, nCheckedC, nExempt, nStripped, nXrepoSkipped }
}

// ---- 报告 -------------------------------------------------------------------
// **判据 C 那一行不能加 `if (nCheckedC)` 守卫。** 第一版加了，于是判据 C
// 因为死代码恒为 0 时，它**连打印都不打印** —— 看起来跟"本仓没有跨仓引用"
// 完全一样，缺陷自己把自己藏了起来。**计数为 0 本身就是需要看见的信息**，
// 何况 0 与"没跑"必须长得不一样（同 `check_py_names.py` 的"不适用 ≠ 失败"）。
function report(nDocs, r) {
  const total = r.nCheckedA + r.nCheckedB + r.nCheckedC
  console.log(`文档引用检查：${nDocs} 份文档`)
  console.log(`  判据 A（带仓库目录前缀的路径）  ${r.nCheckedA} 条通过`)
  console.log(`  判据 B（裸源码/配置文件名）     ${r.nCheckedB} 条通过`)
  console.log(`  判据 C（跨仓引用，姊妹仓库里找）  ${r.nCheckedC} 条通过`)
  if (r.nXrepoSkipped) {
    console.log(`  **跨仓引用未查**：姊妹仓库目录不存在（CI 只 checkout 本仓库）`
      + ` —— 本次有 ${r.nXrepoSkipped} 条含斜杠的引用**没被检查过**`)
  }
  if (r.nExempt) console.log(`  逃生舱（前后两行内标了"已删"或"第三方包"）  ${r.nExempt} 条豁免`)
  if (r.nStripped) console.log(`  带行号的引用（剥掉 \`:275-279\` 后缀后按文件判存在性）  ${r.nStripped} 条`)

  if (r.problems.length === 0) {
    console.log(`\n全部 ${total} 条引用都指向真实存在的文件。`)
    return 0
  }

  console.error(`\n${r.problems.length} 条引用指向**不存在的文件**：\n`)
  for (const p of r.problems) {
    const hint = p.homes ? `\n      ← 它在姊妹仓库里：${p.homes.join(' 或 ')}` : ''
    console.error(`  ${p.doc}:${p.line}  [判据 ${p.rule}]  ${p.tok}${hint}`)
  }
  console.error(`
死链接比"少写一句话"更糟：读者会去找一个不存在的文件，
或者以为那一步没做过而重做一遍。

四种改法，四选一：
  1. 引用改成真实存在的文件（多半是路径写错或文件改名了）；
  2. 那个文件**确实已经删掉** → 在它的前后两行之内写明
     "已删 / 已删除 / 删了 / 已移除 / 不再存在 / 曾经 / 当时的 / deleted"。
     **只留结论，不留死链接。**
  3. 那是**第三方包的源码**（不是本仓库的） → 写明是哪个包，
     用"包内 / 上游 / 源码 / site-packages / 该包"之类的话。
     只写 \`models.py\` 读者不知道去哪个包里翻。
  4. 那是**姊妹仓库**的文件（判据 C） → 写上仓库目录名，写成
     \`scrna-pipeline-skill/tools/check_figures.mjs\`。
     **不写仓名读者不知道去哪个仓找** —— 而三个仓的目录结构很像，
     他多半会在本仓里找一个同名文件然后放弃。
`)
  return 1
}

// ============================================================================
// --selftest：**带内建自检**（AGENTS 规则 30，E-62 / E-63 的产物）
// ============================================================================
// 为什么必须有：判据 C 的第一版就是**死代码** —— 它写在
// `if (!isA && !isB) continue` 之后，于是"带仓名前缀"那条永远走不到。
// 三仓跑下来 `nCheckedC` 恒为 0，而报告里那行有 `if (nCheckedC)` 守卫，
// **连打印都不打印**，看起来与"本仓没有跨仓引用"完全一样。
//
// **正向跑绿证明不了任何事** —— 死代码在干净文档上也是绿的。
// 所以自检必须：
//   ① 调**真代码**（`scanDocs` / `report`，不是自己再写一遍扫描）；
//   ② 在**隔离的 scratch workspace** 里造正反用例；
//   ③ 每类缺陷只让它自己那条响。
//
// 隔离靠 `scanDocs(root, DOCS, siblings)` 三个参数 —— 真仓库的 ROOT 是
// 由脚本位置推出的常量，自检要能指向临时目录，所以不能是全局量。
function selftest() {
  const cases = []
  const tmp = mkdtempSync(join(tmpdir(), 'docrefselftest-'))
  const WSR = join(tmp, 'ws')
  const ME = join(WSR, 'probe-repo')

  const mk = (p, txt) => { mkdirSync(dirname(p), { recursive: true }); writeFileSync(p, txt, 'utf8') }

  // 本仓
  mk(join(ME, 'AGENTS.md'), '# probe\n')
  mk(join(ME, 'scripts', '01_qc.py'), 'pass\n')
  mk(join(ME, 'references', 'methods.md'), '# m\n')
  // 姊妹仓 A（真仓库）
  mk(join(WSR, 'sibling-a', '.git', 'keep'), '')
  mk(join(WSR, 'sibling-a', 'tools', 'check_figures.mjs'), '// x\n')
  mk(join(WSR, 'sibling-a', '.github', 'workflows', 'alpha_analysis.yml'), 'x\n')
  // 姊妹仓 B（真仓库）
  mk(join(WSR, 'sibling-b', '.git', 'keep'), '')
  mk(join(WSR, 'sibling-b', '.github', 'workflows', 'beta_analysis.yml'), 'x\n')
  // 伪装成仓库的**非** git 目录（真 workspace 里 releases/ 就是这样）
  mk(join(WSR, 'releases', '.github', 'workflows', 'gamma_analysis.yml'), 'x\n')

  const run = text => {
    mk(join(ME, 'AGENTS.md'), text)
    const docs = listDocs(ME)
    const sib = buildSiblings(ME)
    return { r: scanDocs(ME, docs, sib), sib }
  }
  const codes = r => r.problems.map(p => p.rule).sort().join(',')

  try {
    // ---- 用例 1：干净文档零死链接 -------------------------------------------
    {
      const { r } = run('# probe\n本仓：`scripts/01_qc.py`、`references/methods.md`。\n'
        + '跨仓：`sibling-a/tools/check_figures.mjs`、'
        + '`sibling-b/.github/workflows/beta_analysis.yml`。\n')
      const ok = r.problems.length === 0 && r.nCheckedA === 2 && r.nCheckedC === 2
      cases.push(['干净文档零死链接，且判据 A=2 / C=2', ok,
        `problems=${r.problems.length} A=${r.nCheckedA} C=${r.nCheckedC}`])
    }

    // ---- 用例 2：**回归：判据 C 带前缀分支不是死代码** ----------------------
    // 这一条就是 E-64 的回归。把 `const isXrepo = ...` 改成 `false` 会让它红。
    {
      const { r } = run('# probe\n跨仓：`sibling-a/tools/check_figures.mjs`。\n')
      cases.push(['回归 1：带仓名前缀的跨仓引用真的进了判据 C', r.nCheckedC === 1,
        `C=${r.nCheckedC}`])
    }

    // ---- 用例 3：跨仓死链接判红 ---------------------------------------------
    {
      const { r } = run('# probe\n跨仓：`sibling-a/tools/does_not_exist.mjs`。\n')
      cases.push(['跨仓死链接判红（判据 C）', codes(r) === 'C', `codes=${codes(r)}`])
    }

    // ---- 用例 4：缺仓名前缀 → 报 C 并给真实姊妹路径 -------------------------
    {
      const { r } = run('# probe\n跨仓：`.github/workflows/alpha_analysis.yml`。\n')
      const p = r.problems[0]
      const ok = codes(r) === 'C' && p && p.homes
        && p.homes.includes('sibling-a/.github/workflows/alpha_analysis.yml')
      cases.push(['缺仓名前缀 → 报 C 且给出真实姊妹路径', !!ok,
        p ? `homes=${JSON.stringify(p.homes)}` : 'no problem'])
    }

    // ---- 用例 5：**非 git 目录不算姊妹仓**（releases/ 那类）-----------------
    // 注意 `.github/workflows/...` 本身是**判据 A**（前缀 `.github/`），所以
    // 这里的期望是 `A` —— 关键在**提示里不能出现 releases/**：
    // 修之前它会拿 basename 在 releases/ 里反查到一条路径并指过去。
    {
      const { r } = run('# probe\n跨仓：`.github/workflows/gamma_analysis.yml`。\n')
      const hinted = JSON.stringify(r.problems.map(p => p.homes || []))
      const ok = codes(r) === 'A' && !hinted.includes('releases')
      cases.push(['非 git 目录（releases/）不算姊妹仓，提示里不出现它', ok,
        `codes=${codes(r)} hint=${hinted}`])
    }

    // ---- 用例 6：判据 A 死链接仍然红（没被判据 C 抢走）----------------------
    {
      const { r } = run('# probe\n本仓：`scripts/99_missing.py`。\n')
      cases.push(['判据 A 死链接仍然判红', codes(r) === 'A', `codes=${codes(r)}`])
    }

    // ---- 用例 7：判据 B 死链接仍然红 ----------------------------------------
    {
      const { r } = run('# probe\n本仓：`99_missing.py`。\n')
      cases.push(['判据 B 死链接仍然判红', codes(r) === 'B', `codes=${codes(r)}`])
    }

    // ---- 用例 8：逃生舱对判据 C 也生效 --------------------------------------
    {
      const { r } = run('# probe\n跨仓：`.github/workflows/alpha_analysis.yml`（**已删**，当年在 sibling-a 里）。\n')
      const ok = r.problems.length === 0 && r.nExempt === 1
      cases.push(['逃生舱（已删）对判据 C 生效', ok,
        `problems=${r.problems.length} exempt=${r.nExempt}`])
    }

    // ---- 用例 9：姊妹仓不存在时不判红，且**如实打印"未查"** -----------------
    {
      const lonely = join(tmp, 'lonely', 'me')
      mk(join(lonely, 'AGENTS.md'), '# probe\n跨仓：`sibling-a/tools/check_figures.mjs`。\n')
      const docs = listDocs(lonely)
      const sib = buildSiblings(lonely)
      const r = scanDocs(lonely, docs, sib)
      const ok = r.problems.length === 0 && r.nXrepoSkipped === 1 && r.nCheckedC === 0
      cases.push(['姊妹仓不存在 → 不判红，且计数为"未查"', ok,
        `problems=${r.problems.length} skipped=${r.nXrepoSkipped}`])
    }

    // ---- 用例 10：报告里判据 C 那一行**计数为 0 也打印**（E-64 回归）--------
    // 第一版用 `if (nCheckedC)` 守卫，于是死代码把自己藏起来了。
    {
      const lonely = join(tmp, 'lonely2', 'me')
      mk(join(lonely, 'AGENTS.md'), '# probe\n本仓：`scripts/01_qc.py`。\n')
      mk(join(lonely, 'scripts', '01_qc.py'), 'pass\n')
      const r = scanDocs(lonely, listDocs(lonely), buildSiblings(lonely))
      // report() 直接打到 stdout，这里用一个**假的 console** 捕获它
      const realLog = console.log
      let buf = ''
      console.log = (...a) => { buf += a.join(' ') + '\n' }
      try { report(1, r) } finally { console.log = realLog }
      const ok = buf.includes('判据 C') && buf.includes('0 条通过')
      cases.push(['回归 2：判据 C 计数为 0 时那一行照样打印', ok,
        ok ? '' : `buf=${JSON.stringify(buf.slice(0, 200))}`])
    }

    // ---- 用例 11：无文档 → 调用方应报"没东西可查"（不静默通过）--------------
    {
      const empty = join(tmp, 'empty', 'me')
      mk(join(empty, 'notes.txt'), 'x\n')
      cases.push(['没有文档时 listDocs 返回空（调用方据此判 exit=2，不静默通过）',
        listDocs(empty).length === 0, ''])
    }
  } finally {
    try { rmSync(tmp, { recursive: true, force: true }) } catch { /* 临时目录清不掉不影响结论 */ }
  }

  let bad = 0
  for (const [name, ok, detail] of cases) {
    console.log(`  [${ok ? 'OK ' : 'FAIL'}] ${name}${ok || !detail ? '' : `  ← ${detail}`}`)
    if (!ok) bad++
  }
  if (bad) {
    console.error(`\n自检未通过（${bad}/${cases.length} 个用例失败）`)
    return 1
  }
  console.log(`自检通过（${cases.length} 个用例，含 2 条 E-64 回归）`)
  return 0
}

// ---- main -------------------------------------------------------------------
function main() {
  const DOCS = listDocs(ROOT)
  if (DOCS.length === 0) {
    console.error('没有找到任何文档（AGENTS.md / README.md / references/*.md）—— 检查器没东西可查')
    return 2
  }
  const siblings = buildSiblings(ROOT)
  return report(DOCS.length, scanDocs(ROOT, DOCS, siblings))
}

process.exit(process.argv.includes('--selftest') ? selftest() : main())
