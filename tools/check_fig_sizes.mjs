#!/usr/bin/env node
/**
 * 图幅检查：PDF 的实际宽度必须 <= W_DOUBLE（183 mm）。
 *
 * 用法：
 *   node tools/check_fig_sizes.mjs <图目录>
 *   例：node tools/check_fig_sizes.mjs results/lymph_node/figures
 *       node tools/check_fig_sizes.mjs results/GSE42568
 *
 * ── 为什么需要它 ──────────────────────────────────────────────
 *
 * AGENTS 规则（geo 25 / scrna 13 / spatial 13）要求**图幅按毫米**，
 * 宽度夹在 W_SINGLE 89 / W_ONE_HALF 136 / W_DOUBLE 183 之内。
 *
 * 但已有的两个检查都看不见"图有多宽"：
 *   - 验收（`check_acceptance()` / `check_acceptance`）：只判断文件**存在**
 *   - `check_figures.mjs`：只判断**有没有墨**
 * 一张 370 mm 的图和一张 183 mm 的图，在这两者眼里完全一样
 * （实测修之前最宽的 `domain_markers_dotplot` 就是 370 mm）。
 *
 * 所以这里直接读 PDF 的 `/MediaBox`，换算成毫米，跟 183 mm 比。
 * PDF 是矢量、尺寸精确，不受 dpi 影响 —— 比量 PNG 像素可靠。
 *
 * ── 判据 ─────────────────────────────────────────────────────
 *
 *   > 183 mm + 0.5 mm 容差   -> FAIL（装不进任何期刊的一页）
 *   <= 183 mm                -> PASS
 *   读不到 /MediaBox         -> FAIL（否则检查会静默通过，等于没检查）
 *   一个 PDF 都没找到         -> FAIL（同上）
 *
 * **不比"是不是正好等于 89 / 136 / 183"。** 比标准栏宽窄不是错误 ——
 * 本仓库有意用 `mm(165)` / `mm(178)` 这类值，那是双栏图没用满宽，
 * 排版上没问题。这里只拦"超出上限"这一件事。
 *
 * ── 已知的坑 ─────────────────────────────────────────────────
 *
 * R 的 `pdf()` 设备把 MediaBox **向下取整到整数点**，所以
 * `W_DOUBLE`（183 mm = 518.74 pt）落盘成 518 pt = 182.74 mm。
 * 容差就是为这类取整留的 —— 不要把 182.74 读成"少画了 0.26 mm"。
 *
 * Python 侧（matplotlib）不做这个取整，183 mm 就是 518.74 pt。
 * 同一个仓库里两种后端给出两种末位，是正常的。
 */

import fs from 'node:fs'
import path from 'node:path'

const PT_PER_MM = 72 / 25.4
const W_SINGLE_MM = 89
const W_ONE_HALF_MM = 136
const W_DOUBLE_MM = 183
const TOLERANCE_MM = 0.5

const target = process.argv[2]

if (!target) {
  console.error('用法: node tools/check_fig_sizes.mjs <图目录>')
  console.error('  例: node tools/check_fig_sizes.mjs results/lymph_node/figures')
  process.exit(2)
}

if (!fs.existsSync(target) || !fs.statSync(target).isDirectory()) {
  console.error(`图目录不存在或不是目录: ${target}`)
  process.exit(2)
}

/** 读 PDF 的 /MediaBox，返回 { widthMM, heightMM }；读不到返回 null。 */
function mediaBoxMM(file) {
  // latin1 保证字节到字符一一对应，不会因为非 UTF-8 内容丢字节。
  const text = fs.readFileSync(file).toString('latin1')
  const m = /\/MediaBox\s*\[\s*([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)\s*\]/.exec(text)
  if (!m) return null
  const x0 = parseFloat(m[1])
  const x1 = parseFloat(m[3])
  const y0 = parseFloat(m[2])
  const y1 = parseFloat(m[4])
  if (![x0, x1, y0, y1].every(Number.isFinite)) return null
  return {
    widthMM: Math.abs(x1 - x0) / PT_PER_MM,
    heightMM: Math.abs(y1 - y0) / PT_PER_MM,
  }
}

// `__<短SHA>` 是 workflow 事后写的带版本号副本，和规范名内容相同 —— 不重复检查。
const pdfs = fs
  .readdirSync(target)
  .filter((f) => f.toLowerCase().endsWith('.pdf') && !f.includes('__'))
  .sort()

if (pdfs.length === 0) {
  console.error(`[FAIL] ${target} 里一个 PDF 都没有 —— 检查没有实际执行，不能算通过`)
  process.exit(1)
}

const rows = []
const unreadable = []

for (const name of pdfs) {
  const mm = mediaBoxMM(path.join(target, name))
  if (!mm) {
    unreadable.push(name)
    continue
  }
  rows.push({ name, ...mm })
}

const over = rows.filter((r) => r.widthMM > W_DOUBLE_MM + TOLERANCE_MM)

// —— 报告 ——
console.log(`图幅检查：${target}`)
console.log(`  读到一个 ${pdfs.length} 个 PDF`)
console.log('')
for (const r of rows.sort((a, b) => b.widthMM - a.widthMM)) {
  const mark = r.widthMM > W_DOUBLE_MM + TOLERANCE_MM ? '[FAIL]' : '[OK]  '
  console.log(`  ${mark} ${r.name.padEnd(40)} ${r.widthMM.toFixed(1)} x ${r.heightMM.toFixed(1)} mm`)
}

// 分布：让日志里能直接看出图幅都落在哪一档
const buckets = new Map()
for (const r of rows) {
  const near = [W_SINGLE_MM, W_ONE_HALF_MM, W_DOUBLE_MM].find(
    (v) => Math.abs(r.widthMM - v) <= TOLERANCE_MM,
  )
  const key = near ? `${near} mm` : '其他宽度'
  buckets.set(key, (buckets.get(key) || 0) + 1)
}
console.log('')
console.log('  宽度分布:')
for (const k of ['89 mm', '136 mm', '183 mm', '其他宽度']) {
  if (buckets.has(k)) console.log(`    ${k.padEnd(10)} ${buckets.get(k)}`)
}

// —— 判定 ——
let failed = false

if (unreadable.length > 0) {
  console.error('')
  console.error(`[FAIL] ${unreadable.length} 个 PDF 读不到 /MediaBox：`)
  for (const n of unreadable) console.error(`    ${n}`)
  console.error('  读不到就必须判红 —— 放过去等于这个检查没生效。')
  failed = true
}

if (over.length > 0) {
  console.error('')
  console.error(`[FAIL] ${over.length} 张图超出 ${W_DOUBLE_MM} mm 上限（容差 ${TOLERANCE_MM} mm）：`)
  for (const r of over) {
    console.error(`    ${r.name}  ${r.widthMM.toFixed(1)} mm  （超出 ${(r.widthMM - W_DOUBLE_MM).toFixed(1)} mm）`)
  }
  console.error('')
  console.error('  修法：宽度随类别数增长的图要夹住上限 ——')
  console.error('      figsize=(min(W_DOUBLE, max(W_ONE_HALF, <按类别数的表达式>)), <高>)')
  console.error('  只夹下限（max(...)）会随类别数一直长下去。')
  failed = true
}

if (failed) process.exit(1)

console.log('')
console.log(`全部 ${rows.length} 张图都在 ${W_DOUBLE_MM} mm 以内`)
