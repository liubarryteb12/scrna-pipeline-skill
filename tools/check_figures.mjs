#!/usr/bin/env node
/**
 * tools/check_figures.mjs — 像素级检查图不是空白
 *
 * **为什么需要它。** 一个"成功生成"的 PNG 完全可能是空白的：坐标轴画了
 * 但数据点没画、图例占了整个画布、或者数据全被过滤掉导致散点为空。
 * 这类问题在 CI 里表现为绿色 job + 一个 0 字节以外的正常文件，
 * 只有打开看才发现是白图。
 *
 * 做法：解 PNG（zlib inflate + 逐行反滤波），统计非背景像素占比。
 * 纯 JS，不依赖任何图形库。
 *
 * 判定：非背景像素占比 < MIN_INK 判为空白。
 * 阈值不能设太高 —— 一个只有坐标轴的图也有少量墨迹，那正是要抓的。
 *
 * 用法: node tools/check_figures.mjs <目录> [更多目录...]
 */

import { readdirSync, readFileSync, existsSync, statSync } from "node:fs";
import { join } from "node:path";
import { inflateSync } from "node:zlib";

// 非背景像素占比低于此值判为空白
const MIN_INK = 0.002;

function decodePng(buf) {
  if (buf.readUInt32BE(0) !== 0x89504e47) throw new Error("不是 PNG");
  let pos = 8;
  let width = 0, height = 0, bitDepth = 0, colorType = 0;
  const idat = [];
  while (pos < buf.length) {
    const len = buf.readUInt32BE(pos);
    const type = buf.toString("ascii", pos + 4, pos + 8);
    const data = buf.subarray(pos + 8, pos + 8 + len);
    if (type === "IHDR") {
      width = data.readUInt32BE(0);
      height = data.readUInt32BE(4);
      bitDepth = data[8];
      colorType = data[9];
    } else if (type === "IDAT") {
      idat.push(data);
    } else if (type === "IEND") {
      break;
    }
    pos += 12 + len;
  }
  if (bitDepth !== 8) throw new Error(`只支持 8 位深度，实际 ${bitDepth}`);
  const channels = { 0: 1, 2: 3, 3: 1, 4: 2, 6: 4 }[colorType];
  if (!channels) throw new Error(`不支持的 colorType ${colorType}`);

  const raw = inflateSync(Buffer.concat(idat));
  const stride = width * channels;
  const out = Buffer.alloc(height * stride);

  // 逐行反滤波（PNG 的 5 种滤波器）
  for (let y = 0; y < height; y++) {
    const ft = raw[y * (stride + 1)];
    const src = raw.subarray(y * (stride + 1) + 1, y * (stride + 1) + 1 + stride);
    const dst = out.subarray(y * stride, (y + 1) * stride);
    const prev = y > 0 ? out.subarray((y - 1) * stride, y * stride) : null;
    for (let x = 0; x < stride; x++) {
      const a = x >= channels ? dst[x - channels] : 0;
      const b = prev ? prev[x] : 0;
      const c = prev && x >= channels ? prev[x - channels] : 0;
      let v = src[x];
      if (ft === 1) v += a;
      else if (ft === 2) v += b;
      else if (ft === 3) v += (a + b) >> 1;
      else if (ft === 4) {
        const p = a + b - c;
        const pa = Math.abs(p - a), pb = Math.abs(p - b), pc = Math.abs(p - c);
        v += pa <= pb && pa <= pc ? a : pb <= pc ? b : c;
      }
      dst[x] = v & 0xff;
    }
  }
  return { width, height, channels, pixels: out };
}

function inkFraction(png) {
  const { width, height, channels, pixels } = png;
  // 取四角的多数色作为背景色（图通常白底，但不要假定）
  const corner = (x, y) => {
    const o = (y * width + x) * channels;
    return [pixels[o], pixels[o + 1], pixels[o + 2]];
  };
  const cs = [corner(0, 0), corner(width - 1, 0), corner(0, height - 1),
              corner(width - 1, height - 1)];
  const bg = [0, 1, 2].map((i) => {
    const v = cs.map((c) => c[i]).sort((a, b) => a - b);
    return v[Math.floor(v.length / 2)];
  });

  let ink = 0;
  const total = width * height;
  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      const o = (y * width + x) * channels;
      const d = Math.abs(pixels[o] - bg[0]) + Math.abs(pixels[o + 1] - bg[1]) +
                Math.abs(pixels[o + 2] - bg[2]);
      if (d > 30) ink++;
    }
  }
  return ink / total;
}

const dirs = process.argv.slice(2);
if (dirs.length === 0) {
  console.error("用法: node tools/check_figures.mjs <目录> [更多目录...]");
  process.exit(1);
}

let n = 0, blank = 0;
for (const dir of dirs) {
  if (!existsSync(dir)) continue;
  for (const e of readdirSync(dir, { withFileTypes: true })) {
    if (!e.isFile() || !e.name.endsWith(".png")) continue;
    const p = join(dir, e.name);
    n++;
    const size = statSync(p).size;
    try {
      const png = decodePng(readFileSync(p));
      const frac = inkFraction(png);
      if (frac < MIN_INK) {
        console.error(`  [空白] ${p}  ${png.width}x${png.height}  ` +
                      `墨迹占比 ${(frac * 100).toFixed(4)}%  (${size} B)`);
        blank++;
      } else {
        console.log(`  [OK]   ${e.name}  ${png.width}x${png.height}  ` +
                    `墨迹 ${(frac * 100).toFixed(2)}%`);
      }
    } catch (err) {
      console.error(`  [错误] ${p}: ${err.message}`);
      blank++;
    }
  }
}

if (n === 0) {
  console.error("没有找到任何 PNG —— 图目录不对或图没生成");
  process.exit(1);
}
if (blank > 0) {
  console.error(`\n${blank}/${n} 张图有问题`);
  process.exit(1);
}
console.log(`\n${n} 张图全部非空白`);
