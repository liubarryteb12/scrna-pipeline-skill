#!/usr/bin/env python3
"""
gse254707_tar_to_h5ad.py — GSE254707 RAW tar（4988 个 单细胞×基因 表）→ 单个 h5ad

GSE254707（银屑病单细胞，SMART-seq 孔板分型）的 GEO 补充包是 4988 个
`GSM*_HLA*.txt.gz`，每个文件 = 一个细胞（列2：gene symbol \t UMI count，
~33218 行）。这不是 10x 格式，skill 仓 00_fetch.py 的三种 kind 都吃不下，
本转换器把它合成一张 稀疏矩阵 h5ad：

  X        = 原始 UMI 计数（uint32 稀疏 CSC，行为基因、列为细胞 → 转置后 obs=细胞）
  obs      = donor / region / diagnosis / group / celltype / class / gsm / title
             （来自 GSE254707_family.soft.gz 的 !Sample_* 行，随仓内嵌映射表）
  var      = gene symbol（重复名 make_unique 后保留原名列）

分组语义（与甲方「分组情况.txt」对齐）：
  group = LS（Psoriasis_PP_*，皮损） / NLS（Psoriasis_PN_*，非皮损） / HC（Healthy_H_*，健康）
  —— region 字段（PP/PN/PP?）是原始区域能级，group 是三值分析分组。

用法（容器内）：
  python3 scripts/gse254707_tar_to_h5ad.py --tar <RAW.tar> --out raw.h5ad \
      [--soft GSE254707_family.soft.gz] [--max-files N]
--soft 缺省时用脚本内嵌的 (title 前缀 → group) 规则 + GSM 序号顺序。
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
from common import log_error, log_info  # noqa: E402

# title 前缀 → 分析分组（甲方分组语义）
GROUP_RULES = [
    (re.compile(r'^Psoriasis_PP_', re.I), 'LS'),
    (re.compile(r'^Psoriasis_PN_', re.I), 'NLS'),
    (re.compile(r'^Healthy_', re.I), 'HC'),
]


def parse_soft(soft_path: Path) -> dict:
    """GSE family soft → {GSM: {title, donor, region, diagnosis, celltype, class}}"""
    out: dict[str, dict] = {}
    cur = None
    op = gzip.open if str(soft_path).endswith('.gz') else open
    with op(soft_path, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('^SAMPLE'):
                cur = line.split('=', 1)[1].strip()
                out[cur] = {}
            elif cur and line.startswith('!Sample_'):
                key, _, val = line.partition('=')
                key = key.strip().removeprefix('!Sample_').strip()
                val = val.strip()
                if key == 'title':
                    out[cur]['title'] = val
                elif key == 'characteristics_ch1':
                    k, _, v = val.partition(':')
                    k = k.strip().lower()
                    v = v.strip()
                    if k in ('donor', 'region', 'diagnosis', 'celltype', 'class'):
                        out[cur][k] = v
    return out


def group_of(meta: dict) -> str:
    title = meta.get('title', '')
    for pat, g in GROUP_RULES:
        if pat.match(title):
            return g
    region = meta.get('region', '')
    if region.upper() in ('PP', 'LS'):
        return 'LS'
    if region.upper() in ('PN', 'NLS'):
        return 'NLS'
    if meta.get('diagnosis', '').lower() in ('healthy', 'healthy control'):
        return 'HC'
    return 'NA'


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tar', required=True, help='GSE254707_RAW.tar 路径')
    ap.add_argument('--out', required=True, help='输出 h5ad 路径')
    ap.add_argument('--soft', default='', help='GSE254707_family.soft.gz（可选，给 obs 元数据）')
    ap.add_argument('--max-files', type=int, default=0, help='调试：只读前 N 个文件')
    a = ap.parse_args()

    import anndata as ad
    import pandas as pd

    tar_path = Path(a.tar)
    out_path = Path(a.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    soft = parse_soft(Path(a.soft)) if a.soft else {}

    names = [m for m in tarfile.open(tar_path, 'r').getnames()
             if m.endswith('.txt.gz')]
    names.sort()
    if a.max_files:
        names = names[:a.max_files]
    n_files = len(names)
    if n_files == 0:
        log_error('tar 里没有 .txt.gz 文件')
        return 1

    gene_index: dict[str, int] = {}
    rows, cols, vals = [], [], []   # COO（行=基因，列=细胞）
    obs_rows = []
    ncols = 0
    skipped_dup_gene = 0

    with tarfile.open(tar_path, 'r') as tf:
        for ci, name in enumerate(names):
            gsm = Path(name).name.split('_')[0]
            meta = dict(soft.get(gsm, {}))
            meta.setdefault('title', Path(name).name.replace('.txt.gz', ''))
            meta['gsm'] = gsm
            meta['group'] = group_of(meta)
            fh = tf.extractfile(name)
            with gzip.open(fh, 'rt', encoding='utf-8', errors='replace') as fh2:
                for line in fh2:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 2:
                        continue
                    g = parts[0]
                    try:
                        v = float(parts[1])
                    except ValueError:
                        continue
                    if v == 0:
                        continue
                    gi = gene_index.setdefault(g, len(gene_index))
                    rows.append(gi)
                    cols.append(ci)
                    vals.append(v)
            obs_rows.append(meta)
            ncols += 1
            if (ci + 1) % 500 == 0:
                log_info(f'  {ci + 1}/{n_files} 细胞...')

    n_genes = len(gene_index)
    log_info(f'合并完成: {n_genes} 基因 x {ncols} 细胞（非零 {len(vals)}）')

    X = sp.coo_matrix((vals, (rows, cols)), shape=(n_genes, ncols)).tocsc()
    # 计数必须是非负整数（SMART-seq UMI）——非整数值直接拒绝
    data = X.data
    if not np.allclose(data, np.round(data), atol=1e-8):
        bad = float(np.mean(~np.isclose(data, np.round(data), atol=1e-8)))
        log_error(f'**非整数计数占比 {bad:.2%} —— 不是 UMI 计数矩阵，拒绝产出**')
        return 2
    X_int = X.astype(np.int32)

    genes = sorted(gene_index, key=gene_index.get)
    var = pd.DataFrame({'gene_symbol': genes}, index=genes)
    obs = pd.DataFrame(obs_rows,
                       index=[r.get('title', f'cell_{i}') for i, r in enumerate(obs_rows)])
    obs.index.name = 'cell_id'
    # pseudobulk 样本单位 = (donor, group)：银屑病供体每人同时有 LS 与 NLS 细胞
    #（配对设计），sample_key=donor 会判 mixed_group；正确样本 = donor_region。
    obs['sample_id'] = obs['donor'].astype(str) + '_' + obs['group'].astype(str)

    adata = ad.AnnData(X=X_int.T.tocsr(), obs=obs, var=var)   # 转置：obs=细胞
    adata.var_names_make_unique()
    adata.raw = None
    adata.write_h5ad(out_path)
    log_info(f'已写出 {out_path}: {adata.n_obs} 细胞 x {adata.n_vars} 基因')
    grp = obs['group'].value_counts().to_dict()
    log_info(f'分组分布: {grp}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
