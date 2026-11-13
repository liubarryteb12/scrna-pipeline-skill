#!/usr/bin/env bash
set -euo pipefail

checkpoint_name="checkpoint_scrna_bridge"
echo "step_executed=${checkpoint_name}"

export LC_ALL=C.UTF-8
export OMP_NUM_THREADS=3
export OPENBLAS_NUM_THREADS=3
export MKL_NUM_THREADS=3

input_path=""
seed_override=""
from_lang=""
to_lang=""
stage="post-dimred"

while [ "$#" -gt 0 ]; do
  case "$1" in
    --input)
      if [ "$#" -lt 2 ] || [ -z "$2" ]; then
        echo "[BRIDGE FAIL] 参数 --input 缺少路径" >&2
        exit 2
      fi
      input_path="$2"
      shift 2
      ;;
    --from)
      if [ "$#" -lt 2 ] || [ -z "$2" ]; then
        echo "[BRIDGE FAIL] 参数 --from 缺少语言" >&2
        exit 2
      fi
      from_lang="$2"
      shift 2
      ;;
    --to)
      if [ "$#" -lt 2 ] || [ -z "$2" ]; then
        echo "[BRIDGE FAIL] 参数 --to 缺少语言" >&2
        exit 2
      fi
      to_lang="$2"
      shift 2
      ;;
    --stage)
      if [ "$#" -lt 2 ] || [ -z "$2" ]; then
        echo "[BRIDGE FAIL] 参数 --stage 缺少阶段" >&2
        exit 2
      fi
      stage="$2"
      shift 2
      ;;
    --seed)
      if [ "$#" -lt 2 ] || [ -z "$2" ]; then
        echo "[BRIDGE FAIL] 参数 --seed 缺少整数值" >&2
        exit 2
      fi
      seed_override="$2"
      shift 2
      ;;
    *)
      echo "[BRIDGE FAIL] 未知参数: $1" >&2
      exit 2
      ;;
  esac
done

if [ -z "$input_path" ]; then
  echo "[BRIDGE FAIL] 缺少必填参数 --input" >&2
  if [ -z "$from_lang" ] || [ -z "$to_lang" ]; then
    echo "[BRIDGE FAIL] 缺少必填参数 --from 与 --to（语言方向必须显式声明）" >&2
    exit 2
  fi
  exit 2
fi

case "$from_lang" in
  Python|R) ;;
  *)
    echo "[BRIDGE FAIL] 未知 --from: $from_lang（仅 Python / R）" >&2
    exit 2
    ;;
esac

case "$to_lang" in
  Python|R) ;;
  *)
    echo "[BRIDGE FAIL] 未知 --to: $to_lang（仅 Python / R）" >&2
    exit 2
    ;;
esac

case "$stage" in
  post-dimred|pre-dimred) ;;
  *)
    echo "[BRIDGE FAIL] 未知 --stage: $stage（仅 post-dimred / pre-dimred）" >&2
    exit 2
    ;;
esac

seed_value="${seed_override:-${BIOANALYS_SEED:-}}"
if [ -z "$seed_value" ]; then
  echo "BIOANALYS_SEED 未注入（引擎契约：seed 必须注入）" >&2
  exit 3
fi
if [[ ! "$seed_value" =~ ^-?[0-9]+$ ]]; then
  echo "BIOANALYS_SEED 不是整数: $seed_value" >&2
  exit 3
fi

if [ ! -f "$input_path" ]; then
  echo "[BRIDGE FAIL] 输入文件不存在: $input_path" >&2
  exit 2
fi

case "$input_path" in
  *.h5ad) ;;
  *)
    echo "[BRIDGE SKIP] $checkpoint_name 非 h5ad 输入，无单细胞桥接断言可做: $input_path"
    exit 0
    ;;
esac

engine_path="${BRIDGE_ENGINE:-/mnt/ssd/bio_data/ops/verify_bridge_checkpoint.py}"
if [ ! -f "$engine_path" ] || [ ! -r "$engine_path" ] || ! command -v python3 >/dev/null 2>&1; then
  echo "[BRIDGE SKIP] $checkpoint_name 桥接引擎不可达: $engine_path"
  exit 0
fi

CHECKPOINT_INPUT="$input_path" CHECKPOINT_STAGE="$stage" python3 - <<'PYSCRIPT'
import os
import re
import sys

import anndata as ad
import numpy as np

ENSEMBL_PATTERN = re.compile(r"^ENS[A-Z]*G[0-9]+(?:\.[0-9]+)?$", re.IGNORECASE)
GSM_PATTERN = re.compile(r"^GSM[0-9]+$", re.IGNORECASE)
NAMED_BARCODE_PATTERN = re.compile(r"^barcode(?:[_-].+)?$", re.IGNORECASE)
TENX_BARCODE_PATTERN = re.compile(
    r"^(?:.+,)?[ACGTN]{8,}(?:-[0-9]+)?$", re.IGNORECASE
)
CELL_NAME_PATTERN = re.compile(r"^cell[_-].+$", re.IGNORECASE)


def fail(message):
    sys.stderr.write("[BRIDGE FAIL] line2_scrna %s\n" % message)
    sys.stderr.flush()
    sys.exit(1)


def warn(message):
    print("[BRIDGE WARN] line2_scrna %s" % message)


def is_cell_identifier(value):
    text = str(value).strip()
    return bool(
        GSM_PATTERN.fullmatch(text)
        or NAMED_BARCODE_PATTERN.fullmatch(text)
        or TENX_BARCODE_PATTERN.fullmatch(text)
        or CELL_NAME_PATTERN.fullmatch(text)
    )


input_path = os.environ["CHECKPOINT_INPUT"]
stage = os.environ["CHECKPOINT_STAGE"]

try:
    adata = ad.read_h5ad(input_path)
except Exception as exc:
    fail("h5ad 无法读取（%s: %s）" % (type(exc).__name__, exc))

if adata.n_obs == 0 or adata.n_vars == 0:
    fail("矩阵维度为 0: obs=%d var=%d" % (adata.n_obs, adata.n_vars))

var_names = [str(value) for value in adata.var_names]
obs_names = [str(value) for value in adata.obs_names]
if not var_names:
    fail("var_names 为空")
if not obs_names:
    fail("obs_names 为空")

ensembl_count = sum(1 for value in var_names if ENSEMBL_PATTERN.fullmatch(value))
if ensembl_count > len(var_names) * 0.5:
    fail(
        "var_names 超过半数是 Ensembl ID（%d/%d），未映射为基因 Symbol"
        % (ensembl_count, len(var_names))
    )

cell_like_obs = sum(1 for value in obs_names if is_cell_identifier(value))
if cell_like_obs <= len(obs_names) * 0.5:
    fail(
        "obs_names 不足半数是 barcode/GSM 细胞标识（%d/%d），疑似矩阵方向反转"
        % (cell_like_obs, len(obs_names))
    )

cell_like_var = sum(1 for value in var_names if is_cell_identifier(value))
if cell_like_var > len(var_names) * 0.5:
    fail(
        "var_names 超过半数是 barcode/GSM（%d/%d），疑似矩阵方向反转为基因 x 细胞"
        % (cell_like_var, len(var_names))
    )

print(
    "[BRIDGE OK] line2_scrna gene symbols: %d, Ensembl=%d"
    % (len(var_names), ensembl_count)
)
print(
    "[BRIDGE OK] line2_scrna matrix direction: cells=%d x genes=%d"
    % (adata.n_obs, adata.n_vars)
)
print(
    "[BRIDGE OK] line2_scrna cell identifiers: %d/%d"
    % (cell_like_obs, len(obs_names))
)

if "X_pca" not in adata.obsm:
    if stage == "pre-dimred":
        warn("obsm 无 X_pca（stage=pre-dimred，早期矩阵降级为告警）")
    else:
        fail("obsm 无 X_pca，降维嵌入缺失")
else:
    pca = np.asarray(adata.obsm["X_pca"])
    if pca.ndim != 2:
        fail("obsm['X_pca'] 不是二维矩阵，实际 ndim=%d" % pca.ndim)
    if pca.shape[0] != adata.n_obs:
        fail(
            "obsm['X_pca'] 行数与 obs 不一致: %d vs %d"
            % (pca.shape[0], adata.n_obs)
        )
    print(
        "[BRIDGE OK] line2_scrna X_pca: shape=%d x %d"
        % (pca.shape[0], pca.shape[1])
    )
PYSCRIPT

set +e
engine_output="$(
  python3 "$engine_path" \
    --input "$input_path" \
    --from "$from_lang" \
    --to "$to_lang" \
    --stage "$stage" 2>&1
)"
engine_rc=$?
set -e

if [ -n "$engine_output" ]; then
  printf '%s\n' "$engine_output"
fi

if [ "$engine_rc" -eq 0 ]; then
  echo "[BRIDGE OK] $checkpoint_name 引擎四项断言通过（$from_lang -> $to_lang, stage=$stage）"
  exit 0
fi

echo "[BRIDGE FAIL] $checkpoint_name 引擎校验失败（退出码 $engine_rc）" >&2
exit "$engine_rc"
