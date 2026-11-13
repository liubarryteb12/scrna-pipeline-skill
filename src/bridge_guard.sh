#!/usr/bin/env bash
# bridge_guard.sh：跨语言桥接四项断言的本机守卫。
#
# 用途：母版脚本在 Python 与 R 之间传 .h5ad / .rds 时，先调用本守卫，
# 由 verify_bridge_checkpoint.py 校验基因名、矩阵方向、细胞名/样本名、
# embedding 四项断言。断言不过，母版中止。
#
# 引擎解析顺序：
#   1. BRIDGE_ENGINE 环境变量指向的 verify_bridge_checkpoint.py
#   2. 与本脚本同目录的 verify_bridge_checkpoint.py
#   3. PATH 上的 verify_bridge_checkpoint.py
#
# 找不到引擎时打印 [BRIDGE SKIP] 并放行，避免离线干跑被基础设施绑架；
# 设 BIO_STRICT_BRIDGE=1 时改为 [BRIDGE FAIL] 硬阻断。
#
# 用法：
#   bridge_guard.sh --input <路径> --from Python|R --to R|Python \
#                   [--script <母版名>] [--stage post-dimred|pre-dimred]
#
# 退出码：0 通过或跳过；非 0 为断言失败或引擎异常。
set -euo pipefail

where="bridge_guard"
input_path=""
from_lang=""
to_lang=""
script_name=""
stage="post-dimred"

while [ "$#" -gt 0 ]; do
  case "$1" in
    --input)
      [ "$#" -ge 2 ] && [ -n "$2" ] || { echo "[BRIDGE FAIL] $where --input 缺少路径" >&2; exit 2; }
      input_path="$2"; shift 2 ;;
    --from)
      [ "$#" -ge 2 ] && [ -n "$2" ] || { echo "[BRIDGE FAIL] $where --from 缺少语言" >&2; exit 2; }
      from_lang="$2"; shift 2 ;;
    --to)
      [ "$#" -ge 2 ] && [ -n "$2" ] || { echo "[BRIDGE FAIL] $where --to 缺少语言" >&2; exit 2; }
      to_lang="$2"; shift 2 ;;
    --script)
      [ "$#" -ge 2 ] && [ -n "$2" ] || { echo "[BRIDGE FAIL] $where --script 缺少名称" >&2; exit 2; }
      script_name="$2"; shift 2 ;;
    --stage)
      [ "$#" -ge 2 ] && [ -n "$2" ] || { echo "[BRIDGE FAIL] $where --stage 缺少取值" >&2; exit 2; }
      stage="$2"; shift 2 ;;
    -h|--help)
      sed -n '2,20p' "$0"; exit 0 ;;
    *)
      echo "[BRIDGE FAIL] $where 未知参数: $1" >&2; exit 2 ;;
  esac
done

case "$from_lang" in Python|R) ;; *) echo "[BRIDGE FAIL] $where --from 只接受 Python 或 R" >&2; exit 2 ;; esac
case "$to_lang" in Python|R) ;; *) echo "[BRIDGE FAIL] $where --to 只接受 Python 或 R" >&2; exit 2 ;; esac
case "$stage" in post-dimred|pre-dimred) ;; *) echo "[BRIDGE FAIL] $where --stage 只接受 post-dimred 或 pre-dimred" >&2; exit 2 ;; esac
[ -n "$input_path" ] || { echo "[BRIDGE FAIL] $where 缺少 --input" >&2; exit 2; }
[ -f "$input_path" ] || { echo "[BRIDGE FAIL] $where 输入不存在: $input_path" >&2; exit 2; }

case "$input_path" in
  *.h5ad|*.rds) ;;
  *) echo "[BRIDGE SKIP] $where 输入不是 .h5ad/.rds，无四项断言对象: $input_path"; exit 0 ;;
esac

engine="${BRIDGE_ENGINE:-}"
if [ -z "$engine" ] || [ ! -f "$engine" ]; then
  for candidate in "$(dirname "$0")/verify_bridge_checkpoint.py" "$(command -v verify_bridge_checkpoint.py 2>/dev/null || true)"; do
    if [ -n "$candidate" ] && [ -f "$candidate" ]; then engine="$candidate"; break; fi
  done
fi

if [ -z "$engine" ] || [ ! -f "$engine" ]; then
  if [ "${BIO_STRICT_BRIDGE:-0}" = "1" ]; then
    echo "[BRIDGE FAIL] $where 未找到 verify_bridge_checkpoint.py（BIO_STRICT_BRIDGE=1）" >&2
    exit 2
  fi
  echo "[BRIDGE SKIP] $where 未找到 verify_bridge_checkpoint.py，设 BRIDGE_ENGINE 指向引擎"
  exit 0
fi

label="${script_name:-$where}"
if output=$(python3 "$engine" --input "$input_path" --from "$from_lang" --to "$to_lang" --stage "$stage" 2>&1); then
  echo "[BRIDGE OK] $label 四项断言通过"
  exit 0
else
  rc=$?
  printf '%s\n' "$output" >&2
  echo "[BRIDGE FAIL] $label 四项断言失败（退出码 $rc）" >&2
  exit "$rc"
fi
