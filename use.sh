#!/usr/bin/env bash
# 从**云端仓库**按需拉取本 skill。
#
# 本 skill 的唯一真源是 GitHub 仓库，不是任何一台机器上的目录。
# 这个脚本的职责就是「要用的时候把它从云端拉下来」——
# 它不复制、不软链你的开发目录，只从远端 clone/fetch。
#
#   ./use.sh                    拉取/更新到缓存目录，打印路径
#   ./use.sh --register         拉取后注册进 skill root（默认 ~/.agents/skills 等）
#   ./use.sh --dir DIR          指定缓存目录（默认 $DSH_SKILL_CACHE 或 ~/.cache/dsh-skills）
#   ./use.sh --ref REF          钉住分支/标签/commit（默认 main）
#   ./use.sh --print-path       只打印路径，便于 `SKILL=$(./use.sh --print-path)`
#   ./use.sh --clean            清掉缓存副本，并撤销注册
#   ./use.sh --unregister       只从 skill root 移除注册，保留缓存
#
# 为什么默认不注册：注册会让本机多一份状态，换机器就要重建。
# 需要在本机被 agent 直接发现时才加 --register。
#
# 关于 set -e：本脚本仍然开着它，但**不依赖它**。实测（bash 5.3，Git for
# Windows）在 `resolved="$(pull)"` 这种「函数在命令替换里」的结构下，
# 函数内部的失败不一定中止外层。所以每个可能失败的步骤都显式 `|| die`，
# 调用点也显式检查返回值 —— 出错的路径必须自己说出来，不能指望 shell 规则。
set -euo pipefail

SKILL_NAME="scrna-pipeline-skill"
REPO_URL="https://github.com/liubarryteb12/scrna-pipeline-skill.git"
DEFAULT_REF="main"

# ---------------------------------------------------------------------------
# 函数先定义，参数解析后执行 —— 顺序反了会把函数调用当成未定义命令。
# ---------------------------------------------------------------------------

die() { printf '错误: %s\n' "$*" >&2; exit 1; }

# 打印文件顶部的注释块作为帮助。
# **不用 sed / awk**：实测在 PATH 不完整的 shell 里（Git bash 被直接调用）
# 它们会 `command not found`，而 help 恰恰是出问题时最需要能跑的那条路径。
print_help() {
  local line body
  while IFS= read -r line; do
    case "$line" in
      '#!'*) continue ;;                 # shebang 不打印
      '#'*)  body="${line#\#}"; body="${body# }"; printf '%s\n' "$body" ;;
      *)     break ;;                    # 第一个非注释行 -> 头注释结束
    esac
  done < "$0"
}

# 本机所有会被 agent 扫描的 skill root
skill_roots() {
  printf '%s\n' "$HOME/.agents/skills" "$HOME/.claude/skills" "$HOME/.codex/skills"
}

# 把工作副本切到 REF。分支和标签要分开处理：
# `git checkout main` 在游离 HEAD 上会失败，`-B` 则能重建分支。
checkout_ref() {
  local repo="$1" ref="$2"
  if git -C "$repo" rev-parse --verify --quiet "origin/$ref" >/dev/null 2>&1; then
    git -C "$repo" checkout --quiet -B "$ref" "origin/$ref" \
      || die "切到分支 $ref 失败"
  else
    git -C "$repo" checkout --quiet --detach "$ref" \
      || die "切到 $ref 失败（既不是 origin 上的分支，也不是有效的 tag/commit）"
  fi
}

# 拉取/更新。成功时把结果路径打到 stdout，过程信息一律 stderr。
pull() {
  if [[ -d "$DEST/.git" ]]; then
    echo "更新 $DEST（ref=$REF）" >&2
    git -C "$DEST" fetch --quiet --tags --prune origin \
      || die "git fetch 失败：$REPO_URL"
    checkout_ref "$DEST" "$REF"
  else
    echo "克隆 $REPO_URL -> $DEST（ref=$REF）" >&2
    mkdir -p "$CACHE_ROOT" || die "无法创建缓存目录 $CACHE_ROOT"
    git clone --quiet "$REPO_URL" "$DEST" \
      || die "git clone 失败：$REPO_URL（仓库不存在 / 无网络 / 需要代理？）"
    checkout_ref "$DEST" "$REF"
  fi

  # 拉下来的东西必须真的是这个 skill，不能是个空壳。
  # **这道校验是必需的，不是保险**：实测 git clone 失败时后面的步骤照样跑了，
  # 最后就是靠这里拦住的。远端改名/换结构也会安静地给一个空目录。
  [[ -f "$DEST/SKILL.md" ]] \
    || die "拉取后 $DEST/SKILL.md 不存在 —— 仓库结构不对"

  # 记下实际拿到的 commit，便于回答「我用的到底是哪一版」
  echo "  commit: $(git -C "$DEST" rev-parse --short HEAD)" >&2
  printf '%s\n' "$DEST"
}

register() {
  local dest="$1" root link
  while IFS= read -r root; do
    mkdir -p "$root" || die "无法创建 skill root $root"
    link="$root/$SKILL_NAME"
    if [[ -e "$link" || -L "$link" ]]; then
      rm -rf "$link" || die "无法移除旧注册 $link"
    fi
    ln -s "$dest" "$link" || die "无法创建软链 $link"
    echo "  已注册 $link -> $dest" >&2
  done < <(skill_roots)
}

unregister() {
  local root link
  while IFS= read -r root; do
    link="$root/$SKILL_NAME"
    if [[ -e "$link" || -L "$link" ]]; then
      rm -rf "$link" || die "无法移除 $link"
      echo "  已移除 $link" >&2
    fi
  done < <(skill_roots)
}

# ---------------------------------------------------------------------------
# 参数解析
# ---------------------------------------------------------------------------
CACHE_ROOT="${DSH_SKILL_CACHE:-$HOME/.cache/dsh-skills}"
DEST="$CACHE_ROOT/$SKILL_NAME"
REF="$DEFAULT_REF"
MODE="pull"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --register)   MODE="register";   shift ;;
    --unregister) MODE="unregister"; shift ;;
    --clean)      MODE="clean";      shift ;;
    --print-path) MODE="print-path"; shift ;;
    --dir)        [[ $# -ge 2 ]] || die "--dir 需要参数"; CACHE_ROOT="$2"; DEST="$CACHE_ROOT/$SKILL_NAME"; shift 2 ;;
    --ref)        [[ $# -ge 2 ]] || die "--ref 需要参数"; REF="$2"; shift 2 ;;
    -h|--help)    print_help; exit 0 ;;
    *) die "未知参数: $1（--help 看用法）" ;;
  esac
done

# ---------------------------------------------------------------------------
# 执行
# ---------------------------------------------------------------------------
case "$MODE" in
  clean)
    if [[ -d "$DEST" ]]; then
      rm -rf "$DEST" || die "无法删除 $DEST"
      echo "已删除 $DEST" >&2
    fi
    unregister
    ;;
  unregister)
    unregister
    ;;
  print-path)
    [[ -f "$DEST/SKILL.md" ]] || die "尚未拉取（$DEST 不存在）。先跑 ./use.sh"
    printf '%s\n' "$DEST"
    ;;
  register)
    # 显式检查返回值，不依赖 set -e 在命令替换里的行为
    if ! resolved="$(pull)"; then
      die "拉取失败，见上面的错误信息"
    fi
    echo "已拉取: $resolved" >&2
    register "$resolved"
    echo >&2
    echo "完成。本机现在能发现这个 skill 了；换机器重跑 ./use.sh --register 即可。" >&2
    ;;
  pull)
    if ! resolved="$(pull)"; then
      die "拉取失败，见上面的错误信息"
    fi
    echo >&2
    echo "skill 已在: $resolved" >&2
    echo "要让本机 agent 直接发现它，加 --register；只当资料读则不用。" >&2
    ;;
esac
