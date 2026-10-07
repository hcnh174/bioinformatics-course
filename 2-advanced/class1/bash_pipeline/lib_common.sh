#!/usr/bin/env bash
# ============================================================
# lib_common.sh — 共通関数 / Shared helper functions
#   各ステップスクリプトから source して使います。
#   Sourced by every step script.
# ============================================================

# タイムスタンプ付きログ / timestamped logging
log()  { printf '[%s] %s\n'        "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
warn() { printf '[%s] [WARN] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
die()  { printf '[%s] [ERROR] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; exit 1; }

# コマンドの存在確認 / ensure an external command is available
require_cmd() {
  command -v "$1" >/dev/null 2>&1 \
    || die "コマンドが見つかりません / command not found: $1  ( 'conda activate genome' を確認してください )"
}

# 入力ファイルの存在確認 / ensure an input file exists
require_file() {
  [[ -f "$1" ]] || die "入力ファイルが見つかりません / input file not found: $1"
}
