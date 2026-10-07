#!/usr/bin/env bash
# ============================================================
# run_all.sh — パイプライン一括実行 / Run the whole pipeline
#   fastp -> BWA -> MarkDuplicates -> BQSR -> Mutect2 -> HaplotypeCaller
#
# 使い方 / Usage:
#   conda activate genome
#   bash run_all.sh                 # 全ステップ実行 / run all steps
#   bash run_all.sh 02_bwa 03_markduplicates   # 指定ステップのみ / only these
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
source "$HERE/lib_common.sh"

log "===== Genome Analysis パイプライン開始 / pipeline start ====="
log "COURSE_DIR = $COURSE_DIR"
log "WORK_DIR   = $WORK_DIR"
log "THREADS    = $THREADS   JAVA_MEM = $JAVA_MEM"

# 実行順 / execution order
steps=(
  01_fastp
  02_bwa
  03_markduplicates
  04_bqsr
  05_mutect2_somatic
  06_haplotypecaller_germline
)

# 引数があればそのステップのみ実行 / if args are given, run only those steps
if [[ $# -gt 0 ]]; then
  steps=("$@")
fi

for s in "${steps[@]}"; do
  script="$HERE/${s}.sh"
  [[ -f "$script" ]] || die "スクリプトが見つかりません / step script not found: $script"
  log ">>>>> $s"
  bash "$script"
done

log "===== 完了 / all done. 出力先 / outputs in: $WORK_DIR ====="
