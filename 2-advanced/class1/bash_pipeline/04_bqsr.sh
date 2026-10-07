#!/usr/bin/env bash
# ============================================================
# 2-4 BQSR — 塩基クオリティスコアの再校正
#         Base Quality Score Recalibration (GATK)
# 入力 / in : $DIR_MARKDUP/{sample}_marked_duplicates.bam
# 出力 / out: $DIR_BQSR/{sample}_bqsr.bam
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
source "$HERE/lib_common.sh"

log "=== 2-4 BaseRecalibrator 開始 / start ==="
require_cmd gatk
require_file "$REF"
require_file "$KNOWN1"
require_file "$KNOWN2"
require_file "$KNOWN3"
require_file "$KNOWN4"
mkdir -p "$DIR_BQSR"

for sample in "${SAMPLES[@]}"; do
  log "BQSR: $sample"

  # 入力 / inputs
  inbam="$DIR_MARKDUP/${sample}_marked_duplicates.bam"
  require_file "$inbam"

  # 出力 / outputs
  recaltable="$DIR_BQSR/${sample}_recalibration.table"
  outbam="$DIR_BQSR/${sample}_bqsr.bam"

  # 1) 再校正テーブルの作成 / build the recalibration table
  gatk --java-options "$JAVA_MEM" BaseRecalibrator \
    -I "$inbam" -R "$REF" \
    --known-sites "$KNOWN1" \
    --known-sites "$KNOWN2" \
    --known-sites "$KNOWN3" \
    --known-sites "$KNOWN4" \
    -O "$recaltable"

  # 2) 再校正の適用 / apply the recalibration
  gatk --java-options "$JAVA_MEM" ApplyBQSR \
    -R "$REF" -I "$inbam" --bqsr-recal-file "$recaltable" \
    -O "$outbam"
done

log "=== 2-4 BaseRecalibrator 完了 / done -> $DIR_BQSR ==="
ls -1 "$DIR_BQSR"
