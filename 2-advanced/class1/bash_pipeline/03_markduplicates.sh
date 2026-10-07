#!/usr/bin/env bash
# ============================================================
# 2-3 MarkDuplicates — PCR 重複リードのマーキング
#                    Mark PCR duplicate reads (Picard)
# 入力 / in : $DIR_BWA/{sample}.bam
# 出力 / out: $DIR_MARKDUP/{sample}_marked_duplicates.bam ( + .bai )
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
source "$HERE/lib_common.sh"

log "=== 2-3 MarkDuplicates 開始 / start ==="
require_cmd picard
mkdir -p "$DIR_MARKDUP"

for sample in "${SAMPLES[@]}"; do
  log "MarkDuplicates: $sample"

  # 入力 / inputs
  inbam="$DIR_BWA/${sample}.bam"
  require_file "$inbam"

  # 出力 / outputs
  outbam="$DIR_MARKDUP/${sample}_marked_duplicates.bam"
  metrics="$DIR_MARKDUP/${sample}_marked_dup_metrics.txt"
  bai="$DIR_MARKDUP/${sample}_marked_duplicates.bai"

  # 重複マーキング / mark duplicates
  picard "$JAVA_MEM" MarkDuplicates -I "$inbam" -O "$outbam" -M "$metrics"

  # 索引作成 / build BAM index
  picard "$JAVA_MEM" BuildBamIndex -INPUT "$outbam" -OUTPUT "$bai"
done

log "=== 2-3 MarkDuplicates 完了 / done -> $DIR_MARKDUP ==="
ls -1 "$DIR_MARKDUP"
