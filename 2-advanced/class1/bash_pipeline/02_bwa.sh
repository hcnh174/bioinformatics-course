#!/usr/bin/env bash
# ============================================================
# 2-2 BWA — 参照ゲノムへのアライメント + 整列・索引付け
#         Read alignment to the reference + sort + index
# 入力 / in : $DIR_FASTP/{sample}_{1,2}_trimmed.fastq.gz
# 出力 / out: $DIR_BWA/{sample}.bam ( + .bai )
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
source "$HERE/lib_common.sh"

log "=== 2-2 BWA 開始 / start ==="
require_cmd bwa
require_cmd samtools
require_file "${BWA_INDEX}.bwt"   # インデックスの存在確認 / check the BWA index prefix
mkdir -p "$DIR_BWA"

for sample in "${SAMPLES[@]}"; do
  log "bwa mem: $sample"

  # 入力 / inputs
  in1="$DIR_FASTP/${sample}_1_trimmed.fastq.gz"
  in2="$DIR_FASTP/${sample}_2_trimmed.fastq.gz"
  require_file "$in1"
  require_file "$in2"

  # 出力 / outputs
  sam="$DIR_BWA/${sample}.sam"
  tmpbam="$DIR_BWA/${sample}_tmp.bam"
  outbam="$DIR_BWA/${sample}.bam"

  # アライメント / align  ( -R でリードグループを付与 / add read group )
  bwa mem -t "$THREADS" \
    -R "@RG\tID:${sample}\tPL:ILLUMINA\tSM:${sample}" \
    "$BWA_INDEX" "$in1" "$in2" > "$sam"

  # SAM -> BAM 変換・整列・索引付け / convert, sort, index
  samtools view -bS "$sam" > "$tmpbam"
  samtools sort -@ "$THREADS" "$tmpbam" -o "$outbam"
  samtools index "$outbam"

  # 中間ファイルを削除 / clean up intermediates ( SAM は確認用に残す / keep SAM for inspection )
  rm -f "$tmpbam"
done

log "=== 2-2 BWA 完了 / done -> $DIR_BWA ==="
ls -1 "$DIR_BWA"
