#!/usr/bin/env bash
# ============================================================
# 2-5 Mutect2 — 体細胞変異の検出（腫瘍 vs 正常）
#            Somatic variant calling (tumor vs normal, GATK)
# 入力 / in : $DIR_BQSR/{tumor,normal}_bqsr.bam
# 出力 / out: $DIR_MUTECT2/tumor_somatic_filtered.vcf.gz
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
source "$HERE/lib_common.sh"

log "=== 2-5 Mutect2 (somatic) 開始 / start ==="
require_cmd gatk
require_file "$REF"
require_file "$GERMLINE_RESOURCE"
require_file "$INTERVALS"
mkdir -p "$DIR_MUTECT2"

# 入力 BAM / input BAMs
tumorbam="$DIR_BQSR/${TUMOR}_bqsr.bam"
normalbam="$DIR_BQSR/${NORMAL}_bqsr.bam"
require_file "$tumorbam"
require_file "$normalbam"

# 出力 / outputs
outvcf="$DIR_MUTECT2/${TUMOR}_somatic.vcf.gz"
filteredvcf="$DIR_MUTECT2/${TUMOR}_somatic_filtered.vcf.gz"

# 1) 体細胞変異の検出 / call somatic variants
#    -normal は正常サンプルのリードグループ名 SM ( ここでは "$NORMAL" ) を指定
#    -normal must match the normal sample's read-group SM name
gatk --java-options "$JAVA_MEM" Mutect2 \
  -R "$REF" -I "$tumorbam" -I "$normalbam" -normal "$NORMAL" \
  --germline-resource "$GERMLINE_RESOURCE" \
  --pcr-indel-model AGGRESSIVE \
  -L "$INTERVALS" \
  -O "$outvcf"

# 2) 偽陽性のフィルタリング / filter likely false positives
gatk --java-options "$JAVA_MEM" FilterMutectCalls \
  -R "$REF" -V "$outvcf" -O "$filteredvcf" \
  --max-events-in-region 5

log "=== 2-5 Mutect2 完了 / done -> $DIR_MUTECT2 ==="
ls -1 "$DIR_MUTECT2"
