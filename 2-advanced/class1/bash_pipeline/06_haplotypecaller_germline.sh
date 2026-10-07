#!/usr/bin/env bash
# ============================================================
# 2-6 HaplotypeCaller — 生殖細胞系列変異の検出（正常サンプル）
#                     Germline variant calling (normal, GATK)
#   HaplotypeCaller -> SNP/INDEL 分割 -> それぞれフィルタ -> 結合
#   call -> split SNP/INDEL -> filter each -> merge
# 入力 / in : $DIR_BQSR/normal_bqsr.bam
# 出力 / out: $DIR_HC/normal_germline_merged_filtered.vcf.gz
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
source "$HERE/lib_common.sh"

log "=== 2-6 HaplotypeCaller (germline) 開始 / start ==="
require_cmd gatk
require_cmd bcftools
require_file "$REF"
require_file "$INTERVALS"
mkdir -p "$DIR_HC"

# 入力 / input
normalbam="$DIR_BQSR/${NORMAL}_bqsr.bam"
require_file "$normalbam"

# 出力 / outputs
outvcf="$DIR_HC/${NORMAL}_germline.vcf.gz"
snpvcf="$DIR_HC/${NORMAL}_germline_snp.vcf.gz"
indelvcf="$DIR_HC/${NORMAL}_germline_indel.vcf.gz"
filteredsnp="$DIR_HC/${NORMAL}_germline_snp_filtered.vcf.gz"
filteredindel="$DIR_HC/${NORMAL}_germline_indel_filtered.vcf.gz"
mergedvcf="$DIR_HC/${NORMAL}_germline_merged_filtered.vcf.gz"

# 1) 生殖細胞系列変異の検出 / call germline variants
gatk --java-options "$JAVA_MEM" HaplotypeCaller \
  -R "$REF" -I "$normalbam" -O "$outvcf" -L "$INTERVALS"

# 2) SNP と INDEL に分割 / split into SNPs and INDELs
gatk --java-options "$JAVA_MEM" SelectVariants \
  -R "$REF" -V "$outvcf" --select-type-to-include SNP -O "$snpvcf"
gatk --java-options "$JAVA_MEM" SelectVariants \
  -R "$REF" -V "$outvcf" --select-type-to-include INDEL -O "$indelvcf"

# 3) それぞれ個別にフィルタリング / filter SNPs and INDELs separately
gatk --java-options "$JAVA_MEM" VariantFiltration \
  -R "$REF" -V "$snpvcf" \
  --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \
  --filter-name 'snp_filter' -O "$filteredsnp"

gatk --java-options "$JAVA_MEM" VariantFiltration \
  -R "$REF" -V "$indelvcf" \
  --filter-expression 'QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0 || QUAL < 30.0' \
  --filter-name 'indel_filter' -O "$filteredindel"

# 4) SNP と INDEL を結合 / merge SNPs and INDELs back together
bcftools concat "$filteredsnp" "$filteredindel" -a -Oz -o "$mergedvcf"

log "=== 2-6 HaplotypeCaller 完了 / done -> $DIR_HC ==="
ls -1 "$DIR_HC"
