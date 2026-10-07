#!/usr/bin/env bash
# ============================================================
# 2-1 fastp — アダプター除去・品質トリミング
#          Adapter removal & quality trimming
# 入力 / in : $RAW_DIR/{tumor,normal}_{1,2}.fq.gz
# 出力 / out: $DIR_FASTP/{sample}_{1,2}_trimmed.fastq.gz  ほか
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/config.sh"
source "$HERE/lib_common.sh"

log "=== 2-1 fastp 開始 / start ==="
require_cmd fastp

# データを作業ディレクトリへ配置 / stage raw data into the work dir
mkdir -p "$DATA_DIR" "$DIR_FASTP"
cp -n "$RAW_DIR"/*.fq.gz "$DATA_DIR"/ 2>/dev/null || true

# 腫瘍・正常を同じ処理でループ / same processing for tumor and normal
for sample in "${SAMPLES[@]}"; do
  log "fastp: $sample"

  # 入力 / inputs
  fq1="$DATA_DIR/${sample}_1.fq.gz"
  fq2="$DATA_DIR/${sample}_2.fq.gz"
  require_file "$fq1"
  require_file "$fq2"

  # 出力 / outputs
  trim1="$DIR_FASTP/${sample}_1_trimmed.fastq.gz"
  trim2="$DIR_FASTP/${sample}_2_trimmed.fastq.gz"
  unp1="$DIR_FASTP/${sample}_1_unpaired.fastq.gz"
  unp2="$DIR_FASTP/${sample}_2_unpaired.fastq.gz"
  html="$DIR_FASTP/${sample}_fastp.html"
  json="$DIR_FASTP/${sample}_fastp.json"

  # fastp 実行 / run fastp  ( -3: 3' 側を品質でトリム / trim 3' by quality )
  fastp --in1 "$fq1" --in2 "$fq2" -3 \
    --out1 "$trim1" --out2 "$trim2" \
    --unpaired1 "$unp1" --unpaired2 "$unp2" \
    -h "$html" -j "$json" \
    -n "$FASTP_N" -l "$FASTP_L" -q "$FASTP_Q" -t "$FASTP_T" -T "$FASTP_TT" -w "$FASTP_W"
done

log "=== 2-1 fastp 完了 / done -> $DIR_FASTP ==="
ls -1 "$DIR_FASTP"
