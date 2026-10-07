#!/usr/bin/env bash
# ============================================================
# config.sh — パイプライン共通設定 / Shared pipeline configuration
# ------------------------------------------------------------
# ここの変数を書き換えれば全ステップに反映されます。
# Edit these variables once; every step picks them up.
# 環境変数で上書きも可能 / override from the shell, e.g.:
#     THREADS=8 WORK_DIR=~/genome2 bash run_all.sh
# ============================================================

# --- ベースディレクトリ / Base directories ---
: "${COURSE_DIR:=$HOME/data}"     # 講義データ（入力・読み取り専用）/ course data (inputs)
: "${WORK_DIR:=$HOME/pipeline_res/bash}"       # 解析結果の出力先 / analysis output root

# --- 生データ / Raw paired-end FASTQ ( *_1.fq.gz / *_2.fq.gz ) ---
: "${RAW_DIR:=$COURSE_DIR/dna}"

# --- 参照ゲノムとインデックス / Reference genome & index ---
: "${REF:=$COURSE_DIR/ref/hg19/hg19.fa}"
: "${BWA_INDEX:=$COURSE_DIR/ref/hg19/bwa_index/hg19_bwa}"

# --- BQSR 用の既知サイト / Known sites for BaseRecalibrator ---
: "${KNOWN1:=$COURSE_DIR/vcf/1000G_omni2.5.hg19.sites.vcf}"
: "${KNOWN2:=$COURSE_DIR/vcf/1000G_phase1.indels.hg19.sites.vcf}"
: "${KNOWN3:=$COURSE_DIR/vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf}"
: "${KNOWN4:=$COURSE_DIR/vcf/dbsnp_138.hg19.vcf}"

# --- Mutect2 用リソース / Mutect2 resources ---
: "${GERMLINE_RESOURCE:=$COURSE_DIR/vcf/af-only-gnomad.raw.sites.chr.vcf}"
: "${INTERVALS:=$COURSE_DIR/dna/braf_hg19_exon15.bed}"

# --- サンプル / Samples ( 腫瘍・正常 / tumor & normal ) ---
: "${TUMOR:=tumor}"
: "${NORMAL:=normal}"
SAMPLES=("$TUMOR" "$NORMAL")

# --- 実行パラメータ / Runtime parameters ---
: "${THREADS:=4}"          # スレッド数 / threads
: "${JAVA_MEM:=-Xmx10g}"   # GATK / Picard の Java ヒープ / Java heap size

# --- fastp パラメータ / fastp parameters ---
: "${FASTP_N:=10}"   # -n : N 塩基数の上限 / N base number limit
: "${FASTP_L:=15}"   # -l : 最小リード長 / min read length
: "${FASTP_Q:=20}"   # -q : クオリティ閾値 / qualified quality
: "${FASTP_T:=1}"    # -t : read1 の 3' トリム / trim tail of read1
: "${FASTP_TT:=1}"   # -T : read2 の 3' トリム / trim tail of read2
: "${FASTP_W:=4}"    # -w : fastp のスレッド数 / fastp worker threads

# --- 各ステップの出力ディレクトリ / Per-step output directories ---
: "${DATA_DIR:=$WORK_DIR/data}"
: "${DIR_FASTP:=$WORK_DIR/01fastp}"
: "${DIR_BWA:=$WORK_DIR/02bwa}"
: "${DIR_MARKDUP:=$WORK_DIR/03markduplicate}"
: "${DIR_BQSR:=$WORK_DIR/04baserecalibrator}"
: "${DIR_MUTECT2:=$WORK_DIR/05variantcall}"
: "${DIR_HC:=$WORK_DIR/06haplotypecaller}"
