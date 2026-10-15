
# Class 2 - Worksheet 2 - Nanopore sequencing
# クラス2 - ワークシート2 - ナノポアシーケンシング

################################################################
# Inspecting and exporting reads from a POD5 file fom a plasmid
# プラスミドのPOD5ファイルからリードを検査およびエクスポートする
################################################################

# Start a screen session so that you can return to the same shell later
# 後で同じシェルに戻れるよう、screenセッションを開始する
screen -S nanopore

# To reconnect later. run this command to restore your session
# 後で再接続するには、このコマンドを実行してセッションを復元してください
# screen -r nanopore

# Create variables pointing to the input and output directories
# 入力ディレクトリおよび出力ディレクトリを指す変数を作成する
indir=~/data/nanopore/pod5s/plasmid_A_pod5
outdir=~/nanopore/plasmid

# Create and move to the output folder
# 出力フォルダを作成し、そこに移動する
mkdir -p $outdir && cd $outdir

# create a variable pointing to the POD5 file containing the raw Nanopore data
# ナノポアの生データを含むPOD5ファイルを指す変数を作成する
pod5file=$indir/FAW10003_pass_barcode02_c970ecee_8ab3e29f_0.pod5

# Use the pod5 tool to inspect the POD5 file (press q to quit)
# pod5ツールを使用してPOD5ファイルを調査します（終了するにはqキーを押してください）
pod5 inspect reads $pod5file | less -S

# Print a summary of the read data
# 読み込んだデータの概要を出力する
pod5 view $pod5file | less -S

# Optional: Export the reads to FASTQ format using dorado (takes a long time and not necessary)
# 任意：doradoを使用してリードをFASTQ形式でエクスポートする（時間がかかるため、必須ではありません）
dorado basecaller hac $pod5file --emit-fastq --max-reads 10 > $outdir/plasmid_reads.fastq

# Align and reads using minimap2 and export the reads to BAM format (takes a long time)
# minimap2を用いてリードをアライメントし、その結果をBAM形式で出力する（長時間かかる）
dorado basecaller hac $pod5file > $outdir/plasmid_reads.bam

# Inspect the BAM file
# BAMファイルを調査する
samtools view -h $outdir/plasmid_reads.bam | head -n 10 | less -S

################################################################
# Variant calling part of chromosome 22 using Clair3
# Clair3を用いた22番染色体におけるバリアントコーリング
################################################################

# Create variables pointing to the input and output directories
# 入力ディレクトリおよび出力ディレクトリを指す変数を作成する
indir=~/data/nanopore/chr20
outdir=~/nanopore/chr20

# Create and move to the output folder
# 出力フォルダを作成し、そこに移動する
mkdir -p $outdir && cd $outdir

# Activate the clair3 conda environment
# clair3 conda環境を有効化する
conda activate clair3

# Define variables to simplify the command
# コマンドを簡略化するために変数を定義する
bamfile=$indir/HG003_chr20_demo.bam
reffile=$indir/GRCh38_no_alt_chr20.fa
model=ont

# inspect the header of the BAM file
# BAMファイルのヘッダーを確認する
samtools view -H $bamfile

# Inspect the entries in the BAM file. Note the messy CIGAR strings!
# BAMファイル内の各エントリーを確認する。CIGAR文字列が乱雑になっている点に注目！
samtools view $bamfile | head -n 40 | less -S

# Run Clair3
# Clair3を実行する
run_clair3.sh --bam_fn=$bamfile --ref_fn=$reffile --platform=ont \
  --model_path=${CONDA_PREFIX}/bin/models/$model --threads=4  --output=$outdir

# Inspect the called variants. Note the high frequency of indels
# 検出されたバリアントを調査する。インデルの頻度が高いことに注目すること。
zcat $outdir/merge_output.vcf.gz | less

# Switch to the nanopore conda environment
# nanoporeのconda環境に切り替える
conda activate nanopore

# Export statistics about the called variants
# 検出されたバリアントに関する統計情報をエクスポートする
bcftools stats $outdir/merge_output.vcf.gz > $outdir/vcf_stats.txt

# Inspect the variant information
# バリアント情報を確認する
grep "^SN" $outdir/vcf_stats.txt
