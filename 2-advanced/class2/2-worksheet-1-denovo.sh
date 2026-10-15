
# Class 2 - Worksheet 1 - De novo assembly
# クラス2 - ワークシート1 - デノボ・アセンブリ

#########################################################################################

# Start a screen session so that you can return to the same shell later
# 後で同じシェルに戻れるよう、screenセッションを開始する
screen -S denovo

# To reconnect later. run this command to restore your session
# 後で再接続するには、このコマンドを実行してセッションを復元してください
# screen -r denovo

# Activate the denovo conda environment
# denovo conda環境を有効化する
conda activate denovo

# Create variables pointing to the input and output directories
# 入力ディレクトリおよび出力ディレクトリを指す変数を作成する
indir=~/data/denovo
outdir=~/denovo

# Create the output and qc folders
# 出力用フォルダと品質管理（QC）用フォルダを作成する
mkdir -p $outdir $outdir/qc

# Move to the denovo output directory
# de novo出力ディレクトリに移動する
cd $outdir

# Generate quality control reports for each fastq file using fastqc
# FastQCを使用して、各FASTQファイルの品質管理レポートを生成する
fastqc $indir/*.fastq -o $outdir/qc/

# Review the output files in your Jupyter notebook file browser: http://133.41.103.64:8000/
# Jupyter Notebookのファイルブラウザで出力ファイルを確認してください: http://133.41.103.64:8000/
# $outdir/qc/SRR40288597_1_fastqc.html $outdir/qc/SRR40288597_2_fastqc.html

# Trim and filter the reads using fastp
# fastpを使用してリードをトリミングおよびフィルタリングします
fastp -i $indir/SRR40288597_1.fastq -I $indir/SRR40288597_2.fastq \
    -o $outdir/phage_1.trim.fq.gz -O $outdir/phage_2.trim.fq.gz \
    --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 50 \
    --json $outdir/qc/fastp.json --html $outdir/qc/fastp.html

# review the fastp output file: $outdir/qc/fastp.html
# fastpの出力ファイルを確認する: $outdir/qc/fastp.html

# Assemble the reads with Spades. Exclude--only-assembler and -k 21,55 should during an actual run
# SPAdesを用いてリードをアセンブルする。実際の実行時には、`--only-assembler` および `-k 21,55` を指定すること。
spades.py -1 $outdir/phage_1.trim.fq.gz -2 $outdir/phage_2.trim.fq.gz \
	-o $outdir/spades --isolate --threads 4 --memory 8 \
    --only-assembler -k 21,55  # Added only to speed up processing.  remove during actual run
    
# Inspect the final scaffolds file. How many separate contigs still remain?
# 最終的なスキャフォールドファイルを調査する。 個別のコンティグはいくつ残っていますか？
cat $outdir/spades/scaffolds.fasta

# Evaluate the quality of the genome assembly using Quast
# Quastを用いてゲノムアセンブリの品質を評価する
quast.py $outdir/spades/contigs.fasta -r $indir/phage_t7.fasta -o $outdir/quast --min-contig 200

# Inspect the output of the report.txt file
cat $outdir/quast/report.txt

# Create and view an image of the assembly graph using Bandage
# Bandageを使用してアセンブリグラフの画像を作成・表示する
Bandage image $outdir/spades/assembly_graph.fastg $outdir/spades/bandage_assembly_graph.png

# Optional: Look for overlapping ends in the contig. Warning: takes a very long time to run!
# 任意：コンティグ内で末端の重複を探す。注意：実行に非常に長い時間がかかります！
# circlator all $outdir/spades/contigs.fasta $outdir/phage_1.trim.fq.gz $outdir/circlator/

# Validate against a known reference sequence using minimap2 to align reads
# minimap2を用いてリードをアライメントし、既知の参照配列に対して検証を行う
minimap2 -ax sr $indir/phage_t7.fasta $outdir/phage_1.trim.fq.gz $outdir/phage_2.trim.fq.gz | \
	samtools sort -o $outdir/mapped.bam
samtools index $outdir/mapped.bam
samtools flagstat $outdir/mapped.bam

# If X11 forwarding is eanbled, inspect the aligned reads in IGV
# X11フォワーディングが有効な場合は、IGVでアライメント済みリードを確認してください。
igv.sh $outdir/mapped.bam --genome $indir/phage_t7.fasta
