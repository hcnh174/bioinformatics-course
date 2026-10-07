# 講義1のワークシート
# *ワークシートは以下からダウンロードできるため、そちらを使ってコピー＆ペーストが推奨
# https://home.hiroshima-u.ac.jp/nhayes/course-site/#/course/bioinformatics/2026

# パイプライン化して実行する

# 1 bash スクリプト

# genome 環境をアクティベート
conda activate genome

# パイプラインディレクトリのコピーを作成する
cp -r ~/data/bash_pipeline ~/

# パイプラインのディレクトリへ移動（配置先に合わせて変更）
cd ~/bash_pipeline

# 全ステップを実行
bash run_all.sh

# 特定のステップのみ実行
bash run_all.sh 05_mutect2_somatic

#####################################################################
#2 Snakemake

# 作成した環境をアクティベート
conda activate snakemake

# パイプラインディレクトリのコピーを作成する
cp -r ~/data/snakemake_pipeline ~/

# パイプラインのディレクトリへ移動
cd ~/snakemake_pipeline

# ドライラン（実行せずに処理の流れ＝DAG を確認）
snakemake -n

# 4 コアで実行（conda 環境を自動構築、mamba を使用）
snakemake -c4 --use-conda --conda-frontend mamba

#####################################################################

# 3 Nextflow
# 作成した環境をアクティベ
ート
conda activate nextflow

# パイプラインディレクトリのコピーを作成する
cp -r ~/data/nextflow_pipeline ~/

# パイプラインのディレクトリへ移動
cd ~/nextflow_pipeline

# 配線の確認（ツール不要）
nextflow run main.nf -stub-run

# conda 環境を自動構築して実行
nextflow run main.nf -profile conda
