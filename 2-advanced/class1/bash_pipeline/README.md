# Genome Analysis パイプライン (bash) / Genome Analysis Pipeline

Lecture 4 Worksheet「Part 2. Genome Analysis」を、**bash スクリプトのみ**で構成した
実行可能なパイプラインです。fastp から体細胞・生殖細胞系列の変異検出までを行います。

A pure-**bash** pipeline reproducing *Part 2. Genome Analysis* of the Lecture 4 worksheet,
from `fastp` through somatic and germline variant calling.

## 処理の流れ / Steps

| # | スクリプト / script | 内容 / step | 主な出力 / main output |
|---|---|---|---|
| 2-1 | `01_fastp.sh` | 品質トリミング / trimming | `01fastp/{sample}_{1,2}_trimmed.fastq.gz` |
| 2-2 | `02_bwa.sh` | アライメント / alignment | `02bwa/{sample}.bam` |
| 2-3 | `03_markduplicates.sh` | PCR 重複マーキング / mark duplicates | `03markduplicate/{sample}_marked_duplicates.bam` |
| 2-4 | `04_bqsr.sh` | 塩基品質再校正 / BQSR | `04baserecalibrator/{sample}_bqsr.bam` |
| 2-5 | `05_mutect2_somatic.sh` | 体細胞変異 / somatic (Mutect2) | `05variantcall/tumor_somatic_filtered.vcf.gz` |
| 2-6 | `06_haplotypecaller_germline.sh` | 生殖細胞系列変異 / germline (HaplotypeCaller) | `06haplotypecaller/normal_germline_merged_filtered.vcf.gz` |

`config.sh` … 全パスとパラメータ / all paths & parameters
`lib_common.sh` … 共通関数（ログ・入力チェック）/ shared helpers
`run_all.sh` … 一括実行 / run everything

## 前提 / Prerequisites

conda 環境 `genome`（ワークシート Part 1 参照）/ the `genome` conda env from Part 1:

```bash
conda create -c conda-forge -c bioconda -n genome \
    fastp bwa gatk4 samtools picard bcftools -y
conda activate genome
```

入力データと参照ファイルは既定で `~/course` 以下にある想定です
（`config.sh` で変更可）/ inputs & references are expected under `~/course`
by default (change them in `config.sh`).

## 使い方 / Usage

```bash
conda activate genome

# 全ステップを順に実行 / run the full pipeline
bash run_all.sh

# 特定のステップだけ実行 / run selected steps only
bash run_all.sh 05_mutect2_somatic

# 設定を一時的に上書き / override settings on the fly
THREADS=8 WORK_DIR=~/genome_run2 bash run_all.sh
```

各スクリプトは単体でも実行できます / each script can be run on its own:

```bash
bash 01_fastp.sh
```

## 備考 / Notes

- 出力は既定で `~/genome/01fastp` … `~/genome/06haplotypecaller` に作られます。
  Outputs default to `~/genome/01fastp` … `~/genome/06haplotypecaller`.
- 各ステップ冒頭で必要なコマンド・入力ファイルの存在を確認し、
  問題があれば分かりやすいメッセージで停止します（`set -euo pipefail`）。
  Every step checks required commands and inputs and stops with a clear message on error.
- パラメータ（fastp のフラグ、`-Xmx`、スレッド数、フィルタ式など）は
  ワークシートの値をそのまま採用しています。
  Parameters follow the worksheet values as-is.
