# Genome Analysis パイプライン (Snakemake / conda + modular)

Lecture 4 Worksheet「Part 2. Genome Analysis」の Snakemake 版です。今回は次の3点に対応しています。

1. **各ルールに `conda:` ディレクティブ** — `--use-conda` でツールごとの環境を自動構築
2. **ルールを別ファイルに分離** — `rules/*.smk` に定義し `Snakefile` から `include`
3. **入力 FASTQ を TSV から取得** — `samples.tsv` に書いたパスを参照

Same Part 2 pipeline, now with per-rule `conda:` envs, rules split into `rules/*.smk`,
and input FASTQ paths taken from a `samples.tsv` sheet.

## ディレクトリ構成 / Layout

```
smk_conda/
├── Snakefile              # 設定読み込み・サンプルシート解析・rule all・include
├── config.yaml            # パス・パラメータ / paths & parameters
├── samples.tsv            # 入力サンプルシート / input sample sheet
├── rules/                 # ルール定義（工程ごと）/ rule definitions per step
│   ├── fastp.smk
│   ├── bwa.smk
│   ├── markduplicates.smk
│   ├── bqsr.smk
│   ├── mutect2.smk
│   └── haplotypecaller.smk
├── envs/                  # conda 環境定義 / per-tool conda envs
│   ├── fastp.yaml
│   ├── bwa.yaml           # bwa + samtools
│   ├── picard.yaml
│   ├── gatk.yaml          # gatk4
│   └── bcftools.yaml
├── rulegraph.png          # 依存関係の図 / rule graph
└── README.md
```

## サンプルシート / Sample sheet (`samples.tsv`)

タブ区切りで、入力 FASTQ のパスを記載します。`type` 列で腫瘍・正常を指定します
（Mutect2 / HaplotypeCaller はこの列を使います）。

Tab-separated. `type` marks tumor vs normal (used by Mutect2 / HaplotypeCaller).

```
sample	type	fq1	fq2
tumor	tumor	~/data/dna/tumor_1.fq.gz	~/data/dna/tumor_2.fq.gz
normal	normal	~/data/dna/normal_1.fq.gz	~/data/dna/normal_2.fq.gz
```

- サンプルを増やすには行を追加するだけ / add a row per sample.
- `~` はホームに展開されます / `~` is expanded to the home directory.
- 別のシートを使う場合は `config.yaml` の `samples:` を変更、または
  `snakemake --config samples=other.tsv` で上書き。

## conda 環境 / Conda environments

各ルールに `conda: "../envs/xxx.yaml"` を指定しています。`--use-conda` を付けると、
Snakemake が各ツールの環境を初回に自動構築して実行します（`genome` 環境を手動で
用意する必要はありません）。

Each rule declares `conda: "../envs/xxx.yaml"`. With `--use-conda`, Snakemake builds each
tool's environment automatically on first run.

## 前提 / Prerequisites

- Snakemake と conda（mamba 推奨）/ Snakemake and conda (mamba recommended):

  ```bash
  conda install -c conda-forge -c bioconda snakemake -y
  ```

- 参照ファイル（`config.yaml` の `~/data` 以下）。参照ゲノムには `.fai` と `.dict`、
  BWA インデックス（`*.bwt` 等）が必要 / references need `.fai`, `.dict`, and the BWA index.

## 使い方 / Usage

```bash
# ドライラン（DAG の確認、実行はしない）/ dry-run
snakemake -n

# conda 環境を自動構築して 4 コアで実行 / run with conda, 4 cores
snakemake -c4 --use-conda -p

# さらに mamba で環境構築を高速化 / speed up env creation with mamba
snakemake -c4 --use-conda --conda-frontend mamba

# 特定の出力だけ / build one target
snakemake -c4 --use-conda ~/genome/05variantcall/tumor_somatic_filtered.vcf.gz

# 依存グラフを描く / draw the rule graph
snakemake --rulegraph | dot -Tpng > rulegraph.png
```

設定の一時上書き / override on the fly:

```bash
snakemake -c8 --use-conda --config threads=8 work_dir=~/genome_run2 samples=other.tsv
```

## 工程 / Steps (rules)

| rule ファイル | worksheet | 内容 |
|---|---|---|
| `rules/fastp.smk` | 2-1 | 品質トリミング（TSV から入力）/ trimming (inputs from TSV) |
| `rules/bwa.smk` | 2-2 | アライメント + 整列 + 索引 |
| `rules/markduplicates.smk` | 2-3 | PCR 重複マーキング |
| `rules/bqsr.smk` | 2-4 | `baserecalibrator` → `applybqsr` |
| `rules/mutect2.smk` | 2-5 | `mutect2` → `filter_mutect`（体細胞）|
| `rules/haplotypecaller.smk` | 2-6 | `haplotypecaller` → `select_snp`/`select_indel` → `filter_snp`/`filter_indel` → `merge_germline`（生殖細胞系列）|

ツール・パラメータ（fastp フラグ、`@RG`、`-Xmx10g`、`--pcr-indel-model AGGRESSIVE`、
VariantFiltration のフィルタ式など）はワークシート通りです。各ルールはログを
`<出力ディレクトリ>/logs/` に書き出します。
