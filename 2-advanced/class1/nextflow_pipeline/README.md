# Genome Analysis パイプライン (Nextflow DSL2 / conda + modular)

Lecture 4 Worksheet「Part 2. Genome Analysis」の Nextflow 版です。Snakemake 版と同じ3点に対応しています。

1. **各プロセスに `conda` ディレクティブ** — `-profile conda` でツールごとの環境を自動構築
2. **プロセスを別ファイルに分離** — `modules/*.nf` に定義し `main.nf` から `include`
3. **入力 FASTQ を TSV から取得** — `samples.tsv` に書いたパスを参照

The Nextflow version of Part 2, matching the Snakemake build: per-process `conda` envs,
processes split into `modules/*.nf`, and inputs taken from a `samples.tsv` sheet.

## 依存関係（DAG）/ Pipeline DAG

`dag.png` を参照。`splitCsv` で TSV を読み、`join` / `branch` / `combine` で
腫瘍・正常を扱い、Mutect2（体細胞）と HaplotypeCaller（生殖細胞系列）に至ります。

## 構成 / Layout

```
nf_pipeline/
├── main.nf                # ワークフロー本体（TSV 読み込み・チャネル配線）
├── nextflow.config        # params・conda/docker プロファイル・レポート
├── samples.tsv            # 入力サンプルシート / input sample sheet
├── modules/               # プロセス定義（工程ごと）/ process definitions
│   ├── fastp.nf
│   ├── bwa.nf
│   ├── markduplicates.nf
│   ├── bqsr.nf            # BASERECALIBRATOR + APPLYBQSR
│   ├── mutect2.nf         # MUTECT2 + FILTER_MUTECT
│   └── haplotypecaller.nf # HAPLOTYPECALLER + SELECT/FILTER (SNP/INDEL) + MERGE
├── envs/                  # conda 環境定義 / per-tool conda envs
│   ├── fastp.yaml  bwa.yaml  picard.yaml  gatk.yaml  bcftools.yaml
├── dag.png                # 依存関係の図 / pipeline DAG
└── README.md
```

## サンプルシート / Sample sheet (`samples.tsv`)

タブ区切り。`type` 列で腫瘍・正常を指定します（Mutect2 / HaplotypeCaller が使用）。
各行が `meta = [id, type]` になり、チャネルを流れます。

```
sample	type	fq1	fq2
tumor	tumor	~/data/dna/tumor_1.fq.gz	~/data/dna/tumor_2.fq.gz
normal	normal	~/data/dna/normal_1.fq.gz	~/data/dna/normal_2.fq.gz
```

- サンプル追加は行を足すだけ / add a row per sample。`~` はホームに展開されます。
- 別シートは `--samples other.tsv` で切替 / switch sheets with `--samples`.

## conda 環境 / Conda environments

各プロセスに `conda "${projectDir}/envs/xxx.yaml"` を指定。`-profile conda` を付けると
Nextflow が各環境を初回に自動構築します（`nextflow.config` で `conda.useMamba = true`）。

Each process declares `conda "${projectDir}/envs/xxx.yaml"`. With `-profile conda`,
Nextflow builds each environment on first run.

## 前提 / Prerequisites

- Nextflow（Java 17+）/ Nextflow (Java 17+):

  ```bash
  curl -s https://get.nextflow.io | bash    # or: conda install -c bioconda nextflow
  ```

- conda / mamba（`-profile conda` を使う場合）
- 参照ファイル（`nextflow.config` の `~/data` 以下）。参照ゲノムには `.fai` と `.dict`、
  BWA インデックス（`*.bwt` 等）が必要 / references need `.fai`, `.dict`, and the BWA index.

## 使い方 / Usage

```bash
# 配線の確認（ツール不要、stub で空ファイルを生成）/ validate wiring with stubs
nextflow run main.nf -stub-run

# conda 環境を自動構築して実行 / run with per-process conda envs
nextflow run main.nf -profile conda

# サンプルシートや設定を上書き / override inputs & params
nextflow run main.nf -profile conda --samples my_samples.tsv --outdir ~/genome_run2 --threads 8

# 中断後の再開（キャッシュ利用）/ resume using the cache
nextflow run main.nf -profile conda -resume
```

出力は `--outdir`（既定 `~/genome`）の `01fastp` … `06haplotypecaller` に `publishDir` で
コピーされます。実行レポート（report/timeline/dag）は `<outdir>/reports/` に出ます。

## プロセス / Processes

| module | worksheet | 内容 |
|---|---|---|
| `modules/fastp.nf` (`FASTP`) | 2-1 | 品質トリミング（TSV から入力）|
| `modules/bwa.nf` (`BWA_MAP`) | 2-2 | アライメント + 整列 + 索引 |
| `modules/markduplicates.nf` (`MARKDUPLICATES`) | 2-3 | PCR 重複マーキング |
| `modules/bqsr.nf` (`BASERECALIBRATOR`, `APPLYBQSR`) | 2-4 | 塩基品質再校正 |
| `modules/mutect2.nf` (`MUTECT2`, `FILTER_MUTECT`) | 2-5 | 体細胞変異（腫瘍 vs 正常）|
| `modules/haplotypecaller.nf` (`HAPLOTYPECALLER`, `SELECT_SNP`, `SELECT_INDEL`, `FILTER_SNP`, `FILTER_INDEL`, `MERGE_GERMLINE`) | 2-6 | 生殖細胞系列変異 |

ツール・パラメータ（fastp フラグ、`@RG`、`-Xmx10g`、`--pcr-indel-model AGGRESSIVE`、
VariantFiltration のフィルタ式など）はワークシート通りです。

## Snakemake 版との対応 / Mapping to the Snakemake version

| Snakemake | Nextflow |
|---|---|
| `rules/*.smk` を `include` | `modules/*.nf` を `include` |
| ルールの `conda:` | プロセスの `conda` ディレクティブ |
| `--use-conda` | `-profile conda` |
| ワイルドカード `{sample}` | `meta` を含むチャネル（`tuple(meta, ...)`）|
| 出力ファイル名で依存解決 | チャネルの受け渡し（`join` / `combine` / `branch`）|
| `snakemake -n` | `nextflow run main.nf -stub-run` |

補足 / note: 参照・既知サイトはワークシート同様、絶対パスで直接参照します（大きく副ファイルも
多いため）。サンプル由来のファイル（FASTQ・BAM・VCF）はチャネルでステージングされます。
`--threads` は既定 4（コア数が少ない環境では下げてください）。
