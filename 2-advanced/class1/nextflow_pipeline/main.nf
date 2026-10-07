#!/usr/bin/env nextflow
// ============================================================
// main.nf — Genome Analysis パイプライン (Nextflow DSL2)
//   Lecture 4 Worksheet「Part 2. Genome Analysis」を再現
//   - プロセスは modules/*.nf に分離 / processes live in modules/*.nf
//   - 各プロセスに conda ディレクティブ / per-process conda envs (-profile conda)
//   - 入力 FASTQ は samples.tsv から取得 / inputs come from a TSV sample sheet
//
//   使い方 / Usage:
//     nextflow run main.nf -stub-run          # 配線確認（ツール不要）/ wiring check
//     nextflow run main.nf -profile conda     # conda 環境を自動構築して実行 / run with conda
//
//   注記 / note: Nextflow の新しい厳格構文（25.04+）ではスクリプト先頭に
//   def などの文を置けません。~ の展開は file() が自動で行うため、そのまま file() に渡します。
//   In the strict language syntax, top-level statements (def ...) are not allowed;
//   file() expands a leading "~" on its own, so we pass paths to file() directly.
// ============================================================

// ---- モジュール読み込み / include processes ----
include { FASTP }                              from './modules/fastp.nf'
include { BWA_MAP }                            from './modules/bwa.nf'
include { MARKDUPLICATES }                     from './modules/markduplicates.nf'
include { BASERECALIBRATOR; APPLYBQSR }        from './modules/bqsr.nf'
include { MUTECT2; FILTER_MUTECT }             from './modules/mutect2.nf'
include { HAPLOTYPECALLER; SELECT_SNP; SELECT_INDEL;
          FILTER_SNP; FILTER_INDEL; MERGE_GERMLINE } from './modules/haplotypecaller.nf'

workflow {

    // ---- サンプルシート (TSV) を読み込み / read the sample sheet ----
    // columns: sample  type(tumor|normal)  fq1  fq2

    def home = System.getenv('HOME')

    ch_reads = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row ->

            def fq1_path = row.fq1.replaceFirst('^~', home)
            def fq2_path = row.fq2.replaceFirst('^~', home)

            tuple(
                [ id: row.sample, type: row.type.toLowerCase() ],
                file(fq1_path),
                file(fq2_path)
            )
        }
    // ---- 2-1 -> 2-4 : サンプルごとに前処理 / per-sample preprocessing ----
    FASTP( ch_reads )
    BWA_MAP( FASTP.out.trimmed )
    MARKDUPLICATES( BWA_MAP.out.bam )
    BASERECALIBRATOR( MARKDUPLICATES.out.bam )

    // marked BAM と recal テーブルを meta で結合 / join BAM + table by meta
    APPLYBQSR( MARKDUPLICATES.out.bam.join( BASERECALIBRATOR.out.table ) )

    // 腫瘍・正常に振り分け / split into tumor and normal
    APPLYBQSR.out.bam
        .branch { meta, bam, bai ->
            tumor : meta.type == 'tumor'
            normal: meta.type == 'normal'
        }
        .set { bqsr }

    // ---- 2-5 : 体細胞変異（腫瘍 vs 正常）/ somatic ----
    ch_pair = bqsr.tumor.combine( bqsr.normal )   // (tmeta,tbam,tbai, nmeta,nbam,nbai)
    MUTECT2( ch_pair )
    FILTER_MUTECT( MUTECT2.out.vcf )

    // ---- 2-6 : 生殖細胞系列変異（正常のみ）/ germline (normal only) ----
    HAPLOTYPECALLER( bqsr.normal )
    SELECT_SNP(   HAPLOTYPECALLER.out.vcf )
    SELECT_INDEL( HAPLOTYPECALLER.out.vcf )
    FILTER_SNP(   SELECT_SNP.out.vcf )
    FILTER_INDEL( SELECT_INDEL.out.vcf )

    // フィルタ済み SNP と INDEL を meta で結合して統合 / join filtered SNP + INDEL, then merge
    MERGE_GERMLINE( FILTER_SNP.out.vcf.join( FILTER_INDEL.out.vcf ) )
}
