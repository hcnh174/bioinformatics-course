// 2-3 MarkDuplicates — PCR 重複マーキング / mark duplicates ( per sample )
process MARKDUPLICATES {
    tag "${meta.id}"
    conda "${projectDir}/envs/picard.yaml"
    publishDir "${params.outdir}/03markduplicate", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_marked_duplicates.bam"), path("${meta.id}_marked_duplicates.bai"), emit: bam
    path "${meta.id}_marked_dup_metrics.txt"

    script:
    """
    picard ${params.java_mem} MarkDuplicates \\
      -I ${bam} -O ${meta.id}_marked_duplicates.bam -M ${meta.id}_marked_dup_metrics.txt
    picard ${params.java_mem} BuildBamIndex \\
      -INPUT ${meta.id}_marked_duplicates.bam -OUTPUT ${meta.id}_marked_duplicates.bai
    """

    stub:
    """
    touch ${meta.id}_marked_duplicates.bam ${meta.id}_marked_duplicates.bai ${meta.id}_marked_dup_metrics.txt
    """
}
