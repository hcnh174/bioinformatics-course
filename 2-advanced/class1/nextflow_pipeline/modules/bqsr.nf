// 2-4 BQSR — 再校正テーブル作成 -> 適用 / recalibrate ( per sample )
process BASERECALIBRATOR {
    tag "${meta.id}"
    conda "${projectDir}/envs/gatk.yaml"
    publishDir "${params.outdir}/04baserecalibrator", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_recalibration.table"), emit: table

    script:
    """
    gatk --java-options "${params.java_mem}" BaseRecalibrator \\
      -I ${bam} -R ${params.ref} ${params.known_sites.collect{ '--known-sites ' + it }.join(' ')} \\
      -O ${meta.id}_recalibration.table
    """

    stub:
    """
    touch ${meta.id}_recalibration.table
    """
}

process APPLYBQSR {
    tag "${meta.id}"
    conda "${projectDir}/envs/gatk.yaml"
    publishDir "${params.outdir}/04baserecalibrator", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(table)

    output:
    tuple val(meta), path("${meta.id}_bqsr.bam"), path("${meta.id}_bqsr.bai"), emit: bam

    script:
    """
    gatk --java-options "${params.java_mem}" ApplyBQSR \\
      -R ${params.ref} -I ${bam} --bqsr-recal-file ${table} \\
      -O ${meta.id}_bqsr.bam
    """

    stub:
    """
    touch ${meta.id}_bqsr.bam ${meta.id}_bqsr.bai
    """
}
