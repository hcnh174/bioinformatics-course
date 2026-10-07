// 2-2 BWA — アライメント + 整列 + 索引 / align, sort, index ( per sample )
process BWA_MAP {
    tag "${meta.id}"
    conda "${projectDir}/envs/bwa.yaml"
    publishDir "${params.outdir}/02bwa", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} \\
      -R "@RG\\tID:${meta.id}\\tPL:ILLUMINA\\tSM:${meta.id}" \\
      ${params.bwa_index} ${r1} ${r2} \\
      | samtools sort -@ ${task.cpus} -o ${meta.id}.bam -
    samtools index ${meta.id}.bam
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bam.bai
    """
}
