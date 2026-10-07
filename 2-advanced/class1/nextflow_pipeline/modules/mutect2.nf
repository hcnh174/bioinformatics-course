// 2-5 Mutect2 — 体細胞変異（腫瘍 vs 正常）/ somatic variant calling
process MUTECT2 {
    tag "${tmeta.id}_vs_${nmeta.id}"
    conda "${projectDir}/envs/gatk.yaml"
    publishDir "${params.outdir}/05variantcall", mode: 'copy'

    input:
    tuple val(tmeta), path(tbam), path(tbai), val(nmeta), path(nbam), path(nbai)

    output:
    // .stats は FilterMutectCalls が必要とするため一緒に受け渡す / carry .stats for FilterMutectCalls
    tuple val(tmeta), path("${tmeta.id}_somatic.vcf.gz"), path("${tmeta.id}_somatic.vcf.gz.tbi"), path("${tmeta.id}_somatic.vcf.gz.stats"), emit: vcf

    script:
    """
    gatk --java-options "${params.java_mem}" Mutect2 \\
      -R ${params.ref} -I ${tbam} -I ${nbam} -normal ${nmeta.id} \\
      --germline-resource ${params.germline_resource} --pcr-indel-model AGGRESSIVE \\
      -L ${params.intervals} -O ${tmeta.id}_somatic.vcf.gz
    """

    stub:
    """
    touch ${tmeta.id}_somatic.vcf.gz ${tmeta.id}_somatic.vcf.gz.tbi ${tmeta.id}_somatic.vcf.gz.stats
    """
}

process FILTER_MUTECT {
    tag "${meta.id}"
    conda "${projectDir}/envs/gatk.yaml"
    publishDir "${params.outdir}/05variantcall", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi), path(stats)

    output:
    tuple val(meta), path("${meta.id}_somatic_filtered.vcf.gz"), path("${meta.id}_somatic_filtered.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk --java-options "${params.java_mem}" FilterMutectCalls \\
      -R ${params.ref} -V ${vcf} -O ${meta.id}_somatic_filtered.vcf.gz \\
      --max-events-in-region 5
    """

    stub:
    """
    touch ${meta.id}_somatic_filtered.vcf.gz ${meta.id}_somatic_filtered.vcf.gz.tbi
    """
}
