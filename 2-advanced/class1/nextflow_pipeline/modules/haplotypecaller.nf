// 2-6 HaplotypeCaller — 生殖細胞系列変異（正常）/ germline variant calling
//   call -> SNP/INDEL 分割 -> それぞれフィルタ -> 結合
process HAPLOTYPECALLER {
    tag "${meta.id}"
    conda "${projectDir}/envs/gatk.yaml"
    publishDir "${params.outdir}/06haplotypecaller", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_germline.vcf.gz"), path("${meta.id}_germline.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk --java-options "${params.java_mem}" HaplotypeCaller \\
      -R ${params.ref} -I ${bam} -O ${meta.id}_germline.vcf.gz -L ${params.intervals}
    """

    stub:
    """
    touch ${meta.id}_germline.vcf.gz ${meta.id}_germline.vcf.gz.tbi
    """
}

process SELECT_SNP {
    tag "${meta.id}"
    conda "${projectDir}/envs/gatk.yaml"
    publishDir "${params.outdir}/06haplotypecaller", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}_germline_snp.vcf.gz"), path("${meta.id}_germline_snp.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk --java-options "${params.java_mem}" SelectVariants \\
      -R ${params.ref} -V ${vcf} --select-type-to-include SNP \\
      -O ${meta.id}_germline_snp.vcf.gz
    """

    stub:
    """
    touch ${meta.id}_germline_snp.vcf.gz ${meta.id}_germline_snp.vcf.gz.tbi
    """
}

process SELECT_INDEL {
    tag "${meta.id}"
    conda "${projectDir}/envs/gatk.yaml"
    publishDir "${params.outdir}/06haplotypecaller", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}_germline_indel.vcf.gz"), path("${meta.id}_germline_indel.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk --java-options "${params.java_mem}" SelectVariants \\
      -R ${params.ref} -V ${vcf} --select-type-to-include INDEL \\
      -O ${meta.id}_germline_indel.vcf.gz
    """

    stub:
    """
    touch ${meta.id}_germline_indel.vcf.gz ${meta.id}_germline_indel.vcf.gz.tbi
    """
}

process FILTER_SNP {
    tag "${meta.id}"
    conda "${projectDir}/envs/gatk.yaml"
    publishDir "${params.outdir}/06haplotypecaller", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}_germline_snp_filtered.vcf.gz"), path("${meta.id}_germline_snp_filtered.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk --java-options "${params.java_mem}" VariantFiltration \\
      -R ${params.ref} -V ${vcf} \\
      --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \\
      --filter-name 'snp_filter' -O ${meta.id}_germline_snp_filtered.vcf.gz
    """

    stub:
    """
    touch ${meta.id}_germline_snp_filtered.vcf.gz ${meta.id}_germline_snp_filtered.vcf.gz.tbi
    """
}

process FILTER_INDEL {
    tag "${meta.id}"
    conda "${projectDir}/envs/gatk.yaml"
    publishDir "${params.outdir}/06haplotypecaller", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}_germline_indel_filtered.vcf.gz"), path("${meta.id}_germline_indel_filtered.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk --java-options "${params.java_mem}" VariantFiltration \\
      -R ${params.ref} -V ${vcf} \\
      --filter-expression 'QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0 || QUAL < 30.0' \\
      --filter-name 'indel_filter' -O ${meta.id}_germline_indel_filtered.vcf.gz
    """

    stub:
    """
    touch ${meta.id}_germline_indel_filtered.vcf.gz ${meta.id}_germline_indel_filtered.vcf.gz.tbi
    """
}

process MERGE_GERMLINE {
    tag "${meta.id}"
    conda "${projectDir}/envs/bcftools.yaml"
    publishDir "${params.outdir}/06haplotypecaller", mode: 'copy'

    input:
    tuple val(meta), path(snp), path(snp_tbi), path(indel), path(indel_tbi)

    output:
    tuple val(meta), path("${meta.id}_germline_merged_filtered.vcf.gz"), emit: vcf

    script:
    """
    bcftools concat ${snp} ${indel} -a -Oz -o ${meta.id}_germline_merged_filtered.vcf.gz
    """

    stub:
    """
    touch ${meta.id}_germline_merged_filtered.vcf.gz
    """
}
