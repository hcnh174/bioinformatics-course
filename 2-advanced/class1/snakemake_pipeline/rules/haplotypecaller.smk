# 2-6 HaplotypeCaller — 生殖細胞系列変異（正常）/ germline variant calling
#   call -> SNP/INDEL 分割 -> それぞれフィルタ -> 結合
rule haplotypecaller:
    input:
        bam       = f"{D_BQSR}/{NORMAL}_bqsr.bam",
        ref       = REF,
        intervals = INTERVALS,
    output:
        vcf = f"{D_HC}/{NORMAL}_germline.vcf.gz",
    log:
        f"{D_HC}/logs/haplotypecaller.log",
    params:
        mem = JAVA_MEM,
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.mem}" HaplotypeCaller \
          -R {input.ref} -I {input.bam} -O {output.vcf} -L {input.intervals} 2> {log}
        """

rule select_snp:
    input:
        vcf = f"{D_HC}/{NORMAL}_germline.vcf.gz",
        ref = REF,
    output:
        vcf = f"{D_HC}/{NORMAL}_germline_snp.vcf.gz",
    log:
        f"{D_HC}/logs/select_snp.log",
    params:
        mem = JAVA_MEM,
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.mem}" SelectVariants \
          -R {input.ref} -V {input.vcf} --select-type-to-include SNP -O {output.vcf} 2> {log}
        """

rule select_indel:
    input:
        vcf = f"{D_HC}/{NORMAL}_germline.vcf.gz",
        ref = REF,
    output:
        vcf = f"{D_HC}/{NORMAL}_germline_indel.vcf.gz",
    log:
        f"{D_HC}/logs/select_indel.log",
    params:
        mem = JAVA_MEM,
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.mem}" SelectVariants \
          -R {input.ref} -V {input.vcf} --select-type-to-include INDEL -O {output.vcf} 2> {log}
        """

rule filter_snp:
    input:
        vcf = f"{D_HC}/{NORMAL}_germline_snp.vcf.gz",
        ref = REF,
    output:
        vcf = f"{D_HC}/{NORMAL}_germline_snp_filtered.vcf.gz",
    log:
        f"{D_HC}/logs/filter_snp.log",
    params:
        mem = JAVA_MEM,
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.mem}" VariantFiltration \
          -R {input.ref} -V {input.vcf} \
          --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \
          --filter-name 'snp_filter' -O {output.vcf} 2> {log}
        """

rule filter_indel:
    input:
        vcf = f"{D_HC}/{NORMAL}_germline_indel.vcf.gz",
        ref = REF,
    output:
        vcf = f"{D_HC}/{NORMAL}_germline_indel_filtered.vcf.gz",
    log:
        f"{D_HC}/logs/filter_indel.log",
    params:
        mem = JAVA_MEM,
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.mem}" VariantFiltration \
          -R {input.ref} -V {input.vcf} \
          --filter-expression 'QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0 || QUAL < 30.0' \
          --filter-name 'indel_filter' -O {output.vcf} 2> {log}
        """

rule merge_germline:
    input:
        snp   = f"{D_HC}/{NORMAL}_germline_snp_filtered.vcf.gz",
        indel = f"{D_HC}/{NORMAL}_germline_indel_filtered.vcf.gz",
    output:
        vcf = f"{D_HC}/{NORMAL}_germline_merged_filtered.vcf.gz",
    log:
        f"{D_HC}/logs/merge_germline.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
        bcftools concat {input.snp} {input.indel} -a -Oz -o {output.vcf} 2> {log}
        """
