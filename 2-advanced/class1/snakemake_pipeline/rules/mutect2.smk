# 2-5 Mutect2 — 体細胞変異（腫瘍 vs 正常）/ somatic variant calling
rule mutect2:
    input:
        tumor     = f"{D_BQSR}/{TUMOR}_bqsr.bam",
        normal    = f"{D_BQSR}/{NORMAL}_bqsr.bam",
        ref       = REF,
        germline  = GERMLINE_RESOURCE,
        intervals = INTERVALS,
    output:
        vcf = f"{D_MUTECT2}/{TUMOR}_somatic.vcf.gz",
    log:
        f"{D_MUTECT2}/logs/mutect2.log",
    params:
        mem    = JAVA_MEM,
        normal = NORMAL,
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.mem}" Mutect2 \
          -R {input.ref} -I {input.tumor} -I {input.normal} -normal {params.normal} \
          --germline-resource {input.germline} --pcr-indel-model AGGRESSIVE \
          -L {input.intervals} -O {output.vcf} 2> {log}
        """

rule filter_mutect:
    input:
        vcf = f"{D_MUTECT2}/{TUMOR}_somatic.vcf.gz",
        ref = REF,
    output:
        vcf = f"{D_MUTECT2}/{TUMOR}_somatic_filtered.vcf.gz",
    log:
        f"{D_MUTECT2}/logs/filter_mutect.log",
    params:
        mem = JAVA_MEM,
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.mem}" FilterMutectCalls \
          -R {input.ref} -V {input.vcf} -O {output.vcf} --max-events-in-region 5 2> {log}
        """
