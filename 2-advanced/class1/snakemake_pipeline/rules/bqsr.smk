# 2-4 BQSR — 再校正テーブル作成 -> 適用 / recalibrate ( per sample )
rule baserecalibrator:
    input:
        bam   = f"{D_MARKDUP}/{{sample}}_marked_duplicates.bam",
        bai   = f"{D_MARKDUP}/{{sample}}_marked_duplicates.bai",
        ref   = REF,
        known = KNOWN,
    output:
        table = f"{D_BQSR}/{{sample}}_recalibration.table",
    log:
        f"{D_BQSR}/logs/{{sample}}_baserecalibrator.log",
    params:
        mem       = JAVA_MEM,
        known_arg = KNOWN_ARG,
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.mem}" BaseRecalibrator \
          -I {input.bam} -R {input.ref} {params.known_arg} -O {output.table} 2> {log}
        """

rule applybqsr:
    input:
        bam   = f"{D_MARKDUP}/{{sample}}_marked_duplicates.bam",
        table = f"{D_BQSR}/{{sample}}_recalibration.table",
        ref   = REF,
    output:
        bam = f"{D_BQSR}/{{sample}}_bqsr.bam",
    log:
        f"{D_BQSR}/logs/{{sample}}_applybqsr.log",
    params:
        mem = JAVA_MEM,
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.mem}" ApplyBQSR \
          -R {input.ref} -I {input.bam} --bqsr-recal-file {input.table} -O {output.bam} 2> {log}
        """
