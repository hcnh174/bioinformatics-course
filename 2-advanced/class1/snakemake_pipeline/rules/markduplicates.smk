# 2-3 MarkDuplicates — PCR 重複マーキング / mark duplicates ( per sample )
rule markduplicates:
    input:
        bam = f"{D_BWA}/{{sample}}.bam",
    output:
        bam     = f"{D_MARKDUP}/{{sample}}_marked_duplicates.bam",
        bai     = f"{D_MARKDUP}/{{sample}}_marked_duplicates.bai",
        metrics = f"{D_MARKDUP}/{{sample}}_marked_dup_metrics.txt",
    log:
        f"{D_MARKDUP}/logs/{{sample}}.log",
    params:
        mem = JAVA_MEM,
    conda:
        "../envs/picard.yaml"
    shell:
        r"""
        picard {params.mem} MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} 2> {log}
        picard {params.mem} BuildBamIndex -INPUT {output.bam} -OUTPUT {output.bai} 2>> {log}
        """
