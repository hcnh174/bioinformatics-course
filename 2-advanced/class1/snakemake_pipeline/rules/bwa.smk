# 2-2 BWA — アライメント + 整列 + 索引 / align, sort, index ( per sample )
rule bwa_map:
    input:
        r1  = f"{D_FASTP}/{{sample}}_1_trimmed.fastq.gz",
        r2  = f"{D_FASTP}/{{sample}}_2_trimmed.fastq.gz",
        idx = f"{BWA_INDEX}.bwt",
    output:
        bam = f"{D_BWA}/{{sample}}.bam",
        bai = f"{D_BWA}/{{sample}}.bam.bai",
    log:
        f"{D_BWA}/logs/{{sample}}.log",
    params:
        rg    = r"@RG\tID:{sample}\tPL:ILLUMINA\tSM:{sample}",
        index = BWA_INDEX,
    threads: THREADS
    conda:
        "../envs/bwa.yaml"
    shell:
        r"""
        bwa mem -t {threads} -R "{params.rg}" {params.index} {input.r1} {input.r2} 2> {log} \
          | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """
