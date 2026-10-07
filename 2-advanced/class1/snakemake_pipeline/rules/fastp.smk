# 2-1 fastp — 品質トリミング / adapter & quality trimming ( per sample )
#   入力 FASTQ は samples.tsv 由来（fastq_input 関数）/ inputs come from the TSV
rule fastp:
    input:
        unpack(fastq_input),
    output:
        t1   = f"{D_FASTP}/{{sample}}_1_trimmed.fastq.gz",
        t2   = f"{D_FASTP}/{{sample}}_2_trimmed.fastq.gz",
        u1   = f"{D_FASTP}/{{sample}}_1_unpaired.fastq.gz",
        u2   = f"{D_FASTP}/{{sample}}_2_unpaired.fastq.gz",
        html = f"{D_FASTP}/{{sample}}_fastp.html",
        json = f"{D_FASTP}/{{sample}}_fastp.json",
    log:
        f"{D_FASTP}/logs/{{sample}}.log",
    params:
        n = FP["n"], l = FP["l"], q = FP["q"],
        t = FP["t"], T = FP["T"], w = FP["w"],
    threads: FP["w"]
    conda:
        "../envs/fastp.yaml"
    shell:
        r"""
        fastp --in1 {input.r1} --in2 {input.r2} -3 \
          --out1 {output.t1} --out2 {output.t2} \
          --unpaired1 {output.u1} --unpaired2 {output.u2} \
          -h {output.html} -j {output.json} \
          -n {params.n} -l {params.l} -q {params.q} -t {params.t} -T {params.T} -w {params.w} \
          2> {log}
        """
