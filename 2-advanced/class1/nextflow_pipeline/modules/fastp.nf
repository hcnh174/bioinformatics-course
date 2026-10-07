// 2-1 fastp — 品質トリミング / adapter & quality trimming ( per sample )
process FASTP {
    tag "${meta.id}"
    conda "${projectDir}/envs/fastp.yaml"
    publishDir "${params.outdir}/01fastp", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_1_trimmed.fastq.gz"), path("${meta.id}_2_trimmed.fastq.gz"), emit: trimmed
    path "${meta.id}_fastp.html"
    path "${meta.id}_fastp.json"

    script:
    """
    fastp --in1 ${r1} --in2 ${r2} -3 \\
      --out1 ${meta.id}_1_trimmed.fastq.gz --out2 ${meta.id}_2_trimmed.fastq.gz \\
      --unpaired1 ${meta.id}_1_unpaired.fastq.gz --unpaired2 ${meta.id}_2_unpaired.fastq.gz \\
      -h ${meta.id}_fastp.html -j ${meta.id}_fastp.json \\
      -n ${params.fastp.n} -l ${params.fastp.l} -q ${params.fastp.q} \\
      -t ${params.fastp.t} -T ${params.fastp.T} -w ${params.fastp.w}
    """

    stub:
    """
    touch ${meta.id}_1_trimmed.fastq.gz ${meta.id}_2_trimmed.fastq.gz \\
          ${meta.id}_fastp.html ${meta.id}_fastp.json
    """
}
