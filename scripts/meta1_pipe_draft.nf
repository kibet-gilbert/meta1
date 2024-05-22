nextflow.enable.dsl=2

process fastqc {
    tag 'FastQC'

    input:
    tuple val(sample), path(reads)

    output:
    path "results/fastqc/*"

    script:
    """
    mkdir -p ./results/fastqc
    fastqc -t 4 -o ./results/fastqc/ $reads
    """
}

process fastp {
    tag 'Fastp'

    input:
    tuple val(sample), path(reads)

    output:
    path "results/fastp/*"

    script:
    """
    mkdir -p ./results/fastp
    fastp --in1 $reads[0] \
          --in2 $reads[1] \
          --out1 ./results/fastp/${sample}-trim_R1.fastq.gz \
          --out2 ./results/fastp/${sample}-trim_R2.fastq.gz \
          --html ./results/fastp/${sample}.fastp.html \
          --failed_out ./results/fastp/${sample}_failed.fastq.gz \
          --thread 4 \
          --detect_adapter_for_pe \
          --dedup \
          |& tee ./results/fastp/${sample}.fastp.log
    """
}

workflow {
    samples = Channel.of(
        tuple("SRR28624259", ["./data/fastq/SRR28624259_R1.fastq.gz", "./data/fastq/SRR28624259_R2.fastq.gz"])
    )

    samples |>
    fastqc |>
    fastp
}

