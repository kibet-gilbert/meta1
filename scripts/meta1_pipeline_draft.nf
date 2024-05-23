#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads='/home/gkibet/bioinformatics/training/meta1/data/fastq/test-data/*_R{1,2}.fastq.gz'
//params.reads='/home/gkibet/bioinformatics/training/meta1/data/fastq/raw-data/*_L001_R{1,2}_001.fastq.gz'
reads = Channel.fromFilePairs(params.reads, checkIfExists:true)
println(reads.view())
params.kraken2_db_human = '/home/gkibet/bioinformatics/training/meta1/data/databases/kraken2_human_db/'
kraken2_db_human = Channel.fromPath(params.kraken2_db_human, checkIfExists:true)

params.kraken2_db_viral = './data/databases/k2_viral_20240112/'
kraken2_db_viral = Channel.fromPath(params.kraken2_db_viral, checkIfExists:true)

params.taxonomy='/home/gkibet/bioinformatics/training/meta1/data/databases/kronaDB/taxonomy/'
taxonomy = Channel.fromPath(params.taxonomy, checkIfExists:true)




// params.krona_db_taxonomy = './data/databases/kronaDB/taxonomy'
// params.outputDir = './results'

process FastQC {
    // publishDir "${params.outputDir}/fastqc", mode: 'copy'
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"), emit: zip

    script:
    """
    fastqc -t $task.cpus -o . ${reads}
    """
}

process FastP {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*-trim*"), emit: trimmed_reads
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.log"), emit: log


    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${sample_id}-trim_R1.fastq.gz \\
        --out2 ${sample_id}-trim_R2.fastq.gz \\
        --html ${sample_id}.fastp.html \\
        --failed_out ${sample_id}_failed.fastq.gz \\
        --thread $task.cpus \\
        --detect_adapter_for_pe \\
        --dedup \\
        |& tee ${sample_id}.fastp.log
    """
}

process Kraken2Host {
    // publishDir "${params.outputDir}/kraken2", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(trimmed_reads)
    path(kraken2_db_human)

    output:
    tuple val(sample_id), path("*.nohost{.,_}*"), emit: nohost_reads
    tuple val(sample_id), path("*.report.txt"), emit: kraken_hostmap_report1
    tuple val(sample_id), path("*.kraken2.out"), emit: kraken_hostmap_report2

    script:
    """
    kraken2 --db ${kraken2_db_human} --threads 4 \\
            --unclassified-out ${sample_id}.nohost#.fastq \\
            --classified-out ${sample_id}.host#.fastq \\
            --report ${sample_id}.kraken2.report.txt \\
            --output ${sample_id}.kraken2.out \\
            --gzip-compressed --report-zero-counts \\
            --paired $trimmed_reads
    """
}

process Kraken2Taxonomy {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(nohost_reads)
    path(kraken2_db_viral)

    output:
    tuple val(sample_id), path("*.classified{.,_}*"), emit: classified_reads
    tuple val(sample_id), path("*_kreport.txt"), emit: kraken_viral_report1
    tuple val(sample_id), path("*_kraken2.out"), emit: kraken_viral_report2

    script:
    """
    kraken2 --db ${kraken2_db_viral} --threads 4 \\
            --unclassified-out ${sample_id}.unclassified#.fastq \\
            --classified-out ${sample_id}.classified#.fastq \\
            --report ${sample_id}_tax_kreport.txt \\
            --output ${sample_id}_tax_kraken2.out \\
            --report-zero-counts --paired $nohost_reads
    """
}

process KronaTools {
    tag "$sample_id"  
    
    input:
    tuple val(sample_id), path(kraken_viral_report2)
    path(taxonomy)
    
    output:
    tuple val(sample_id), path("*.html"), emit: krona_html
    // path "${params.outputDir}/krona/SRR28624259_taxonomy.krona.html"

    script:
    //wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O ${params.krona_db_taxonomy}/taxdump.tar.gz
    //ktUpdateTaxonomy.sh --only-build ${params.krona_db_taxonomy}
    """    
    cat $kraken_viral_report2 | cut -f 2,3 > ${sample_id}_kraken2.krona
    ktImportTaxonomy -tax ${taxonomy} \\
                     -o ${sample_id}_taxonomy.krona.html \\
                     ${sample_id}_kraken2.krona
    """
}

process MetaSPAdes {
    tag "$sample_id"     
    
    input:
    tuple val(sample_id), path(nohost_reads)

    output:
    tuple val(sample_id), path("*_contigs.fasta"), emit: contigs
    tuple val(sample_id), path("*_scaffolds.fasta"), emit: scaffolds
    tuple val(sample_id), path("*_graph.gfa"), emit: graphgfa
    tuple val(sample_id), path("*.log"), emit: spadeslog

    script:
    """
    spades.py \\
        --meta \\
        --threads 4 \\
        --memory 8 \\
        -1 ${nohost_reads[0]} \\
        -2 ${nohost_reads[1]} \\
        -o .
    mv assembly_graph_with_scaffolds.gfa SPAdes-${sample_id}_graph.gfa
    mv scaffolds.fasta SPAdes-${sample_id}_scaffolds.fasta
    mv contigs.fasta SPAdes-${sample_id}_contigs.fasta
    mv spades.log SPAdes-${sample_id}.log
    """
}

process Kraken2Contigs {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(contigs)
    path(kraken2_db_viral)

    output:
    tuple val(sample_id), path("*_kreport.txt"), emit: kraken_viral_report1
    tuple val(sample_id), path("*_kraken2.out"), emit: kraken_viral_report2

    script:
    """
    kraken2 --db ${kraken2_db_viral} --threads 4 \\
            --report ${sample_id}_taxContigs_kreport.txt \\
            --output ${sample_id}_taxContigs_kraken2.out \\
            --report-zero-counts $contigs
    """
}

process KronaContigs {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(kraken_viral_report2)
    path(taxonomy)

    output:
    tuple val(sample_id), path("*.html"), emit: kronacontigs_html

    script:
    """
    cat $kraken_viral_report2 | cut -f 2,3 > ${sample_id}_taxContigs_kraken2.krona
    ktImportTaxonomy -tax ${taxonomy} \
                     -o ${sample_id}_taxContigs_taxonomy.krona.html \
                     ${sample_id}_taxContigs_kraken2.krona
    """
}

// Workflow definition
workflow {
    FastQC(reads)
    FastP(reads)
    Kraken2Host(FastP.out.trimmed_reads, kraken2_db_human)
    Kraken2Taxonomy(Kraken2Host.out.nohost_reads, kraken2_db_viral)
    KronaTools(Kraken2Taxonomy.out.kraken_viral_report2, taxonomy)
    MetaSPAdes(Kraken2Host.out.nohost_reads)
    // println(Kraken2Host.out.nohost_reads.view())
    // println(Kraken2Host.out.nohost_reads.collect().view())
    Kraken2Contigs(MetaSPAdes.out.contigs,kraken2_db_viral)
    KronaContigs(Kraken2Contigs.out.kraken_viral_report2,taxonomy)
}

