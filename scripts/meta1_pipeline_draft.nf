#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads='/home/gkibet/bioinformatics/training/meta1/data/fastq/test-data/*_R{1,2}.fastq.gz'
reads = Channel.fromFilePairs(params.reads, checkIfExists:true)
println(reads.view())

// params.kraken2_db_human = './data/databases/kraken2_human_db/'
// params.kraken2_db_viral = './data/databases/k2_viral_20240112/'
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
        --in1 $reads[1] \\
        --in2 $reads[2] \\
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

// process Kraken2Host {
//     publishDir "${params.outputDir}/kraken2", mode: 'copy'
    
//     input:
//     path trimmed_fastq1 from FastP.out[0]
//     path trimmed_fastq2 from FastP.out[1]

//     output:
//     path "${params.outputDir}/kraken2/SRR28624259.nohost_1.fastq"
//     path "${params.outputDir}/kraken2/SRR28624259.nohost_2.fastq"

//     script:
//     """
//     kraken2 --db ${params.kraken2_db_human} --threads 4 \
//             --unclassified-out ${params.outputDir}/kraken2/SRR28624259.nohost#.fastq \
//             --classified-out ${params.outputDir}/kraken2/SRR28624259.host#.fastq \
//             --report ${params.outputDir}/kraken2/SRR28624259.kraken2.report.txt \
//             --output ${params.outputDir}/kraken2/SRR28624259.kraken2.out \
//             --gzip-compressed --report-zero-counts \
//             --paired $trimmed_fastq1 $trimmed_fastq2
//     """
// }

// process Kraken2Taxonomy {
//     publishDir "${params.outputDir}/kraken2", mode: 'copy'
    
//     input:
//     path nohost_fastq1 from Kraken2Host.out[0]
//     path nohost_fastq2 from Kraken2Host.out[1]

//     output:
//     path "${params.outputDir}/kraken2/SRR28624259.classified_1.fastq"
//     path "${params.outputDir}/kraken2/SRR28624259.classified_2.fastq"

//     script:
//     """
//     kraken2 --db ${params.kraken2_db_viral} --threads 4 \
//             --unclassified-out ${params.outputDir}/kraken2/SRR28624259.unclassified#.fastq \
//             --classified-out ${params.outputDir}/kraken2/SRR28624259.classified#.fastq \
//             --report ${params.outputDir}/kraken2/SRR28624259_tax_kreport.txt \
//             --output ${params.outputDir}/kraken2/SRR28624259_tax_kraken2.out \
//             --report-zero-counts --paired $nohost_fastq1 $nohost_fastq2
//     """
// }

// process KronaTools {
//     publishDir "${params.outputDir}/krona", mode: 'copy'
    
//     input:
//     path kraken2_output from Kraken2Taxonomy.out[2]

//     output:
//     path "${params.outputDir}/krona/SRR28624259_taxonomy.krona.html"

//     script:
//     """
//     mkdir -p ${params.krona_db_taxonomy}
//     wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O ${params.krona_db_taxonomy}/taxdump.tar.gz
//     ktUpdateTaxonomy.sh --only-build ${params.krona_db_taxonomy}
//     cat $kraken2_output | cut -f 2,3 > ${params.outputDir}/krona/SRR28624259_kraken2.krona
//     ktImportTaxonomy -tax ${params.krona_db_taxonomy} \
//                      -o ${params.outputDir}/krona/SRR28624259_taxonomy.krona.html \
//                      ${params.outputDir}/krona/SRR28624259_kraken2.krona
//     """
// }

// process MetaSPAdes {
//     publishDir "${params.outputDir}/metaspades", mode: 'copy'
    
//     input:
//     path classified_fastq1 from Kraken2Taxonomy.out[0]
//     path classified_fastq2 from Kraken2Taxonomy.out[1]

//     output:
//     path "${params.outputDir}/metaspades/"

//     script:
//     """
//     spades.py --isolate --threads 4 --memory 8 \
//               -1 $classified_fastq1 -2 $classified_fastq2 \
//               -o ${params.outputDir}/metaspades
//     mv ${params.outputDir}/metaspades/assembly_graph_with_scaffolds.gfa ${params.outputDir}/metaspades/SPAdes-SRR28624259_graph.gfa
//     mv ${params.outputDir}/metaspades/scaffolds.fasta ${params.outputDir}/metaspades/SPAdes-SRR28624259_scaffolds.fasta
//     mv ${params.outputDir}/metaspades/contigs.fasta ${params.outputDir}/metaspades/SPAdes-SRR28624259_contigs.fasta
//     mv ${params.outputDir}/metaspades/spades.log ${params.outputDir}/metaspades/SPAdes-SRR28624259.log
//     """
// }

// process Kraken2Contigs {
//     publishDir "${params.outputDir}/kraken2", mode: 'copy'
    
//     input:
//     path contigs_fasta from MetaSPAdes.out[2]

//     output:
//     path "${params.outputDir}/kraken2/SRR28624259_taxContigs_kraken2.out"

//     script:
//     """
//     kraken2 --db ${params.kraken2_db_viral} --threads 4 \
//             --report ${params.outputDir}/kraken2/SRR28624259_taxContigs_kreport.txt \
//             --output ${params.outputDir}/kraken2/SRR28624259_taxContigs_kraken2.out \
//             --report-zero-counts $contigs_fasta
//     """
// }

// process KronaContigs {
//     publishDir "${params.outputDir}/krona", mode: 'copy'
    
//     input:
//     path kraken2_contigs_output from Kraken2Contigs.out

//     output:
//     path "${params.outputDir}/krona/SRR28624259_taxContigs_taxonomy.krona.html"

//     script:
//     """
//     cat $kraken2_contigs_output | cut -f 2,3 > ${params.outputDir}/krona/SRR28624259_taxContigs_kraken2.krona
//     ktImportTaxonomy -tax ${params.krona_db_taxonomy} \
//                      -o ${params.outputDir}/krona/SRR28624259_taxContigs_taxonomy.krona.html \
//                      ${params.outputDir}/krona/SRR28624259_taxContigs_kraken2.krona
//     """
// }

// Workflow definition
workflow {
    FastQC(reads)
    FastP(reads)
    // Kraken2Host()
    // Kraken2Taxonomy()
    // KronaTools()
    // MetaSPAdes()
    // Kraken2Contigs()
    // KronaContigs()
}

