#!/usr/bin/bash env

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

# Fastqc:
mkdir 
fastqc -t 4 \
	-o ./data/fastqc/ \
	./data/fastq/sample01_R1.fastq.gz \
	./data/fastq/sample01_R2.fastq.gz

# Fastp:
fastp --in1 ./data/fastq/sample01_R1.fastq.gz \
	--in2 ./data/fastq/sample01_R2.fastq.gz \
	--out1 ./data/fastp/sample01_R1.trim.fastq.gz \
	--out2 ./data/fastp/sample01_R2.trim.fastq.gz \
	--json ./data/fastp/sample01.fastp.json \
	--html ./data/fastp/sample01.fastp.html \
	--failed_out ./data/fastp/sample01_fail.fastq.gz \
	--thread 4 \
	-5 -3 -r \
	--detect_adapter_for_pe \
	--qualified_quality_phred 20 \
	--cut_mean_quality 20 \
	--length_required 15 \
	--dedup \
	|& tee ./data/fastp/sample01.fastp.log

# Fastqc:
fastqc -t 4 \
	-o ./data/fastqc/ \
	./data/fastq/sample01_R1.fastq.gz \
	./data/fastq/sample01_R2.fastq.gz

# Kraken2 (host genome):
kraken2 -db ./data/database/host_db/kraken2_human_db \
	--threads 4 \
	--unclassified-out ./data/kraken/sample01.unclassified#.fastq \
	--classified-out ./data/kraken/sample01.classified#.fastq \
	--report ./data/kraken/sample01.kraken2.report.txt \
	--output ./data/kraken/sample01.kraken2.out \
	--gzip-compressed \
	--report-zero-counts \
	--paired ./data/fastp/sample01_R1.trim.fastq.gz \
	./data/fastp/sample01_R2.trim.fastq.gz

# kraken2 (taxonomic classification):
kraken2 -db ./data/database/kraken2/k2_pluspf_16gb_20221209 \
	--threads 4 \
	--unclassified-out ./data/kraken/sample01.all_unclassified#.fastq \
	--classified-out ./data/kraken/sample01.all_classified#.fastq \
	--report ./data/kraken/sample01_kreport.txt \
	--output ./data/kraken/sample01_kraken2.out \
	--bzip2-compressed \
	--report-zero-counts \
	--paired ./data/kraken/sample01.unclassified_1.fastq.bz2 \
	./data/kraken/sample01.unclassified_2.fastq.bz2

# metaspades:
metaspades.py \
            $args \
            --threads "${task.cpus}" \
            --memory $maxmem \
            --pe1-1 ${reads[0]} \
            --pe1-2 ${reads[1]} \
            -o spades
        mv spades/assembly_graph_with_scaffolds.gfa SPAdes-${meta.id}_graph.gfa
        mv spades/scaffolds.fasta SPAdes-${meta.id}_scaffolds.fasta
        mv spades/contigs.fasta SPAdes-${meta.id}_contigs.fasta
        mv spades/spades.log SPAdes-${meta.id}.log
        gzip "SPAdes-${meta.id}_contigs.fasta"
        gzip "SPAdes-${meta.id}_graph.gfa"
        gzip -c "SPAdes-${meta.id}_scaffolds.fasta" > "SPAdes-${meta.id}_scaffolds.fasta.gz"
# maxbin2:
#mkdir input/ && mv $contigs input/
#    run_MaxBin.pl \\
#        -contig input/$contigs \\
#        $associate_files \\
#        -thread $task.cpus \\
#        $args \\
#        -out $prefix
#
#    gzip *.fasta *.noclass *.tooshort *log *.marker

# metabat2:
metabat2 \\
        $args \\
        -i $fasta \\
        $depth_file \\
        -t $task.cpus \\
        --saveCls \\
        -o ${prefix}

# Das tool (Bin-refinement Optional):
DAS_Tool \\
        $args \\
        $proteins_pred \\
        $db_dir \\
        -t $task.cpus \\
        -i $bin_list \\
        -c $clean_contigs \\
        -o $prefix
# checkm:
checkm \\
        qa \\
        --threads ${task.cpus} \\
        --file ${prefix}.${suffix} \\
        $marker_file \\
        $analysis_dir \\
        $coverage \\
        $exclude \\
        $args
# kraken2:
## visualization:
# kronatools:
ktImportTaxonomy -tax ./data/database/krona/taxonomy \
	-o ./data/krona/sample01_taxonomy.krona.html \
	./data/krona/sample01-kraken2.krona

# pavian:

## Optional:
# bbtools:
# CAT: taxonomic classification of bins after binning
#
