#!/usr/bin/nexflow env

# Activate mamba environment:
#source ~/bioinformatics/anaconda3/etc/profile.d/conda.sh
#mamba activate meta1

# Fastqc:
mkdir -p ./results/{fastqc,fastp,kraken2,metaspades,maxbin,checkm,krona,pavian}
fastqc -t 4 \
	-o ./results/fastqc/ \
	./data/fastq/SRR28624259_R1.fastq.gz \
	./data/fastq/SRR28624259_R2.fastq.gz

## Fastp:
fastp --in1 ./data/fastq/SRR28624259_R1.fastq.gz \
	--in2 ./data/fastq/SRR28624259_R2.fastq.gz \
	--out1 ./results/fastp/SRR28624259-trim_R1.fastq.gz \
	--out2 ./results/fastp/SRR28624259-trim_R2.fastq.gz \
	--html ./results/fastp/SRR28624259.fastp.html \
	--failed_out ./results/fastp/SRR28624259_failed.fastq.gz \
	--thread 4 \
	--detect_adapter_for_pe \
	--dedup \
	|& tee ./results/fastp/sample01.fastp.log

# Kraken2 (host genome):
kraken2 --db ./data/databases/kraken2_human_db/ \
	--threads 4 \
	--unclassified-out ./results/kraken2/SRR28624259.nohost#.fastq \
	--classified-out ./results/kraken2/SRR28624259.host#.fastq \
	--report ./results/kraken2/SRR28624259.kraken2.report.txt \
	--output ./results/kraken2/SRR28624259.kraken2.out \
	--gzip-compressed \
	--report-zero-counts \
	--paired ./results/fastp/SRR28624259-trim_R1.fastq.gz \
	./results/fastp/SRR28624259-trim_R2.fastq.gz

## kraken2 (taxonomic classification):
kraken2 -db ./data/databases/k2_viral_20240112/ \
	--threads 4 \
	--unclassified-out ./results/kraken2/SRR28624259.unclassified#.fastq \
	--classified-out ./results/kraken2/SRR28624259.classified#.fastq \
	--report ./results/kraken2/SRR28624259_tax_kreport.txt \
	--output ./results/kraken2/SRR28624259_tax_kraken2.out \
	--report-zero-counts \
	--paired ./results/kraken2/SRR28624259.nohost_1.fastq \
	./results/kraken2/SRR28624259.nohost_2.fastq

# kronatools:
# Download taxonomy
mkdir -p data/databases/kronaDB/taxonomy
wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O ./data/databases/kronaDB/taxdump.tar.gz
ktUpdateTaxonomy.sh --only-build ./data/databases/kronaDB/taxonomy
cat ./results/kraken2/SRR28624259_tax_kraken2.out | cut -f 2,3 > ./results/krona/SRR28624259_kraken2.krona
ktImportTaxonomy -tax ./data/databases/kronaDB/taxonomy \
	-o ./data/krona/SRR28624259_taxonomy.krona.html \
	./results/krona/SRR28624259_kraken2.krona


# metaspades:
spades.py --isolate \
	--threads 4 \
	--memory 8 \
	-1 ./results/kraken2/SRR28624259.classified_1.fastq \
	-2 ./results/kraken2/SRR28624259.classified_2.fastq \
	-o ./results/metaspades
        mv ./results/metaspades/assembly_graph_with_scaffolds.gfa ./results/metaspades/SPAdes-SRR28624259_graph.gfa
        mv ./results/metaspades/scaffolds.fasta ./results/metaspades/SPAdes-SRR28624259_scaffolds.fasta
        mv ./results/metaspades/contigs.fasta ./results/metaspades/SPAdes-SRR28624259_contigs.fasta
        mv ./results/metaspades/spades.log ./results/metaspades/SPAdes-SRR28624259.log

## kraken2 (contigs classification):
kraken2 -db ./data/databases/k2_viral_20240112/ \
	--threads 4 \
	--report ./results/kraken2/SRR28624259_taxContigs_kreport.txt \
	--output ./results/kraken2/SRR28624259_taxContigs_kraken2.out \
	--report-zero-counts \
	./results/metaspades/SPAdes-SRR28624259_contigs.fasta

### visualization:
## kronatools:
cat ./results/kraken2/SRR28624259_taxContigs_kraken2.out | cut -f 2,3 > ./results/krona/SRR28624259_taxContigs_kraken2.krona
ktImportTaxonomy -tax ./data/databases/kronaDB/taxonomy \
	-o ./results/krona/SRR28624259_taxContigs_taxonomy.krona.html \
	./results/krona/SRR28624259_taxContigs_kraken2.krona

