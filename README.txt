Analysis of metagenomics data (shotgun)
Steps:
step1: Project dir set up:
	mkdir -p ./results{fastqc,fastp,kraken2,metaspades,maxbin,checkm,krona,pavian}

step2: Environment set-up
	- mamba installation:
	- packages: from Bioconda
code:
	$ mamba create -n meta1
	$ mamba install -c bioconda fastp fastqc kraken2 krona checkm maxbini2 metabat2 megahit spades cat das_tool

step3: Data retrival - fastq, kraken2 databases (human/viral), checkm database, krona database

step4: create bash script step-by-step and test (scripts/meta1_pipeline.sh)

step5: package mamba environment to .yaml
	$ mamba env export -n meta1 > meta1_environment.yaml

step6: create a singularity/apptainer definition file (scripts/meta1_apptainer.def)

step7: build an singularity/apptainer image
	$ apptainer build ./scripts/meta1_apptainer.sif ./scripts/meta1_apptainer.def

step8: Deploy the image with different test data:
	$ ~/Downloads/Programs/nextflow run scripts/meta1_pipeline_draft.nf -with-apptainer ./meta1_apptainer.sif -resume
	$ ~/Downloads/Programs/nextflow run scripts/meta1_pipeline_draft.nf -with-apptainer ./meta1_apptainer.sif -resume --reads='data/fastq/raw-data/*_L001_R{1,2}_001.fastq.gz'
	$ ~/Downloads/Programs/nextflow run scripts/meta1_pipeline_draft.nf -with-apptainer ./meta1_apptainer.sif -resume --reads='data/fastq/raw-data2/*_R{1,2}.fastq.gz'

Comments: Notable challenges or errors
- 

