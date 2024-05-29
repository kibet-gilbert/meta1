Analysis of viral metagenomics shotgun data
-------------------------------------------

pipeline name: magviral
workflow manager: nextflow
environment manager: conda/singularity(apptainer)

Steps:
-----
step1: Initiated a github repo: https://github.com/kibet-gilbert/meta1/
	- make a project directory: 'meta1'
	$ mkdir meta1 && cd meta1
	$ git init
	- create a README.md
	$ touch README.md
	- add, commit and push
	$ git add -A
	$ git commit -m "Initial commit"
	$ git push

step2: Project dir development & set up:
	$ mkdir -p ./data/{fastq,databases}
	$ mkdir -p ./results/{fastqc,fastp,kraken2,metaspades,maxbin,checkm,krona,pavian}

step3: Environment set-up
	- mamba installation: install miniconda and then install mamba:
	$ conda install mamba
	- install packages from Bioconda
	$ mamba create -n magviral
	$ mamba activate magviral
	$ mamba install -c bioconda fastp fastqc kraken2 krona checkm maxbin2 metabat2 megahit spades cat das_tool

step4: Data retrival - fastq, kraken2 databases (human/viral), checkm database, krona database
	$ mkdir -p ./data/fastq/raw_data && cd ./data/fastq/raw_data/
	$ wget https://github.com/kibet-gilbert/meta1/raw/main/data/fastq/raw-data/MTW-364_S2_L001_R1_001.fastq.gz
	$ wget https://github.com/kibet-gilbert/meta1/raw/main/data/fastq/raw-data/MTW-364_S2_L001_R2_001.fastq.gz
	$ wget https://github.com/kibet-gilbert/meta1/raw/main/data/fastq/raw-data/MTW-4167_S1_L001_R1_001.fastq.gz
	$ wget https://github.com/kibet-gilbert/meta1/raw/main/data/fastq/raw-data/MTW-4167_S1_L001_R2_001.fastq.gz
	$ mkdir -p ../test_data && cd ../test_data/
	$ wget https://github.com/kibet-gilbert/meta1/raw/main/data/fastq/raw-data/SRR28624259_R1.fastq.gz
	$ wget https://github.com/kibet-gilbert/meta1/raw/main/data/fastq/raw-data/SRR28624259_R2.fastq.gz
	$ cd ../../databases/
	$ wget https://ndownloader.figshare.com/files/23567780 -O kraken2_human_db.tar.gz
	$ wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz 
	$ wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
	$ tar -xvf *.tar.gz
	$ cd kronaDB/
	$ wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
	$ ktUpdateTaxonomy.sh --only-build ./taxonomy
	$ cd ../../../

step5: create bash script step-by-step and test (scripts/meta1_pipeline.sh)
	- Test one-liner commands and add that which works to the script

step6: convert the bash script to a nextflow workflow (scripts/magviral_draft.nf)
	- Step-by-step create processes based on the bash script commands/steps
	- For each process provide: $tag, input channels, output channels, script
	- For each process provide the directives in the config file (nextflow.config)
	- Also provide params: declare default params

step7: package mamba environment to .yaml
	$ mamba env export -n magviral > magviral.yaml

step8: create a singularity/apptainer definition file (scripts/magviral_apptainer.def)
	- Design the definition file paying attention to sections: https://apptainer.org/docs/user/main/definition_files.html#sections

step9: build an singularity/apptainer image
	$ apptainer build ./scripts/magviral_apptainer.sif ./scripts/magviral_apptainer.def

step10: Deploy the image with different test data:
	$ ~/Downloads/Programs/nextflow run ./scripts/magviral.nf -with-apptainer ./scripts/magviral_apptainer.sif -resume
	$ ~/Downloads/Programs/nextflow run ./scripts/magviral.nf -with-apptainer ./scripts/magviral_apptainer.sif -resume --reads='data/fastq/raw-data/*_L001_R{1,2}_001.fastq.gz'
	$ ~/Downloads/Programs/nextflow run ./scripts/magviral.nf -with-apptainer ./scripts/magviral_apptainer.sif -resume --reads='data/fastq/raw-data2/*_R{1,2}.fastq.gz'
	$ ~/Downloads/Programs/nextflow run ./scripts/magviral.nf -with-apptainer ./scripts/magviral_apptainer.sif -resume --reads='data/fastq/hev_data/*_S*_con_R{1,2}_001.fastq.gz'

Comments: Notable challenges or errors
--------------------------------------
- Dataset sizes:
	- Fastq file sizes in metagenomics are large (>500Mb) - We eneded up using smaller sample isolate shotgun data
	- Classification databases are large and needs fast internet speeds - used small viral-only database 

- Some processes fail due to different reasons:
	- Binning fails where contigs are too few, or taxonomic diverity in the data is too low
	- krona visualization can fail due to "too many taxIds"
	- kraken2 filteration of host reads fails though very rarely due to unknown reasons***

- Some processes take so much resources (memory/time)
	- spades assembly - MAG (metagenomic assembly of Genomes) needs lots of RAM and takes long made more challenging by fastq file size.
	- kraken2 taxonomic classification can also take a lot of RAM where the database is large. memory mapping can help resolve it but is slow

