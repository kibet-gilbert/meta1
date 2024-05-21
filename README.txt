Analysis of metagenomics data (shotgun)
Steps:
step1: Project dir set up:
	mkdir -p ./results{fastqc,fastp,kraken2,metaspades,maxbin,checkm,krona,pavian}
step2: Environment set-up
	- mamba installation
	- packages: from Bioconda
code: 	mamba create -n meta1
	mamba install -c bioconda fastp fastqc kraken2 checkm maxbin spades
