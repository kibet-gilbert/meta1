# Meta1 [Metagenomics Diagnosis of Clinical Specimens]
## Analysis of shotgun metagenomics data from clinical samples for pathogen diagnosis

### attributes
- pipeline name: meta1
- workflow manager: nextflow
- environment manager: conda/singularity

### workflow outline
- fastqc [read QC checks]
- fastp [read trimming]
- kraken2 [dehosting, read classification and contig classification]
- spades [de novo assembly]
- krona [biom tables vizualization]

### step1: Project dir development & set up:
```
mkdir -p ./data/{fastq,databases}

mkdir -p ./results/{fastqc,fastp,kraken2,metaspades,maxbin,checkm,krona}                                    
```

### step2: Environment set-up
```
mamba create env -f meta1_env.yml

mamba activate meta1
```

### step3: Setup data and databases
`
mkdir -p ./data/fastq/raw_data && cd ./data/fastq/raw_data/
`

- download the raw data from the following gdrive link: "https://shorturl.at/nRDJg"
- unzip and move it to the "./data/fastq/raw_data/" folder

```
cd ../../databases/

wget https://ndownloader.figshare.com/files/23567780 -O kraken2_human_db.tar.gz

wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz 

tar -xvf *.tar.gz

cd kronaDB/

wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

ktUpdateTaxonomy.sh --only-build ./taxonomy
```

### step4: Running the pipeline
- create bash script step-by-step and test (scripts/meta1_pipeline.sh)
- Test one-liner commands and add that which works to the script
- convert the bash script to a nextflow workflow (scripts/meta1_draft.nf)
- Step-by-step create processes based on the bash script commands/steps
- For each process provide: $tag, input channels, output channels, script
- For each process provide the directives in the config file (nextflow.config)
- Also provide params: declare default params

### step5: package mamba environment to .yaml
`
mamba env export -n meta1 > meta1.yaml
`

### step6: create a singularity/apptainer definition file (scripts/meta1_apptainer.def)
- Design the definition file paying attention to sections: https://apptainer.org/docs/ user/main/definition_files.html#sections

### step7: build an singularity/apptainer image
`
apptainer build ./scripts/meta1_apptainer.sif ./scripts/meta1_apptainer.def
`

### step8: Deploy the image with different test data:
`
nextflow run ./scripts/meta1.nf -with-apptainer ./scripts/meta1_apptainer.sif -resume --reads='data/fastq/raw-data/*_L001_R{1,2}_001.fastq.gz'
`

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