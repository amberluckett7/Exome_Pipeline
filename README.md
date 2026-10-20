# Exome pipeline 

Author: **Amber Luckett** 

Email: **A.Luckett3@exeter.ac.uk**

This exome pipeline is developed from R04 Genome pipeline developed by James Russ-Silsby (https://github.com/JamesR-S/Genome_Pipeline).
The pipeline also runs in nextflow and makes use of containers for many of the processing steps.

This guide will go through the set up for the pipeline and how to run it.

## Setting up the pipeline

### Installing Dependencies

The easiest way to install the dependencies is through micromamba. This is **strongly recommended** (over conda) on ISCA, as conda may interfere with centrally managed conda installations.

This can be installed with the following command:
```
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

Once micromamba has been installed and loaded, you will need the `environment.yml` file before creating the environment.
You can download and install the environment with the following commands:
```
wget https://raw.githubusercontent.com/amberluckett7/Exome_Pipeline/main/environment.yml
micromamba env create -f environment.yml
```

If you encounter issues downloading the file, or if the repo is private, clone the repository instead:
```
git clone git@github.com:amberluckett7/Exome_Pipeline.git
cd Exome_Pipeline
micromamba env create -f environment.yml
```

Alternatively, if you need to manually specify packages, you can run (note that this is not recommended unless you have issues with the environment.yml file):
```
micromamba create \
	-p exome_pipeline_env \
	-c conda-forge \
	-c bioconda \
	conda-forge::apptainer=1.4.1 \
	mitosalt \
	bioconda::nextflow=25.4.2
```

To activate the environment, you then run:
```
micromamba activate exome_pipeline_env
```

### Setting up environment

In addition to installing the dependencies you will also need to add the `APPTAINER_TMPDIR` and `APPTAINER_CACHEDIR` variable to your environment. This is the directory to which temporary files that are part of the container download process will be saved. This can be added to the `.bashrc` file so that it gets loaded to the environment on start up as follows:

```
echo 'export APPTAINER_TMPDIR="/lustre/projects/Research_Project-MRC147594/exome_sequencing/NF_temp/"' \
echo 'export APPTAINER_CACHEDIR=/lustre/projects/Research_Project-MRC147594/genome_sequencing/Genome_Pipeline_Containers/"' \
export 
>> ~/.bashrc
```

After this you will need to log out and back in again to the server or run `source ~/.bashrc` to restart your environment.

If you experience issues with this, these variables can be specified in later run_exomes.sh script.

## Running the pipeline 

### Using tmux (optional, but highly recommended)

Because the pipeline takes 1-3 days to run (depending on batch size) it is highly recommended to run it within a `tmux` session. This means that if you get disconnected from the server the pipeline will continue to run and you can return to the terminal session it is running in at any time.

To do this first install tmux with `micromamba install tmux` then start a session with the command `tmux`. The session can be left at any time without killing the job by using the key combination `ctrl` + `b`, followed by `d` and rejoined at any time by running `tmux attach`. To kill the the session simply run `exit`. For more information on running tmux run `man tmux`.

### Setting up the control file

The pipeline will look for a `control` file with in the batch directory with sample information.

This file must contain MALE and FEMALE lines containing the sexes of all samples in the batch, FAMILY lines for each family/singleton in the batch, and FASTQ lines for each fastq file pair in the batch. 

Here is an example:

```
MALE  EX_WES_865
FAMILY  EX2609030
FASTQ EX2609030 ILLUMINA  23N3GWLT3_L1  TwEx2_EX2609030_23N3GWLT3_S82_L001_R1_001.fastq.gz  TwEx2_EX2609030_23N3GWLT3_S82_L001_R2_001.fastq.gz
FASTQ EX2609030 ILLUMINA  23N3GWLT3_L2  TwEx2_EX2609030_23N3GWLT3_S82_L002_R1_001.fastq.gz  TwEx2_EX2609030_23N3GWLT3_S82_L002_R2_001.fastq.gz
FASTQ EX2609030 ILLUMINA  23N3GWLT3_L3  TwEx2_EX2609030_23N3GWLT3_S82_L003_R1_001.fastq.gz  TwEx2_EX2609030_23N3GWLT3_S82_L003_R2_001.fastq.gz
```
The python script make_control.py can generate this file for you. Input requires a .txt file with columns ID and Sex (coded as Male or Female).

The pipeline can also be run using bam files where control.tsv will look like:

```
MALE	TwEx2_EX2502154
FAMILY	TwEx2_EX2502154
BAM	TwEx2_EX2502154	ILLUMINA	223NL2LT1	TwEx2_EX2502154_223NL2LT1.bam
```

## Important pipeline file requirements
### File inputs
The bash script to run the exome pipeline requires specific file inputs and locations. If not correct, the pipeline will error.

The Exome pipeline requires that fastq files be present in a folder called 01_reads. Therefore, the --batchdir path in the run_exomes.sh should be structured as: /path/to/fastq/01_reads

The pipeline also requires a .bed file which defines the targetted WES regions, which is called in to use with --target_bed

--referenceFasta should point towards the .fasta file appropriate to your data. 

As the pipeline uses bwa-mem1 and bwa-mem2, you will a folder in the directory where --referenceFasta points to called bwa-mem1.
This needs to contain fasta index files: .amb, .ann, .bwt, .pac, .sa files.

In the --referenceFasta path, you will need the fasta index files for bwa-mem2 which are: .fasta, .0123, .amb, .ann, .bwt.2bit.64, .pac.

You will also need a samtools fasta index .fai and .dict for the GATK process in the --referenceFasta path.
The same is needed for the chr6 specific fasta file (--chr6_fasta).

--cnvref points towards the required .cnn file for CNVkit which can be generated using:

```
cnvkit.py reference \
  --fasta /path/to/fasta.fasta \
  --targets /path/to/bed.bed \
  --output ./output.cnn
```

--cytomegaloFasta is leftover from the Genome Pipeline but acts as a pipeline check.

### Variant screening
--genePanelTsv
--genePanelScript
These point to the files required to identify specific genes within your exome cohort. 
The monogenic_gene_panel.tsv and gene_panel_report.py files are designed to identify and filter specific genes for the Monogenic research group.
These can be edited to look for other genes with different filtering criteria or if not needed, set up a .tsv with NA in the fields. 


## Launching the pipeline

To begin load the conda environment with `micromamba activate exome_pipeline_env`.

Once this is done the pipeline can then be launched with the following sbatch command followed by a bash script (replacing with your actual paths to each directory):

```
sbatch run_exomes.sh
```
Where an example of run_exomes.sh is:

```
#!/bin/bash
#SBATCH -J X0129-Plate_26_pipeline
#SBATCH -p mrcq
#SBATCH -A Research_Project-MRC147594
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH -o /lustre/projects/Research_Project-MRC147594/exome_sequencing/X0129-Plate_26/pipeline_%j.log
#SBATCH -e /lustre/projects/Research_Project-MRC147594/exome_sequencing/X0129-Plate_26/pipeline_%j.err

cd /lustre/projects/Research_Project-MRC147594/exome_sequencing

export APPTAINER_TMPDIR=/lustre/projects/Research_Project-MRC147594/exome_sequencing/NF_temp
export APPTAINER_CACHEDIR=/lustre/projects/Research_Project-MRC147594/genome_sequencing/Genome_Pipeline_Containers

nextflow run /lustre/projects/Research_Project-MRC147594/exome_sequencing/Exome_Pipeline/Exome_Pipeline_test/mainv2.nf \
  -profile isca \
  --control /ennis/projects/sequencing_service/X0129-Plate_26/control.tsv \
  --batchDir /ennis/projects/sequencing_service/X0129-Plate_26/01_reads \
  --target_bed /lustre/projects/Research_Project-MRC147594/exome_sequencing/Exome_Pipeline/Hg38_exome_comp_spikein_v2.0.2.bed.gz \
  --referenceFasta /lustre/projects/Research_Project-MRC147594/resources/grch38-v2-pipeline/Homo_sapiens_assembly38.fasta \
  --cnvref /lustre/projects/Research_Project-MRC147594/exome_sequencing/hg38_exome_comp_spikein_v2.0.2.cnn \
  --chr6_fasta /lustre/projects/Research_Project-MRC147594/resources/grch38-v2-pipeline/GRCh38_chr6.fasta \
  --cytomegaloFasta /lustre/projects/Research_Project-MRC147594/resources/grch38-v2-pipeline/cytomegalovirus/cytomegalovirus.fasta \
  --genePanelTsv /lustre/projects/Research_Project-MRC147594/resources/grch38-v2-pipeline/gene_panels/monogenic_gene_panel.tsv \
  --genePanelScript /lustre/projects/Research_Project-MRC147594/exome_sequencing/Exome_Pipeline/Exome_Pipeline_test/gene_panel_report.py
```

If the pipeline fails for any reason it can also be resumed by running the above with the addition of the `-resume` flag.

As the pipeline runs it will write the outputs back to the batch directory within subdirectories with an `r04_` prefix.

Once the pipeline finishes successfully it will automatically run a cleanup step that clears the `work/` directory that contains all temporary files. 
However please always check that this runs correctly and manually wipe any remaining files if any are left behind.

