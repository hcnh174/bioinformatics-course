# Lecture 1 Worksheet
# The worksheet can be downloaded from the link below; copying and pasting from there is recommended.
# https://home.hiroshima-u.ac.jp/nhayes/course-site/#/course/bioinformatics/2026

# Running the Pipeline

# 1 bash script
# Activate the “genome” conda environment
conda activate genome

# Make a copy of the pipeline directory
cp -r ~/data/bash_pipeline ~/

# Move to the pipeline directory (adjust to where you placed it)
cd ~/bash_pipeline

# Run all steps
bash run_all.sh

# Run a specific step only
bash run_all.sh 05_mutect2_somatic

#####################################################################

# 2 Snakemake
# Activate the “snakemake” conda environment
conda activate snakemake

# Make a copy of the pipeline directory
cp -r ~/data/snakemake_pipeline ~/

# Move to the pipeline directory
cd ~/snakemake_pipeline

# Dry run (check the flow / DAG without executing)
snakemake -n

# Run with 4 cores (auto-build conda envs, using mamba)
snakemake --cores 4 --use-conda --conda-frontend mamba

#####################################################################

# 3 Nextflow
# Activate the “nextflow” conda environment
conda activate nextflow

# Make a copy of the pipeline directory
cp -r ~/data/nextflow_pipeline ~/

# Move to the pipeline directory
cd ~/nextflow_pipeline

# Check the wiring (no tools needed)
nextflow run main.nf -stub-run

# Run, auto-building conda environments
nextflow run main.nf -profile conda
