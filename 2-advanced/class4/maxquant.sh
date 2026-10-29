# https://maxquant.org/
# https://a3s.fi/proteomics/MaxQuant_tutorial.tar.gz
#MAXQUANT_HOME=/mnt/tools/MaxQuant_v2.8.1.0
MAXQUANT_HOME=/mnt/tools/MaxQuant_v.2.9.0

# Print help message to make sure software is working
dotnet $MAXQUANT_HOME/bin/MaxQuantCmd.dll --help

# Create and move to an output directory in your home directory
DATA_DIR=~/maxquant
mkdir $DATA_DIR
cd $DATA_DIR

## Examine the parameter file
nano mqpar.xml

# Avoid error by generating empty sdrf.tsv file
#mkdir -p combined/sdrf/
#touch combined/sdrf/sdrf.tsv

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# Do a dry run to show the planned steps
dotnet $MAXQUANT_HOME/bin/MaxQuantCmd.dll mqpar.xml -n

# Run the tool
dotnet $MAXQUANT_HOME/bin/MaxQuantCmd.dll mqpar.xml 2>&1 | tee maxquant.log

# Inspect the output
cd combined/txt
ls
nano proteinGroups.txt
