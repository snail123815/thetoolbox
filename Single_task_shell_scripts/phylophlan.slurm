#!/bin/bash
#SBATCH --job-name=phyl
#SBATCH --output=phyl.out
#SBATCH  --error=phyl.err
#SBATCH --ntasks=20
#SBATCH --mem=150G
#SBATCH --partition=cpu-medium
#SBATCH --time=23:59:00
#SBATCH --mail-user="c.du@biology.leidenuniv.nl"
#SBATCH --mail-type="END"

CPU=20
baseDir=/data1/duc/Proj_tree
sourceDir=ready_genomes
outputDir="$baseDir"/"$sourceDir"-treeBuild

configFile="$baseDir"/phylophlan-dna-dna.cfg
maas="$baseDir"/phylophlan_subsititution_models.tsv
mkdir -p $outputDir
CONDA_ENV="/data1/projects/pi-vriesendorpb/condaEnvs/phylophlan"
source /home/duc/.bashrc

########################################################
# Command build

MKCFG="phylophlan_write_config_file "
MKCFG+=" -o "$configFile
MKCFG+=" -d n"
MKCFG+=" --db_dna makeblastdb"
MKCFG+=" --map_dna blastn"
MKCFG+=" --msa muscle"
MKCFG+=" --trim trimal"
MKCFG+=" --tree1 fasttree"
MKCFG+=" --tree2 raxml"
MKCFG+=" --overwrite"

CMD_CLEAN="phylophlan"
CMD_CLEAN+=" -c ""$baseDir"/"$sourceDir"
CMD_CLEAN+=" -d scoCoreGenes"
CMD_CLEAN+=" --maas "$maas
CMD_CLEAN+=" --diversity low"
CMD_CLEAN+=" --nproc ""$CPU"
CMD_CLEAN+=" --output_folder ""$outputDir"
CMD_CLEAN+=" -f ""$configFile"

CMD="phylophlan"
CMD+=" -i ""$baseDir"/"$sourceDir"
CMD+=" -d scoCoreGenes"
CMD+=" --maas "$maas
CMD+=" --diversity low"
CMD+=" --nproc ""$CPU"
CMD+=" --output_folder ""$outputDir"
CMD+=" -f ""$configFile"
#CMD+=" --configs_folder ""$configDir" # config folder needs to be at 
## $CONDA_ENV/lib/python3.1/site-packages/phylophlan/phylophlan_configs/
## phylophlan_write_config_file $CONDA_ENV/lib/python3.1/site-packages/phylophlan/phylophlan_configs/
## Problem is the `python3.1` in the path.

CMD+=" --verbose 2>&1 |"
CMD+=" tee ""$outputDir"/phylophlan_output.log

########################################################

# Adding [$SHELL] will give distinction between echo and script output
echo "[$SHELL] #### Starting Job"
# The SLURM_JOB_USER and SLURM_JOB_ID is automatically obtained from sbatch command
echo "[$SHELL] This is $SLURM_JOB_USER, job ID $SLURM_JOB_ID"
echo "[$SHELL] Node scratch: ""$SCRATCH"
# Enter the designed environment
echo "[$SHELL] Activate env at $CONDA_ENV"

micromamba activate $CONDA_ENV
echo "[$SHELL] "$(which phylophlan)

# Get start time
START_TIME=$(date +%s)
echo "[$SHELL] Started at: $(date)"
# Current working dir
echo "[$SHELL] CWD: ""$CWD"
# Change language settings in case perl complains
export LANGUAGE=en_US.UTF-8
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# Start the work
#echo "[$SHELL] Copy target to scratch."
#cp -r "$TARGET" "$SCRATCH"
#echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
#cd "$SCRATCH" || exit

echo "[$SHELL] Run make config file command:"
echo "[$SHELL] \$ $MKCFG"
eval "$MKCFG"

echo "[$SHELL] Run cleaning command:"
echo "[$SHELL] \$ $CMD_CLEAN"
eval "$CMD_CLEAN"

echo "[$SHELL] Run tree building command:"
echo "[$SHELL] \$ $CMD"
eval "$CMD"

#if [ -d "$TEMPRES" ]; then
#    echo "[$SHELL] Command execution finished."
#    echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
#    echo "[$SHELL] Save to $OUTPUTDIR"
#    mkdir -p "$OUTPUTDIR"
#    cp -a "$TEMPRES" "$OUTPUTDIR"
#    echo "[$SHELL] #### DONE"
#else
#    echo "[$SHELL] Output dir $TEMPRES not found, command execution failed."
#fi

echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
