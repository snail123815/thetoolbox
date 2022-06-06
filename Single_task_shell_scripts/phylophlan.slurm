#!/bin/bash
#SBATCH   --job-name=phyl
#SBATCH --output=phyl.out
#SBATCH  --error=phyl.err
#SBATCH --ntasks=15
#SBATCH --mem=60G
#SBATCH --partition=cpu-short
#SBATCH --time=03:59:00
#SBATCH --mail-user="c.du@biology.leidenuniv.nl"
#SBATCH --mail-type="ALL"

CPU=15
baseDir=/data1/projects/pi-vriesendorpb/duc/Proj_tree
sourceDir=ncbi-ftp-download-Streptomyces-15
configFile="$baseDir"/ncbi-ftp-download-Streptomyces-phylophanDB.cfg
configDir="$baseDir"/phylophlan_configs
outputDir="$baseDir"/"$sourceDir"-tree-on-phylophlanDB
rm -rf $outputDir
mkdir -p $outputDir
CONDA_ENV="/data1/projects/pi-vriesendorpb/condaEnvs/phylophlan"
source /home/duc/.bashrc

########################################################
# Command build

CMD="phylophlan"
CMD+=" -i ""$baseDir"/"$sourceDir"
CMD+=" -d phylophlan"
CMD+=" --diversity low"
CMD+=" --nproc ""$CPU"
CMD+=" --output_folder ""$outputDir"
CMD+=" -f ""$configFile"
CMD+=" --configs_folder ""$configDir"
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

echo "[$SHELL] Run command:"
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
#
#echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
