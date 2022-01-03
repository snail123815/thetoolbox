#!/bin/bash
#SBATCH --job-name=prok_MBT
#SBATCH --output=prokka_MBT.out
#SBATCH  --error=prokka_MBT.err
#SBATCH --ntasks=8
#SBATCH --mem=10G
#SBATCH --partition=cpu-short
#SBATCH --time=01:00:00
#SBATCH --mail-user="c.du@biology.leidenuniv.nl"
#SBATCH --mail-type="ALL"

NCPU="8"

TARGET="MBT.fa"
GENUS="actinobacteria"
SPECIES="sp."
STRAIN="MBT"
LOCUSTAG="LOCUS_TAG_MBT"
PREFIX="prokka_MBT" # Filename output prefix [auto] (default '')
OUTPUT_NAME=$PREFIX"_$STRAIN""_$LOCUSTAG"_$(date "+%Y%m%d_%H%M%S")

GRAM="POS"
CENTRE="MBT"
GCODE="11"

CONDA_ENV="/data/projects/pi-vriesendorpb/tools/anaconda3/envs/prokka_with_signalp"
source /home/duc/.bashrc

SCRATCH=/scratchdata/${SLURM_JOB_USER}/${SLURM_JOB_ID} # fresh for every job
TEMPRES=$SCRATCH/$OUTPUT_NAME
CWD=$(pwd)
OUTPUTDIR=$CWD/$OUTPUT_NAME

##################################################
# Build command
CMD="prokka --compliant --addgenes --mincontiglen 200 --rfam"
CMD+=" --centre "$CENTRE
CMD+=" --genus "$GENUS
CMD+=" --species "$SPECIES
CMD+=" --strain "$STRAIN
CMD+=" --locustag "$LOCUSTAG
CMD+=" --gram "$GRAM
CMD+=" --prefix "$PREFIX
CMD+=" --gcode "$GCODE
CMD+=" --cpu "$NCPU
CMD+=" "$TARGET
##################################################

# Adding [$SHELL] will give distinction between echo and script output
echo "[$SHELL] #### Starting Job"
# The SLURM_JOB_USER and SLURM_JOB_ID is automatically obtained from sbatch command
echo "[$SHELL] This is $SLURM_JOB_USER, job ID $SLURM_JOB_ID"
echo "[$SHELL] Node scratch: ""$SCRATCH"
# Enter the designed environment
echo "[$SHELL] Activate AntiSMASH at $CONDA_ENV"
micromamba activate $CONDA_ENV
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
echo "[$SHELL] Copy target to scratch."
cp -r "$CWD/$TARGET" "$RUNDIR"
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
cd "$RUNDIR" || exit

echo "[$SHELL] Run: Prokka"
echo "[$SHELL] Command: \$ $CMD"
eval "$CMD"

if [ -d "$TEMPRES" ]; then
    echo "[$SHELL] Command execution finished."
    echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
    echo "[$SHELL] Save to $OUTPUTDIR"
    mkdir -p "$OUTPUTDIR"
    cp -a "$TEMPRES" "$OUTPUTDIR"
    echo "[$SHELL] #### DONE"
else
    echo "[$SHELL] Output dir $TEMPRES not found, command execution failed."
fi

echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
