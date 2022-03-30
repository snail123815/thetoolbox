#!/bin/bash
#SBATCH --job-name=pgJJ1683
#SBATCH   --output=pgJJ1683.out
#SBATCH    --error=pgJJ1683.err
#SBATCH --ntasks=12
#SBATCH --mem=36G
#SBATCH --partition=cpu-short
#SBATCH --time=03:59:00
#SBATCH --mail-user="c.du@biology.leidenuniv.nl"
#SBATCH --mail-type="ALL"

NCPU="8"

TARGET="assembly.fasta"
GENUS_SPECIES="Acidithiobacillus"
STRAIN="JJ9999"

COMMENT="Annoted with PGAP on Alice"
BIOPROJECT="PRJNA999999"
BIOSAMPLE="SAMN99999999"
LOCUSTAG="TELOTA"

OUTPUT_NAME="pgap_"_$(date "+%Y%m%d_%H%M%S")
TAXCHECK_ONLY=1
TAXCHECK=1
IGNORE_ERRORS=1

TARGETFILE=$(basename $TARGET)
OUTPUTDIR=$(dirname $TARGET)


CONDA_ENV="/home/duc/apps/micromamba"


source /home/duc/.bashrc

PGAP_PATH="/data/projects/pi-vriesendorpb/tools/pgap"
SCRATCH=/scratchdata/${SLURM_JOB_USER}/${SLURM_JOB_ID} # fresh for every job
INPUT=$SCRATCH"/input.yaml"
SUBMOL=$SCRATCH"/submol.yaml"

TEMPRES=$SCRATCH/$OUTPUT_NAME/
CWD=$(pwd)

###################################################
# prepare yaml files

echo "fasta:
  class: File
  location: $TARGETFILE
submol:
  class: File
  location: submol.yaml" > $INPUT
echo "[$SHELL] input.yaml:"
echo "$(cat $INPUT)"

echo "organism: # REQUIRED
    genus_species: '$GENUS_SPECIES' # full name 'Streptomyces avermitilis', species name can leave empty like 'Streptomyces', mast be valid in [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)
    strain: '$STRAIN' # optional
comment: '$COMMENT' # optional, free text
bioproject: '$BIOPROJECT' # requried if submit, place holder
biosample: '$BIOSAMPLE' # requried if submit, place holder
locus_tag_prefix: '$LOCUSTAG' # required if submit, place holder" > $SUBMOL
echo "[$SHELL] submol.yaml:"
echo "$(cat $SUBMOL)"


##################################################
# Build command
CMD="python3 pgap/pgap.py -n"
CMD+=" -o "$TEMPRES
CMD+=" --container-name pgap-container"
CMD+=" --no-internet"
CMD+=" -D singularity"
if [ "$TAXCHECK_ONLY" -eq "1" ]; then; CMD+=" --taxcheck-only"
else; if [ "$TAXCHECK" -eq "1" ]; then; CMD+=" --taxcheck"; fi
if [ "$IGNORE_ERRORS" -eq "1" ]; then; CMD+=" --ignore-all-errors"; fi
CMD+=" input.yaml"
##################################################

# Adding [$SHELL] will give distinction between echo and script output
echo "[$SHELL] #### Starting Job"
# The SLURM_JOB_USER and SLURM_JOB_ID is automatically obtained from sbatch command
echo "[$SHELL] This is $SLURM_JOB_USER, job ID $SLURM_JOB_ID"
echo "[$SHELL] Node scratch: ""$SCRATCH"
# Enter the designed environment
echo "[$SHELL] Activate python3 and Singularity at $CONDA_ENV"
micromamba activate $CONDA_ENV
module load Singularity
echo "[$SHELL] $(python3 -V)"
echo "[$SHELL] $(module list Singularity 2>&1 | grep '1)')"
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
echo "[$SHELL] Copy data to scratch (aprox. 4 min)"
cp -r $PGAP_PATH $SCRATCH
cp $TARGET $INPUT $SUBMOL $SCRATCH 
echo "[$SHELL] SCRATCH dir: "
echo "$(ls -al $SCRATCH/)"
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
cd "$SCRATCH" || exit

echo "[$SHELL] Run: PGAP"
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
