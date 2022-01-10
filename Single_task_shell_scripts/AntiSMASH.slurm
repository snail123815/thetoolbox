#!/bin/bash
#SBATCH   --job-name=antsm_
#SBATCH --output=antismash_.out
#SBATCH  --error=antismash_.err
#SBATCH --ntasks=8
#SBATCH --mem=60G
#SBATCH --partition=cpu-short
#SBATCH --time=00:59:00
#SBATCH --mail-user="c.du@biology.leidenuniv.nl"
#SBATCH --mail-type="ALL"

TARGET="GCF_011751035.2_ASM1175103v2_genomic.gbff"
CPU=8
TITLE="MBT27"
DESCRIPTION="refseq_pacbio_only"
OUTPUT_NAME=antismash_"$TITLE"_"$DESCRIPTION"_$(date "+%Y%m%d_%H%M%S")
GFF3_FILE="" # Leave empty if not used

CONDA_ENV="/data/projects/pi-vriesendorpb/tools/anaconda3/envs/antismash"
source /home/duc/.bashrc

SCRATCH=/scratchdata/${SLURM_JOB_USER}/${SLURM_JOB_ID} # fresh for every job
TEMPRES=$SCRATCH/$OUTPUT_NAME
CWD=$(pwd)
OUTPUTDIR=$CWD/$OUTPUT_NAME
TARGET=$(realpath $TARGET)

########################################################
# ANTISMASH command build
CMD="antismash --cpus $CPU"
CMD+=" --taxon bacteria" # optional (bacteria, fungi)
#CMD+=" --fullhmmer"           # Run a whole-genome HMMer analysis.
#CMD+=" --cassis"              # Motif based prediction of SM gene cluster regions. (eukaryotic)
CMD+=" --clusterhmmer"        # Run a cluster-limited HMMer analysis.
CMD+=" --tigrfam"             # Annotate clusters using TIGRFam profiles.
CMD+=" --smcog-trees"         # Generate phylogenetic trees of sec. met. cluster orthologous groups.
#CMD+=" --tta-threshold 0.65"  # Lowest GC content to annotate TTA codons at (default: 0.65).
CMD+=" --cb-general"          # Compare identified clusters against a database of antiSMASH-predicted clusters.
CMD+=" --cb-subclusters"      # Compare identified clusters against known subclusters responsible for synthesising precursors.
CMD+=" --cb-knownclusters"    # Compare identified clusters against known gene clusters from the MIBiG database.
CMD+=" --asf"                 # Run active site finder analysis.
CMD+=" --pfam2go"             # Run Pfam to Gene Ontology mapping module.
CMD+=" --rre"                 # Run RREFinder precision mode on all RiPP gene clusters. Needs fimo
CMD+=" --cc-mibig"            # Run a comparison against the MIBiG dataset
CMD+=" --genefinding-tool none" # {glimmerhmm,prodigal,prodigal-m,none,error}
                              # Specify algorithm used for gene finding: GlimmerHMM,
                              # Prodigal, Prodigal Metagenomic/Anonymous mode, or
                              # none. The 'error' option will raise an error if
                              # genefinding is attempted. The 'none' option will not
                              # run genefinding. (default: error).
if [ "$GFF3_FILE" != "" ]; then
    CMD+=" --genefinding-gff3 "$GFF3_FILE # Specify GFF3 file to extract features from.
fi

CMD+=" --output-dir "$TEMPRES
CMD+=" --html-title "$TITLE
if [ "$DESCRIPTION" != "" ]; then
    CMD+=" --html-description "$DESCRIPTION
fi
CMD+=" "$(basename $TARGET)
########################################################

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
cp -r "$TARGET" "$SCRATCH"
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
cd "$SCRATCH" || exit

echo "[$SHELL] Run: AntiSMASH"
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
