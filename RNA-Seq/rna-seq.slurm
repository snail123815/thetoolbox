#!/bin/bash
#SBATCH --job-name=RSeq_
#SBATCH --output=RNASeq_.out
#SBATCH  --error=RNASeq_.err
#SBATCH --ntasks=8
#SBATCH --mem=10G
#SBATCH --partition=cpu-medium
#SBATCH --time=06:59:00
#SBATCH --mail-user="c.du@biology.leidenuniv.nl"
#SBATCH --mail-type="ALL"

NCPU="8"
RAW="/data/duc/RNA-SeqProj/rawreads"
OUT_ALIGN="/data/duc/RNA-SeqProj/aligned"
OUT_COUNT="/data/duc/RNA-SeqProj/countTable"
GENOME="/data/duc/RNA-SeqProj/genome.gbk" # genbank file with annotation required
TARGET_FEATURE="locus_tag"
GROUP_FEATURE=$TARGET_FEATURE

# Pair-end settings:
PAIREND=1
RAWEXT=""
PESUFFIX=""
SNAMES=""

# featureCounts style
FRACTION_COUNTING=0
LOOSE_PE_COUNTING=0

SCRATCH=/scratchdata/${SLURM_JOB_USER}/${SLURM_JOB_ID} # fresh for every job
SCRIPTS_DIR="/home/duc/toolbox-du/RNA-Seq"
CONDA_ENV="/data/projects/pi-vriesendorpb/tools/anaconda3/envs/RNASeqBact"
CWD=$(pwd)
TEMP_RAW=$SCRATCH/rawreads
TEMP_ALIGN=$SCRATCH/aligned
TEMP_COUNT=$SCRATCH/countTable
#######################
# Build commands align
CMD_ALIGN="python3 $SCRIPTS_DIR/Align_RNASeq_bowtie2.py"
CMD_ALIGN+=" --ncpu "$NCPU              # number of cpu to use
CMD_ALIGN+=" --raw "$TEMP_RAW                # path(s) to raw data
CMD_ALIGN+=" --out "$TEMP_ALIGN          # path to output alignment
CMD_ALIGN+=" --genome "$GENOME          # path to genome file(s), also supports indexed genome (*.bt2)
if [ $PAIREND = 1 ]; then
    CMD_ALIGN+=" --isPe"
fi
if [ "$RAWEXT" != "" ]; then
    CMD_ALIGN+=" --rawExt "$RAWEXT      # extension ".fastq.gz" or ".fq.gz"
fi
if [ "$PESUFFIX" != "" ]; then
    CMD_ALIGN+=" --pesuffix "$PESUFFIX  # suffix to pairend file name. eg. a_1_.fq.gz and a_2_.fq.gz,
                                        # set this to --pesuffix _1_ _2_
fi
if [ "$SNAMES" != "" ]; then
    CMD_ALIGN+=" --sampleNames "$SNAMES # sample names if file names contain random generated string which
                                        # will unpair the samples after removing pairend suffix
fi

########################
# Build commands count
CMD_COUNT="python3 $SCRIPTS_DIR/CountReads_RNASeq_featureCounts.py"
CMD_COUNT+=" --input "$TEMP_ALIGN    # alignment folder (bam or sam files)
CMD_COUNT+=" --output "$TEMP_COUNT   # gene counts output folder
CMD_COUNT+=" --gbk "$GENOME          # genbank file with annotation
CMD_COUNT+=" --ncpu "$NCPU           # number of cpu to use
if [ $PAIREND = 1 ]; then
    CMD_COUNT+=" --isPe"             # set if pairend
fi
CMD_COUNT+=" -t "$TARGET_FEATURE     # Specify feature type(s) in a GTF annotation. If multiple
                                     # types are provided, they should be separated by ',' with
                                     # no space in between. 'exon' by default. Only rows in the
                                     # annotation with a matched feature will be extracted and
                                     # used for read mapping.
CMD_COUNT+=" -g "$GROUP_FEATURE      # Specify attribute type in GTF annotation. 'gene_id' by
                                     # default. Meta-features used for read counting will be
                                     # extracted from annotation using the provided value.
                                     # This attribute tyep is usually the gene identifier.
                                     # This argument is useful for the meta-feature level summarization.
if [ $FRACTION_COUNTING = 1 ]; then
    CMD_COUNT+=" --fractionCounting" # will add -M --fraction -O if True
fi
if [ $LOOSE_PE_COUNTING = 1 ]; then
    CMD_COUNT+=" --peLoose"          # will use loose configuration in pairend counting, set it if you do not
                                     # want the features: -P, Check validity of paired-end distance when
                                     # counting read pairs. Use -d and -D to set thresholds; -B, only count
                                     # read pairs that have both ends aligned. (must set togeter with -P)
fi
########################

echo "[$SHELL] #### Starting Job"
echo "[$SHELL] This is $SLURM_JOB_USER, job ID $SLURM_JOB_ID"
echo "[$SHELL] Node scratch: ""$SCRATCH"
echo "[$SHELL] Activate environment at $CONDA_ENV"
source /home/duc/.bashrc
micromamba activate $CONDA_ENV
echo "[$SHELL] Current environment python: $(which python)"
START_TIME=$(date +%s)
echo "[$SHELL] Started at: $(date)"
echo "[$SHELL] CWD: ""$CWD"
# Change language settings in case perl complains
export LANGUAGE=en_US.UTF-8
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# Start the work
echo "[$SHELL] Copy data to scratch."
mkdir -p "$TEMP_RAW"
cp -r "$RAW"/* "$TEMP_RAW"
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
cd "$SCRATCH" || exit

# Alignment
echo "[$SHELL] Run: alignment"
echo "[$SHELL] Command: \$ $CMD_ALIGN"
eval "$CMD_ALIGN"

if [ -d "$TEMP_ALIGN" ]; then
    echo "[$SHELL] Command execution finished."
    echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
    echo "[$SHELL] Save to $OUT_ALIGN"
    mkdir -p "$OUT_ALIGN"
    cp -r "$TEMP_ALIGN"/* "$OUT_ALIGN"
    echo "[$SHELL] #### DONE ALIGNMNET"
else
    echo "[$SHELL] Output dir $TEMP_ALIGN not found, command execution failed."
    echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
    exit
fi

# Count
echo "[$SHELL] Run: count"
echo "[$SHELL] Command: \$ $CMD_COUNT"
eval "$CMD_COUNT"

if [ -d "$TEMP_COUNT" ]; then
    echo "[$SHELL] Command execution finished."
    echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
    echo "[$SHELL] Save to $OUT_COUNT"
    mkdir -p "$OUT_COUNT"
    cp -r "$TEMP_COUNT"/* "$OUT_COUNT"
    echo "[$SHELL] #### DONE COUNT"
else
    echo "[$SHELL] Output dir $TEMP_COUNT not found, command execution failed."
    echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"
    exit
fi

echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - START_TIME))" +%T) (HH:MM:SS)"