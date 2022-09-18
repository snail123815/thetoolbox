import os
import argparse
import subprocess
import shutil
import logging
from tempfile import TemporaryDirectory
from multiprocessing import Pool


acceptedAntismashSwitches = (
    'bacteria', 'fungi',
    'fullhammer', 'cassis', 'clusterhmmer', 'tigrfam',
    'smcog-trees',
    'cb-general', 'cb-subclusters', 'cb-knownclusters',
    'asf', 'pfam2go', 'rre', 'cc-mibig',
)

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'p',
        help='Input folder containing gbk files (can be gz or xz)'
    )
    parser.add_argument(
        '--processes',
        type=int, default=8,
        help="Number of antismashes run in parallel"
    )
    parser.add_argument(
        '--threads', type=int, default=4,
        help="Number of threads for each antismash run"
    )
    parser.add_argument(
        '--switches', help=f'Comma seperated switches, allowed are {acceptedAntismashSwitches}',
        default='bacteria,cb-general,cb-subclusters,cb-knownclusters'
    )
    parser.add_argument(
        '--dry', action='store_true'
    )
    parser.add_argument(
        '--overwrite', action='store_true'
    )
    parser.add_argument(
        '--forceGeneFinding', action='store_true'
    )
    parser.add_argument(
        '--noProgress', action='store_true'
    )
    return parser.parse_args()


def cleanFiles(files, targetExts):
    newFiles = []
    for f in files:
        fn, ext = os.path.splitext(f)
        if ext in ['.gz', '.xz']:
            fn, formatExt = os.path.splitext(fn)
        else:
            formatExt = ext
        if formatExt in targetExts:
            newFiles.append(f)
    return newFiles

def runAntismash(
        root, file, outputRoot, threads=4,
        overwrite=False,
        switches=['bacteria', 'cb-general', 'cb-subclusters', 'cb-knownclusters'],
        forceGeneFinding=False,
        allowLongHeaders=True,
        dryrun=False,
    ):

    inputFile = os.path.join(root, file)
    for switch in switches:
        assert switch in acceptedAntismashSwitches, f'"{switch} not in {acceptedAntismashSwitches}'

    fn, formatExt = os.path.splitext(file)
    if formatExt in ['.gz', '.xz']:
        prog = ('gzip' if formatExt == '.gz' else 'xz')
        tempdir = TemporaryDirectory()
        unzippedFile = os.path.join(tempdir.name, fn)
        assert os.path.isfile(inputFile)
        if not dryrun:
            unzip = subprocess.run(f'{prog} -dc {inputFile} > {unzippedFile}',
                    shell=True)
            if unzip.returncode != 0:
                logging.error(f'unzip "{inputFile}" failed')
                return 'unzip failed', file
        inputFile = unzippedFile
        fn, formatExt = os.path.splitext(fn)
    outputDir = os.path.join(outputRoot, fn)
    
    # check if it has already been done
    if overwrite or not os.path.isfile(os.path.join(outputDir, fn+".zip")):
        if not dryrun and os.path.isdir(outputDir):
            shutil.rmtree(outputDir)
    else:
        return 'previous result exists', file
    if not dryrun:
        os.makedirs(outputDir)
    
        
    cmd = f"antismash --cpus {threads}"
    if 'bacteria' in switches and 'fungi' in switches:
        logging.exception('Only one of "bacteria" and "fungi" allowed')
        return 'argument error', file
    elif 'bacteria' in switches:
        cmd += " --taxon bacteria"
    elif 'fungi' in switches:
        cmd += " --taxon fungi"
    else:
        logging.exception('One of "bacteria" and "fungi" needs to be set with "switches"')
        return 'argument error', file
    if 'fullhammer' in switches:
        cmd += " --fullhmmer"           # Run a whole-genome HMMer analysis.
    if 'cassis' in switches:
        cmd += " --cassis"              # Motif based prediction of SM gene cluster regions. (eukaryotic)
    if 'clusterhmmer' in switches:
        cmd += " --clusterhmmer"        # Run a cluster-limited HMMer analysis.
    if 'tigrfam' in switches:
        cmd += " --tigrfam"             # Annotate clusters using TIGRFam profiles.
    if 'smcog-trees' in switches:
        cmd += " --smcog-trees"         # Generate phylogenetic trees of sec. met. cluster orthologous groups.
    if 'cb-general' in switches:
        cmd += " --cb-general"          # Compare identified clusters against a database of antiSMASH-predicted clusters.
    if 'cb-subclusters' in switches:
        cmd += " --cb-subclusters"      # Compare identified clusters against known subclusters responsible for synthesising precursors.
    if 'cb-knownclusters' in switches:
        cmd += " --cb-knownclusters"    # Compare identified clusters against known gene clusters from the MIBiG database.
    if 'asf' in switches:
        cmd += " --asf"                 # Run active site finder analysis.
    if 'pfam2go' in switches:
        cmd += " --pfam2go"             # Run Pfam to Gene Ontology mapping module.
    if 'rre' in switches:
        cmd += " --rre"                 # Run RREFinder precision mode on all RiPP gene clusters. Needs fimo
    if 'cc-mibig' in switches:
        cmd += " --cc-mibig"            # Run a comparison against the MIBiG dataset
    if formatExt in ['.gbk','.gb','.gbff'] and not forceGeneFinding:
        cmd += " --genefinding-tool none"
    else:
        cmd += " --genefinding-tool prodigal" # {glimmerhmm,prodigal,prodigal-m,none,error}
                                  # Specify algorithm used for gene finding: GlimmerHMM,
                                  # Prodigal, Prodigal Metagenomic/Anonymous mode, or
                                  # none. The 'error' option will raise an error if
                                  # genefinding is attempted. The 'none' option will not
                                  # run genefinding. (default: error).

    cmd += f" --output-dir {outputDir}"
    cmd += f" --html-title {fn}"
    if allowLongHeaders:
        cmd += " --allow-long-headers"
    cmd += f" {inputFile}"
    
    if not dryrun:
        antismash = subprocess.run(cmd, shell=True)
        if antismash.returncode != 0:
            logging.error('\nNot succeed on command:')
            logging.error('\t'+cmd)
            return 'Antismash error', file
        return 'success', file
    else:
        logging.info(cmd)
        return 'dry', file


if __name__ == "__main__":
    args = parseArgs()
    pathIn = os.path.realpath(args.p.strip())
    r, p = os.path.split(pathIn)
    pathOut = os.path.join(r, p+'-antismash')
    os.makedirs(pathOut, exist_ok=True)
    logging.basicConfig(filename=os.path.join(pathOut, 'antismash.log'), filemode='a', level='INFO')

    files = cleanFiles(os.listdir(pathIn), ['.gbk', '.gbff'])

    runnerPool = Pool(args.processes)
    results = []
    for f in files:
        results.append(runnerPool.apply_async(runAntismash, (pathIn, f, pathOut),
            kwds={
                'switches': args.switches.split(','),
                'threads': args.threads,
                'dryrun': args.dry,
                'forceGeneFinding': args.forceGeneFinding,
                'overwrite': args.overwrite,
            }
        ))
    runnerPool.close()
    returns = []
    if args.noProgress:
        for res in results:
            returns.append(res.get())
    else:
        from tqdm import tqdm
        for res in tqdm(results):
            returns.append(res.get())
    logging.info('#'*100)
    returns.sort(key=lambda x: x[0])
    for ret in returns:
        logging.info(ret)