import argparse
import os
from multiprocessing import Pool
from glob import glob
from tqdm import tqdm
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
import re
from typing import Literal, TypedDict
from pyBioinfo_modules.wrappers.antismash \
    import findClusterNumberStr, clusterGbkGlobTxt
from pyBioinfo_modules.wrappers.mash \
    import calculate_medoid, mashSketchFiles, mashDistance
from pyBioinfo_modules.wrappers.bigscape \
    import runBigscape
from tempfile import TemporaryDirectory
import pickle


class ClusterInfo(TypedDict):
    gbkFile: Path
    gcProducts: str
    organism: str
    fromSequence: str
    coreRelativeLocs: list[FeatureLocation]
    joinedProteinFastaFile: Path
    fastaId: str


def parseClusterGbk(
    infile: Path,
    proteinFastaDir: Path,
    id: int | None = None,
    nflank: int = 0,
) -> ClusterInfo:
    """Parses the genbank files for DNA, protein, cluster, organism
    [Rewrote of parsegbkcluster() from BiG-MAP]
    parameters
    ----------
    infile
        Path of .gbk file
    proteinFastaDir
        Path of protein fastas to store joined proteins for each cluster
    nflank
        Number of CDS to flank the core regions
    returns
    ----------
    DNA = DNA sequence
    proteins = protein sequence
    clustername = name of the cluster
    organism = name of the organism
    #cluster_enzymes = {loc:gene_kind}
    """
    proteins = []
    gcProducts = []
    cdsIndexs = []
    coreIndexs = []
    coreRelativeLocs = []
    records = list(SeqIO.parse(infile, "genbank"))
    assert len(records) == 1
    record = records[0]
    for i, feature in enumerate(record.features):
        # Parsing the regionname, part of antismash output
        if "region" in feature.type:
            if "product" in feature.qualifiers:
                for product in feature.qualifiers["product"]:
                    gcProducts.append(product)
        # Parsing the protein sequence
        if feature.type == "CDS":
            cdsIndexs.append(i)
            # remembering every CDS gene index
            try:
                proteins.append(feature.qualifiers['translation'][0])
            except KeyError:
                pass
            # Parsing the relative core locations
            try:
                if feature.qualifiers['gene_kind'][0] == "biosynthetic":
                    coreIndexs.append(i)
            except KeyError:
                pass

    if len(coreIndexs) > 0:
        if cdsIndexs.index(min(coreIndexs)) - nflank < 0 or \
                cdsIndexs.index(max(coreIndexs)) + nflank + 1 > len(cdsIndexs):
            print(
                "!!!flank_genes (-f) is higher than the number of flanking " +
                f"genes in the cluster of file: {infile}, " +
                "using whole gene cluster instead!!!"
            )
        if nflank == 0:
            coreRelativeLocs = [
                record.features[i].location for i in coreIndexs]
        else:
            if cdsIndexs.index(min(coreIndexs)) - nflank >= 0:
                core_region = cdsIndexs[
                    cdsIndexs.index(min(coreIndexs)) -
                    nflank:cdsIndexs.index(max(coreIndexs)) + nflank + 1
                ]
            else:
                core_region = cdsIndexs[0:cdsIndexs.index(
                    max(coreIndexs)) + nflank + 1]
            coreRelativeLocs = [
                record.features[i].location for i in core_region]

    organism = record.description
    # Parsing organism name
    # The attr. description is from "DEFINATION" line for gbk record
    # biopython/Bio/GenBank/__init__.py line 764
    # in class _FeatureConsumer(_BaseGenBankConsumer)
    if " " in organism:
        organism = "_".join(organism.split(",")[0].split())
    organism = re.sub(f'strain', '', organism)
    organism = re.sub(f'[()"]+', '', organism)
    organism = re.sub(f'[\\/]+', '-', organism)
    organism = re.sub(f'[.]+', '.', organism)
    organism = re.sub(f'[_]+', '_', organism)

    # write to protein fasta file
    joinProteins = Seq(''.join(proteins))
    regionAndNumber = findClusterNumberStr(infile)
    fromSequence = infile.name.split(regionAndNumber)[0][:-1]
    if id is None:
        proteinFastaFileName = \
            f"GC_PROT-{fromSequence}-{regionAndNumber}" + ".fasta"
    else:
        proteinFastaFileName = f"P{id}"
    fastaId = f"{organism}|{fromSequence}|" + \
        f"GC_PROT--{regionAndNumber}" + f"--Entryname={':'.join(gcProducts)}"
    proteinFastaFile = proteinFastaDir / proteinFastaFileName
    SeqIO.write(SeqRecord(joinProteins, id=fastaId, description=""),
                proteinFastaFile, 'fasta')
    assert proteinFastaFile.is_file()

    return ClusterInfo(
        gbkFile=infile,
        gcProducts=":".join(gcProducts),
        organism=organism,
        fromSequence=fromSequence,
        coreRelativeLocs=coreRelativeLocs,
        joinedProteinFastaFile=proteinFastaFile,
        fastaId=fastaId
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('p', type=Path,
                        help='Input dir')
    parser.add_argument('--cpus', type=int,
                        help='Processes number', default=4)
    parser.add_argument('--tmp', type=Path,
                        help='temporary files folder', default=None)
    parser.add_argument('--mashClusterRatio', type=float, default=0.8,
                        help='Ratio of match/total for two BGC to be clustered')
    parser.add_argument('--verboseBigscape', type=bool,
                        help='add --verbose to bigscape run', default=False)
    parser.add_argument('--outputPath', type=Path,
                        help='Bigscape results', required=True)

    args = parser.parse_args()
    inputPath: Path = args.p

    # 1. get gbk files from input path
    # 2. get cluster info for each gbk files
    # 2.1 print protein sequences to fasta
    # 3. make sketch and calculate distance
    # 4. get medoid from distance
    # 5. move gbk files together for BiG-SCAPE
    # 6. run BiG-SCAPE

    # Prepare dirs
    if args.tmp is None:
        fastaTempD = TemporaryDirectory()
        mashTempD = TemporaryDirectory()
        representativeGbksTempD = TemporaryDirectory()
        proteinFastaDir = Path(fastaTempD.name)
        proteinMashDir = Path(mashTempD.name)
        representativeGbksDir = Path(representativeGbksTempD.name)
    else:
        if not args.tmp.exists():
            args.tmp.mkdir(parents=True)
        proteinFastaDir = (args.tmp / 'proteinFastas').resolve()
        proteinMashDir = args.tmp / 'mashSketchAndDist'
        representativeGbksDir = args.tmp / 'representativeGbks'
        proteinFastaDir.mkdir(exist_ok=True)
        proteinMashDir.mkdir(exist_ok=True)
        representativeGbksDir.mkdir(exist_ok=True)

    sketchFile = (proteinMashDir / 'GC_PROT.msh').resolve()
    distanceTableFile = proteinMashDir / 'mash_output_GC.tab'
    mashTableFinishedFlagFile = proteinMashDir / 'distFinished'
    gbksListFile = proteinMashDir / 'gbks_list.pickle'
    clusterInfoDictFile = proteinMashDir / 'clusters_info_dict.pickle'
    familyMedoidDictFile = proteinMashDir / 'family_medoid_dict.pickle'
    familyDistMatrixDictFile = proteinMashDir / 'family_dist_matrix_dict.pickle'
    familyMedoidFinishedFlagFile = proteinMashDir / 'medoidFinished'

    allRepresentativeGbksMovedFlagFile = proteinMashDir / 'ALLIN_rep_gbks'

    args.outputPath.mkdir(exist_ok=True)
    bigscapeOutput = args.outputPath / 'bigscape'

    try:
        if gbksListFile.exists():
            with gbksListFile.open('rb') as fh:
                gbks = pickle.load(fh)
        else:
            gbks = list(inputPath.glob(clusterGbkGlobTxt))
            for d in tqdm([d for d in inputPath.iterdir() if d.is_dir()],
                          desc=(f'Gathering gbk files from {inputPath.name}')):
                gbks.extend(Path(gbk) for gbk in
                            glob(str(d / ('**/' + clusterGbkGlobTxt)),
                            recursive=True))
            if args.tmp is not None:
                with gbksListFile.open('wb') as fh:
                    pickle.dump(gbks, fh)

        if clusterInfoDictFile.exists():
            with clusterInfoDictFile.open('rb') as fh:
                clusterInfoDict = pickle.load(fh)
        else:
            with Pool(args.cpus) as readGbkPool:
                results = [
                    readGbkPool.apply_async(
                        parseClusterGbk, (gbk, proteinFastaDir, i))
                    for i, gbk in enumerate(gbks)
                ]
                # clusterInfos: list[ClusterInfo] = [
                #     r.get() for r in tqdm(results, desc='Parsing cluster info')
                # ]
                clusterInfoDict = {}
                for r in tqdm(results, desc='Parsing cluster info'):
                    ci = r.get()
                    clusterInfoDict[ci['joinedProteinFastaFile'].name] = ci
            if args.tmp is not None:
                with clusterInfoDictFile.open('wb') as fh:
                    pickle.dump(clusterInfoDict, fh)

        if not mashTableFinishedFlagFile.exists():
            cwd = Path().resolve()
            os.chdir(proteinFastaDir)
            print('Making sketch...')
            mashSketchFiles(
                proteinFastaDir.glob('P*'),
                sketchFile.with_suffix(''),
                kmer=16,
                sketch=5000,
                nthreads=args.cpus,
                molecule='protein'
            )
            os.chdir(cwd)
            print("Calculationg distance...")
            mashDistance(
                sketchFile,
                distanceTableFile,
                nthreads=args.cpus,
            )
            mashTableFinishedFlagFile.touch()
        assert distanceTableFile.exists()

        if familyMedoidFinishedFlagFile.exists():
            with familyMedoidDictFile.open('rb') as fh:
                dict_medoids = pickle.load(fh)
            with familyDistMatrixDictFile.open('rb') as fh:
                family_distance_matrice = pickle.load(fh)
        else:
            print('Gather families and calculating medoid...')
            dict_medoids, family_distance_matrice = calculate_medoid(
                distanceTableFile,
                args.mashClusterRatio
            )
            if args.tmp is not None:
                with familyMedoidDictFile.open('wb') as fh:
                    pickle.dump(dict_medoids, fh)
                with familyDistMatrixDictFile.open('wb') as fh:
                    pickle.dump(family_distance_matrice, fh)
                familyMedoidFinishedFlagFile.touch()

        if not allRepresentativeGbksMovedFlagFile.exists():
            for name in dict_medoids.keys():
                repGbkFile = clusterInfoDict[name]['gbkFile']
                moveToFile = representativeGbksDir / repGbkFile.name
                moveToFile.symlink_to(repGbkFile.resolve())
            allRepresentativeGbksMovedFlagFile.touch()
        assert len(dict_medoids) == len([
            f for f in representativeGbksDir.iterdir() if f.is_file()
        ])

        runBigscape(
            representativeGbksDir, bigscapeOutput,
            cpus=args.cpus,
            cutoffs=[0.2, ],
            silent=False,
            verbose=args.verboseBigscape
        )

        print('FINISH')  # break point
    finally:
        if args.tmp is None:
            fastaTempD.cleanup()
            mashTempD.cleanup()
            representativeGbksTempD.cleanup()


if __name__ == "__main__":
    main()
