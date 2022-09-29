import argparse
import shutil
from multiprocessing import Pool
from glob import glob
from tqdm import tqdm
from pathlib import Path
from turtle import distance
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


class ClusterInfo(TypedDict):
    dnaSeq: Seq
    joinProteins: Seq
    gcProducts: str
    organism: str
    coreRelativeLocs: list[FeatureLocation]
    gbkFilePath: Path


def parseClusterGbk(infile: Path, nflank: int = 0) -> ClusterInfo:
    """Parses the genbank files for DNA, protein, cluster, organism
    [Rewrote of parsegbkcluster() from BiG-MAP]
    parameters
    ----------
    infile
        string, name of the input .gbk file
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
    organism = re.sub(f'[._]+', '_', organism)

    return ClusterInfo(
        dnaSeq=record.seq,
        joinProteins=Seq(''.join(proteins)),
        gcProducts=":".join(gcProducts),
        organism=organism,
        coreRelativeLocs=coreRelativeLocs,
        gbkFilePath=infile
    )


def writeFasta(
    sequence: Seq,
    seqstype: Literal['GC_PROT', 'GC_DNA', 'HG_PROT', 'HG_DNA'],
    gcProducts: str,
    organism: str,
    infile: Path,
    outdir: Path
) -> tuple[Path, Path, str, str]:
    """Writes the fasta file for each sequence
    sequences(dnaSeq or joinProteins), gcProducts, organism = output of parseClusterGbk()
    infile = 000001.1.region001.gbk (no consecutive dot ... in file name)
    Kitasatospora_acidiphila_MMS16-CNU292_NZ_VIGB01000001.1.region001.gbk

    returns
    ----------
    outfile = the name of the written fasta file
    fromSequence = ID for the input organism genome
    fastaId = header of the cluster dna sequence
    """
    regionAndNumber = findClusterNumberStr(infile)
    fromSequence = infile.name.split(regionAndNumber)[0][:-1]
    fileName = (seqstype + gcProducts if 'HG' in seqstype else seqstype)
    fileName += f"-{fromSequence}-{regionAndNumber}"
    fastaId = f"{organism}|{fromSequence}|{seqstype}--Region{regionAndNumber}" \
        + f"--Entryname={gcProducts}"
    outfile = outdir / (fileName + '.fasta')
    SeqIO.write(SeqRecord(sequence), outfile, 'fasta')
    assert outfile.is_file()
    return infile, outfile, fromSequence, fastaId


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('p', type=Path, help='Input dir')
    parser.add_argument('--cpus', type=int, help='Processes number', default=4)
    parser.add_argument('--outputPath', type=Path, help='Bigscape results',
                        required=True)

    args = parser.parse_args()
    inputPath: Path = args.p

    # 1. get gbk files from input path
    # 2. get cluster info for each gbk files
    # 2.1 print protein sequences to fasta
    # 3. make sketch and calculate distance
    # 4. get medoid from distance
    # 5. move gbk files together for BiG-SCAPE
    # 6. run BiG-SCAPE

    gbks = list(inputPath.glob(clusterGbkGlobTxt))
    for d in tqdm([d for d in inputPath.iterdir() if d.is_dir()],
                  desc=(f'Gathering gbk files from {inputPath.name}')):
        gbks.extend(Path(gbk) for gbk in
                    glob(str(d / ('**/' + clusterGbkGlobTxt)), recursive=True))

    with Pool(args.cpus) as readGbkPool:
        results = [
            readGbkPool.apply_async(parseClusterGbk, (gbk,))
            for gbk in gbks
        ]
        clusterInfos: list[ClusterInfo] = [
            r.get() for r in tqdm(results, desc='Reading gbk files')
        ]

    fastaDir = TemporaryDirectory()
    mashDir = TemporaryDirectory()
    gbkFamiliesDir = TemporaryDirectory()
    try:
        fastaDirPath = Path(fastaDir.name)
        mashDirPath = Path(mashDir.name)
        gbkFamiliesDirPath = Path(gbkFamiliesDir.name)
        with Pool(args.cpus) as writeFastaPool:
            results = [
                writeFastaPool.apply_async(
                    writeFasta,
                    (
                        ci['joinProteins'],
                        'GC_PROT',
                        ci['gcProducts'],
                        ci['organism'],
                        ci['gbkFilePath'],
                        fastaDirPath
                    )
                ) for ci in clusterInfos
            ]
            fastaReturns = [r.get() for r in tqdm(
                results,
                desc='Writing protein fasta files...'
            )]
        fastaToGbk = {str(ret[1]): ret[0] for ret in fastaReturns}
        sketchFile = mashDirPath / 'GC_PROT.msh'
        print('Making sketch...')
        mashSketchFiles(
            fastaDirPath.glob('GC_PROT*.fasta'),
            sketchFile.with_suffix(''),
            kmer=16,
            sketch=5000,
            nthreads=args.cpus,
            molecule='protein'
        )
        distanceTable = mashDirPath / 'mash_output_GC.tab'
        print("Calculationg distance...")
        mashDistance(
            sketchFile,
            distanceTable,
            nthreads=args.cpus,
        )
        print('Gather families and calculating medoid...')
        dict_medoids, family_distance_matrice = \
            calculate_medoid(distanceTable, 0.8)

        familyGbksDirs = []
        for fasta, members in dict_medoids.items():
            if len(members) > 1:
                assert fasta in members
                name = Path(fasta).with_suffix('').name
                targetDir = gbkFamiliesDirPath / name
                targetDir.mkdir(exist_ok=True)
                for m in members:
                    gbk = fastaToGbk[m]
                    target = targetDir / gbk.name
                    target.symlink_to(gbk.resolve())
                familyGbksDirs.append(targetDir)

        args.outputPath.mkdir(exist_ok=True)
        cpuDist = (1,1)
        with Pool(round(args.cpus * cpuDist[0]/sum(cpuDist))) as bigscapePool:
            results = [
                bigscapePool.apply_async(
                    runBigscape,
                    (dir, args.outputPath / dir.name),
                    kwds={
                        'cpus': round(args.cpus * cpuDist[1]/sum(cpuDist)),
                        'cutoffs': [0.2, ]
                    }
                ) for dir in familyGbksDirs
            ]
            bigscapeResults = [
                r.get() for r
                in tqdm(results, desc='BiGSCAPE runs')
            ]

        print('finish')  # break point
    finally:
        fastaDir.cleanup()
        mashDir.cleanup()
        gbkFamiliesDir.cleanup()


if __name__ == "__main__":
    main()
