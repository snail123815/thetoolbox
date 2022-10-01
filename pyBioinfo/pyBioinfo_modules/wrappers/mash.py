import sys
import subprocess
from pathlib import Path
from typing import Literal
import numpy as np
from tqdm import tqdm
from tempfile import NamedTemporaryFile

from pyBioinfo_modules.wrappers._environment_settings \
    import CONDAEXE, SHELL, MASH_ENV, withActivateEnvCmd


def mashSketchFiles(
    inputFiles: list[Path],
    output: Path,
    kmer: int,
    sketch: int,
    nthreads: int = 1,
    mashEnv=MASH_ENV,
    condaExe=CONDAEXE,
    shell=SHELL,
    molecule: Literal['DNA', 'protein'] = "DNA"
) -> Path:
    """
    Calculates the distance between the query fasta files
    stored in the sketch file by using mash.
    """
    fileList = NamedTemporaryFile()
    with open(fileList.name, 'w') as fl:
        for f in inputFiles:
            fl.write(f'{f.resolve()}\n')
    cmd = f"mash sketch -o {output} -k {kmer} -p {nthreads} -s {sketch}"
    cmd += (' -a' if molecule == 'protein' else '')
    cmd += f' -l {fileList.name}'
    mashSketchRun = subprocess.run(
        withActivateEnvCmd(cmd, mashEnv, condaExe, shell),
        shell=True, capture_output=True, check=True)
    fileList.close()
    outputMsh = Path(str(output) + '.msh')
    assert outputMsh.is_file()
    return outputMsh


def mashDistance(
    inputMsh: Path,
    outputFile: Path,
    nthreads: int = 1,
    mashEnv=MASH_ENV,
    condaExe=CONDAEXE,
    shell=SHELL,
) -> Path:
    """
    Calculates the distance between the query fasta files
    stored in the sketch file by using mash.
    Parameters
    ----------
    outdir
        string, the path to output directory
    returns
        Path for the output file
   The output fields are
reference-ID  query-ID     distance  p-value  shared-hashes
genome1.fna   genome3.fna  0         0        1000/1000
genome2.fna   genome3.fna  0.022276  0        456/1000
    ----------
    """
    cmd = f"mash dist -p {nthreads} {inputMsh} {inputMsh} > {outputFile}"
    mashDistRun = subprocess.run(
        withActivateEnvCmd(cmd, mashEnv, condaExe, shell),
        shell=True, check=True)
    assert outputFile.exists
    return outputFile

######################################################################
# Extracting LSH clusters (buckets) using cut-off & calculate medoid
######################################################################


def calculate_medoid(
    inputDistanceTablePath: Path,  # output file (return) of mashDistance()
    cutOff: float,  # default 0.8
    med: dict[str, list[str]] = {}
) -> tuple[dict[str, list[str]], dict[str, list[float]]]:
    """
    calculates the GCFs based on similarity threshold
    parameters and calculates the medoid of that GCF
    ----------
    outdir
        string, the path to output directory
    cut_off
        float, between 0 and 1
    returns
    ----------
    dict_medoids = {fasta file of medoid: similar fasta files}
    """
    # Parse the input into a dictionary of gene families
    family: dict[str, str] = {}
    # family dict, key = family members, values = family names
    familyFiltered = {}
    family_members: dict[str, list[str]] = {}
    # family_members dict, key = family name, values = list of members
    family_distance_matrices: dict[str, list[float]] = {}
    dict_medoids: dict[str, list[str]] = med
    pbar = tqdm(total=inputDistanceTablePath.stat().st_size,
                bar_format=r"{l_bar}{bar}| {n:,.0f}/{total:,.0f} {unit} " +
                r"[{elapsed}<{remaining}, {rate_fmt}{postfix}]",
                unit_scale=1 / 1048576, unit='MB',
                desc="Reading mash dist table and generate families")
    with inputDistanceTablePath.open('r') as input:
        readSize = 0
        for idx, line in enumerate(input):
            readSize += sys.getsizeof(line) - 50  # length of '\n'
            if idx % 10000 == 0:
                pbar.update(readSize)
                readSize = 0
            if line.startswith('#') or line.strip() == "":
                continue
            # Split into tab-separated elements
            refId, queryId, distanceStr, pValueStr, nHashesStr \
                = line.strip().split('\t')
            sharedNhashes, totalNhashes = (int(n) for n in
                                           nHashesStr.split("/"))
            shareRatio = sharedNhashes / totalNhashes
            distance = float(distanceStr)
            pValue = float(pValueStr)
            # Look up the family of the first gene
            # Each family is named by the first gene of the family
            if queryId in family.keys():
                familyName = family[queryId]
            else:  # init a new family
                familyName = queryId
                family[queryId] = familyName
                family_members[familyName] = []
                family_distance_matrices[familyName] = []
            # filter genes that are already part of another family
            if queryId == family[queryId]:
                familyFiltered[queryId] = family[queryId]
            # If refId and queryId don't overlap, then there are two options:
            # 1. refId and queryId do belong to the same family, and the
            #    "not overlap" is "odd" in this case we accept queryId into
            #    our family
            # 2. refId and queryId actually belong to different families,
            #    and we might have to create that family
            if shareRatio <= cutOff:
                # refId doesn't overlap at all, so put that one into a separate
                # family
                if refId in family.keys():
                    if family[refId] == familyName:
                        # refId is in our family, so record the distance
                        add_to_distance_matrix(
                            family_distance_matrices[familyName],
                            family_members[familyName],
                            queryId,
                            refId,
                            distance
                        )
                    else:
                        # gene is above cut off or doesn't belong to the family
                        pass
                else:
                    # refId doesn't have a family yet, make one
                    gene1_family_name = refId
                    family[refId] = gene1_family_name
                    family_members[gene1_family_name] = []
                    family_distance_matrices[gene1_family_name] = []
                    # insert refId into that family as only member, with a
                    # distance of 0
                    add_to_distance_matrix(
                        family_distance_matrices[gene1_family_name],
                        family_members[gene1_family_name],
                        refId,
                        refId,
                        distance
                    )
            else:
                # There is some overlap, and we want refId also in this family
                family[refId] = familyName
                add_to_distance_matrix(
                    family_distance_matrices[familyName],
                    family_members[familyName],
                    queryId,
                    refId,
                    distance
                )
    pbar.close()
    # For each family: Build a distance matrix, and then work out the medoid
    for familyName in familyFiltered.keys():
        # Calculate the medoid from the distances
        np_array = np.asarray(family_distance_matrices[familyName])
        medoid_index = np.argmin(np_array.sum(axis=0))
        # Create a dictionary using the medoid as key and
        # the family_members as values
        dict_medoids[family_members[familyName][medoid_index]] \
            = family_members[familyName]
    return dict_medoids, family_distance_matrices


def add_to_distance_matrix(distance_matrix, idList,
                           refId, queryId, distance):
    """
    Adds the distance of refId-queryId to the distance matrix and index
    ----------
    distance_matrix
        {fasta file name: distance matrix}
    idList
        list, fasta file name
    refId
        fasta file name
    queryId
        fasta file name
    distance
        float, between 0 and 1
    returns
    ----------
    """
    index1 = add_new_gene(distance_matrix, idList, refId)
    index2 = add_new_gene(distance_matrix, idList, queryId)
    distance_matrix[index1][index2] = distance
    distance_matrix[index2][index1] = distance
    return ()


def add_new_gene(distance_matrix, idList, id) -> int:
    """
    Adds a distance matrix
    ----------
    distance_matrix
        {fasta file name: distance matrix}
    idList
        list, fasta file name
    id
        fasta file name
    returns
    ----------
    idList.index(gene)
        int, index number of the gene
    """
    if id not in idList:
        # Add to list
        idList.append(id)
        idList.index(id)
        # Extend distance matrix: One new row, one new column
        for row in distance_matrix:
            row.append(0)
        distance_matrix.append([0] * len(idList))
    return idList.index(id)
