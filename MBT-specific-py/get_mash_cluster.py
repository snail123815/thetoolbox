# get result from cluster_BGCs_mash_bigscape.py
# Use the temporary file generated:
# clusters_info_dict.pickle
# family_dist_matrix_dict.pickle
# family_medoid_dict.pickle

# get which family, get number of members, get distance matrix
# make dendrogram based on distance matrix

from operator import index
import pickle
import argparse
from pathlib import Path
from pyBioinfo_modules.wrappers.antismash import ClusterInfo
import pandas as pd
from typing import cast
from pyBioinfo_modules.wrappers.bigscape \
    import runBigscape


def findFamily(
    dict_medoids: dict[str, list[str]],
    clusterInfoDict: dict,
    clustergbkname: str
) -> tuple[str, ClusterInfo, str, dict[str, ClusterInfo]]:
    for id, info in clusterInfoDict.items():
        if info['gbkFile'].name == clustergbkname:
            clusterId: str = id
            clusterInfo: ClusterInfo = info
            break
    assert "clusterId" in locals(), f'Cluster {clustergbkname} not found.'

    if clusterId in dict_medoids.keys():
        medoidId = clusterId
    else:
        for id in dict_medoids:
            if clusterId in dict_medoids[id]:
                medoidId = id
                break
    familyClusterDict = dict([(id, info) for id, info
                              in clusterInfoDict.items()
                              if id in dict_medoids[medoidId]])

    return clusterId, clusterInfo, medoidId, familyClusterDict


def getBigscapeNetwork(
    bigscapeNetworkDir: Path
) -> dict[str, dict[str, dict[str, Path | None]]]:
    """
    bigscapeNetworkDir:
    eg. bigscape/network_files/2022-10-10_01-46-35_hybrids_auto

    bigscapeNetworkPathsDict data structure to hold file Paths,
    dict{
        bgcType,
        dict{
            cutoff1,
            dict{
                "network": type_cutoff1.network,
                "clustering": type_clustering_cutoff1.tsv, # None if not found
            },
            cutoff2,
            dict{
                "network": type_cutoff2.network,
                "clustering": type_clustering_cutoff2.tsv,
            },
            "clans",
            dict{ # empty if not found
                "0.30_0.70": type_clans_0.30_0.70.tsv,
            }
        },
        bgcType,
        dict{
            ...
        }
    }
    """
    bigscapeNetworkPathsDict: dict[str, dict[str, dict[str, Path | None]]] = {}

    bgcNetworkByTypesDirs: list[Path] = [
        d for d in bigscapeNetworkDir.iterdir() if d.is_dir()
    ]
    for bgcNetworkDir in bgcNetworkByTypesDirs:
        bgcType = bgcNetworkDir.name
        networkFiles = [f for f in bgcNetworkDir.iterdir()
                        if f.suffix == '.network']
        bigscapeNetworkPathsDict[bgcType] = {}
        for networkFile in networkFiles:
            cutoff: str = networkFile.stem.split('_')[-1]
            assert cutoff.startswith('c')
            clusteringFile: Path | None = None
            try:
                clusteringFile = next(
                    bgcNetworkDir.glob(f'*_clustering_{cutoff}.tsv')).resolve()
            except StopIteration:
                pass
            networkFile = networkFile.resolve()

            bigscapeNetworkPathsDict[bgcType][cutoff] = {
                "network": networkFile,
                "clustering": clusteringFile
            }
        clanFiles = list(bgcNetworkDir.glob(f'*_clans_*_*.tsv'))
        bigscapeNetworkPathsDict[bgcType]['clans'] = {}
        for clanFile in clanFiles:
            cutoff = '_'.join(clanFile.stem.split('_')[-2:])
            bigscapeNetworkPathsDict[bgcType]['clans'][cutoff] = clanFile.resolve(
            )
    return bigscapeNetworkPathsDict


def bigscapeFamilyNetworkByNr(
    familyNr: int,
    clusteringDf: pd.DataFrame,
    networkDf: pd.DataFrame
) -> pd.DataFrame:
    familyBgcs = clusteringDf[clusteringDf.iloc[:, 0]
                              == familyNr].index.to_list()
    connectedDf = networkDf[networkDf.index == familyBgcs[0]]
    connectedDf = pd.concat(
        [connectedDf, networkDf[networkDf.iloc[:, 0] == familyBgcs[0]]], axis=0)
    if len(familyBgcs) > 1:
        for bgc in familyBgcs[1:]:
            connectedDf = pd.concat(
                [connectedDf, networkDf[networkDf.index == bgc]], axis=0)
            connectedDf = pd.concat(
                [connectedDf, networkDf[networkDf.iloc[:, 0] == bgc]], axis=0)
    connectedDf = connectedDf.reset_index()\
        .drop_duplicates().set_index(connectedDf.index.name)
    return connectedDf


def bigscapeRelatedBgcNetworkByName(
    bigscapeNetworkPathsDict: dict,
    targetBgc: str,
    notRealClanFilter: tuple[int, float] = (50, 0.9)
) -> dict[str, pd.DataFrame]:
    if targetBgc.endswith('.gbk'):
        targetBgc = targetBgc[:-4]
    bigscapeRelatedBgcsNetworkDict = {}

    for bgcType, networkByCutoff in bigscapeNetworkPathsDict.items():
        maxCutoff = sorted(list(networkByCutoff.keys()))[-2]
        networkDf = pd.read_csv(networkByCutoff[maxCutoff]['network'],
                                sep='\t', header=0, index_col=0)
        jointClusteringDf = pd.DataFrame()

        for cutoff, networkDict in networkByCutoff.items():
            if cutoff == 'clans':
                continue
            clusteringFile = networkDict['clustering']

            connectedDf = pd.DataFrame()
            if clusteringFile is not None:
                familyNr: int | None = None
                clusteringDf = pd.read_csv(
                    clusteringFile, sep='\t', header=0, index_col=0)
                jointClusteringDf = pd.concat([
                    jointClusteringDf, clusteringDf
                ])
                try:
                    familyNr = list(clusteringDf.loc[targetBgc, :])[0]
                except KeyError:
                    pass
                if familyNr is not None:
                    connectedDf = pd.concat(
                        [
                            connectedDf,
                            bigscapeFamilyNetworkByNr(
                                familyNr, clusteringDf, networkDf)
                        ],
                        axis=0
                    )
                bigscapeRelatedBgcsNetworkDict[
                    '_'.join([bgcType, cutoff])] = connectedDf.copy()

        for clanCutoff, clansFile in networkByCutoff['clans'].items():
            clanNr: int | None = None
            clansDf = pd.read_csv(clansFile, sep='\t', index_col=0, header=0)
            try:
                clanNr = clansDf.loc[targetBgc, :].iloc[0]
            except KeyError:
                pass

            clanConnectedDf = pd.DataFrame()
            if clanNr is not None:
                familyNrs = clansDf[clansDf.iloc[:, 0] == clanNr]\
                    .iloc[:, 1].unique()
                clanBgcCount = clansDf.iloc[:, 0].value_counts()[clanNr]
                if clanBgcCount > notRealClanFilter[0] and \
                        len(familyNrs) / len(clansDf.iloc[:, 1].unique()) \
                        > notRealClanFilter[1]:
                    # The clan is possibly a place holder.
                    print(f'''
                    Clan {clanNr} is for not clanned
                    (filter: BGC count > {notRealClanFilter[0]} and
                    families in this clan > {notRealClanFilter[1]:.1%} of total
                    number of families ({len(clansDf.iloc[:, 1].unique())}))
                    ''')
                    familyNrs = None

                if familyNrs is not None:
                    # process jointCLusteringDf first
                    jointClusteringDf = jointClusteringDf\
                        .reset_index().drop_duplicates()\
                        .set_index(jointClusteringDf.index.name)
                    for familyNr in familyNrs:
                        clanConnectedDf = pd.concat(
                            [
                                clanConnectedDf,
                                bigscapeFamilyNetworkByNr(
                                    cast(int, familyNr),
                                    jointClusteringDf, networkDf)
                            ],
                            axis=0
                        )
                    clanConnectedDf = clanConnectedDf.reset_index()\
                        .drop_duplicates().set_index(clanConnectedDf.index.name)
            bigscapeRelatedBgcsNetworkDict[
                '_'.join([bgcType, f'clan{clanNr}', clanCutoff])
            ] = clanConnectedDf

    return bigscapeRelatedBgcsNetworkDict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--family_medoid_dict', type=Path
    )
    parser.add_argument(
        '--family_dist_matrix_dict', type=Path
    )
    parser.add_argument(
        '--cluster_info_dict', type=Path
    )
    parser.add_argument(
        '--target_cluster', type=str
    )
    parser.add_argument(
        '--bigscape_network_files', type=Path
    )
    parser.add_argument(
        '--bgcsRoot', type=Path
    )
    parser.add_argument(
        '--out', type=Path
    )
    args = parser.parse_args()

    with open(args.cluster_info_dict, 'rb') as fh:
        clusterInfoDict: dict = pickle.load(fh)
    with open(args.family_dist_matrix_dict, 'rb') as fh:
        family_distance_matrice: dict[str, list[list[float]]] = pickle.load(fh)
    with open(args.family_medoid_dict, 'rb') as fh:
        dict_medoids: dict[str, list[str]] = pickle.load(fh)

    # From mash results
    clusterId, clusterInfo, medoidId, familyClusterDict = findFamily(
        dict_medoids, clusterInfoDict, args.target_cluster
    )
    distance_matrix = family_distance_matrice[medoidId]

    medoidBgcName = familyClusterDict[medoidId]['gbkFile'].name

    # From BiGSCAPE
    bigscapeNetworkPathsDict = getBigscapeNetwork(args.bigscape_network_files)
    collectedBgcPaths: list[Path] = []
    for info in familyClusterDict.values():
        targetBgc = info['gbkFile'].name
        bigscapeRelatedBgcsNetworkDict = bigscapeRelatedBgcNetworkByName(
            bigscapeNetworkPathsDict, targetBgc)
        if not any(df.shape != (0, 0)
                   for df in bigscapeRelatedBgcsNetworkDict.values()):
            continue

        # gather related BGC paths
        collectedBgcNetworkDf = pd.concat(
            [bigscapeRelatedBgcsNetworkDict.values()], axis=0).reset_index()\
            .drop_duplicates().set_index(collectedBgcNetworkDf.index.name)
        collectedBgcs = pd.concat([
            collectedBgcNetworkDf.index.to_series(),
            collectedBgcNetworkDf.iloc[:, 0]
        ], axis=0).drop_duplicates().values
        collectedBgcPaths.extend([
            args.bgcsRoot.parent / bgcInfo['gbkFile']
            for bgcInfo in clusterInfoDict.values()
            if bgcInfo['gbkFile'].stem in collectedBgcs
        ])
        assert all(bgcPath.is_file() for bgcPath in collectedBgcPaths)
        print(
            f'Got {len(collectedBgcs)} bgcs related with {targetBgc}')

    assert len(collectedBgcPaths) > 0, \
        "Target cluster not found in BiGSCAPE results"

    print(f'Running BiGSCAPE')
    collectedBgcGatherPath: Path = args.out / 'bgcs'
    collectedBgcGatherPath.mkdir(parents=True, exist_ok=True)
    bigscapeOut = args.out / 'bigscape'
    for f in collectedBgcPaths:
        target = collectedBgcGatherPath / f.name
        if not target.is_file():
            target.symlink_to(f)
    runBigscape(
        collectedBgcGatherPath,
        bigscapeOut,
        cpus=20
    )

    print('FINISH')


if __name__ == "__main__":
    main()
