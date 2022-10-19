# get result from cluster_BGCs_mash_bigscape.py
# Use the temporary file generated:
# clusters_info_dict.pickle
# family_dist_matrix_dict.pickle
# family_medoid_dict.pickle

# get which family, get number of members, get distance matrix
# make dendrogram based on distance matrix

import pickle
import argparse
from pathlib import Path
import pandas as pd
from pyBioinfo_modules.wrappers.bigscape \
    import runBigscape
from pyBioinfo_modules.mash_cluster_bgc_bigscape.process_bigscape_result \
    import getBigscapeNetwork, bigscapeRelatedBgcNetworkByName, findFamily


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
