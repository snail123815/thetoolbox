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
from pyBioinfo_modules.wrappers.antismash import ClusterInfo


def findFamily(
    dict_medoids: dict[str, list[str]],
    clusterInfoDict: dict,
    clustergbkname: str
) -> tuple[str, ClusterInfo, dict[str, ClusterInfo]]:
    for id, info in clusterInfoDict.items():
        if info['gbkFile'].name == clustergbkname:
            clusterId: str = id
            clusterInfo: ClusterInfo = info
            break
    assert "clusterId" in locals(), f'Cluster {clustergbkname} not found.'
    familyClusterDict = dict([(id, info) for id, info
                              in clusterInfoDict.items()
                              if id in dict_medoids[clusterId]])
    return clusterId, clusterInfo, familyClusterDict


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
    args = parser.parse_args()

    with open(args.cluster_info_dict, 'rb') as fh:
        clusterInfoDict: dict = pickle.load(fh)
    with open(args.family_dist_matrix_dict, 'rb') as fh:
        family_distance_matrice: dict[str, list[list[float]]] = pickle.load(fh)
    with open(args.family_medoid_dict, 'rb') as fh:
        dict_medoids: dict[str, list[str]] = pickle.load(fh)

    clusterId, clusterInfo, familyClusterDict = findFamily(
        dict_medoids, clusterInfoDict, args.target_cluster
    )
    distance_matrix = family_distance_matrice[clusterId]

    print(distance_matrix)


if __name__ == "__main__":
    main()
