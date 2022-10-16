import argparse
import os
from multiprocessing import Pool
from glob import glob
from tqdm import tqdm
from pathlib import Path
from pyBioinfo_modules.wrappers.antismash \
    import clusterGbkGlobTxt, parseClusterGbk
from pyBioinfo_modules.wrappers.mash \
    import calculate_medoid, mashSketchFiles, mashDistance
from pyBioinfo_modules.wrappers.bigscape \
    import runBigscape
from tempfile import TemporaryDirectory
import pickle


def globFilesSafely(
    sourceDir: Path, globPattern: str,
    resultPickle: Path | None = None,
    reload: bool = False, showProgress: bool = False
) -> list[Path]:
    """Make sure glob files will follow symlinks
    In pathlib.Path, the glob() function do not follow symlink:
    [pathlib.Path.glob does not follow symlinks · Issue #77609 · python/cpython
    (github.com)](https://github.com/python/cpython/issues/77609)
    But this is hardly a bug because this is a system behaviour.
    However this will happily follow symlinks:
    ```python
    from glob import glob
    glob(str, recursive=True)
    ```
    """
    if resultPickle is not None and resultPickle.is_file() and not reload:
        with resultPickle.open('rb') as fh:
            return pickle.load(fh)
    else:
        files = list(sourceDir.glob(globPattern))
        if showProgress:
            for d in tqdm(
                [d for d in sourceDir.iterdir() if d.is_dir()],
                desc=(f'Gathering "{globPattern}" files from {sourceDir.name}')
            ):
                files.extend(Path(f) for f in
                             glob(str(d / ('**/' + globPattern)),
                                  recursive=True))
        else:
            for d in [d for d in sourceDir.iterdir() if d.is_dir()]:
                files.extend(Path(f) for f in
                             glob(str(d / ('**/' + globPattern)),
                                  recursive=True))
        if resultPickle is not None:
            with resultPickle.open('wb') as fh:
                pickle.dump(files, fh)
    return files


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
        gbks = globFilesSafely(
            inputPath, clusterGbkGlobTxt,
            (gbksListFile if args.tmp is not None else None),
            showProgress=True
        )

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
