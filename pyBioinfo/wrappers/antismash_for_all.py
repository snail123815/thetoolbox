import argparse
from importlib.metadata import files
import logging
from pathlib import Path
from multiprocessing import Pool
from tqdm import tqdm

from bioSequences.bio_seq_file_extensions \
    import GBK_EXTENSIONS, FNA_EXTENSIONS
from decompress import IMPLEMENTED_COMPRESSION_FORMATS
from antismash import runAntismash


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'p',
        help='Input folder containing gbk files (can be gz or xz)'
    )
    parser.add_argument(
        '--parrllel',
        type=int, default=8,
        help="Number of antismashes run in parallel"
    )
    parser.add_argument(
        '--threads', type=int, default=4,
        help="Number of threads for each antismash run"
    )
    parser.add_argument(
        '--taxon', type=str,
        choices=['bacteria', 'fungi'],
        default='bacteria'
    )
    parser.add_argument(
        '--completeness', type=int,
        choices=[1, 2, 3, 10],
        help=f'Completeness of antismash run',
        default=2
    )
    parser.add_argument(
        '--dry', action='store_true'
    )
    parser.add_argument(
        '--overwrite', action='store_true'
    )
    parser.add_argument(
        '--geneFinding', type=str,
        choices=[
            'glimmerhmm', 'prodigal', 'prodigal-m', 'auto', 'error'
        ],
        help= '--geneFinding-tool for antismash, except "auto".'
        + 'It means "none" when the input is gbk file'
        + ' which should contain annotation.\n'
        + 'The antismash logic based on the fact that gbk file canbe both'
        + 'with or without annotation. But sometimes it miss interprate that.'
        + 'Implementation of this program is that:'
        + 'If input is genbank file, use "none" if set to "auto"'
        + 'If input is DNA sequence file, use "prodigal" if set to "auto"',
        default='auto'
    )
    args = parser.parse_args()

    pathIn = Path(args.p.strip())
    pathOut = pathIn.parent / (pathIn.name + '-antismash')
    pathOut.mkdir(exist_ok=True)
    logging.basicConfig(filename=pathOut / 'antismash.log',
                        filemode='a', level='INFO')
    logging.info('#' * 100)

    files = getValidFiles(pathIn,
                          allowGeneFinding=(args.geneFinding != 'error'))

    runnerPool = Pool(args.parrllel)
    results = []
    for f in files:
        if f.suffix in IMPLEMENTED_COMPRESSION_FORMATS:
            output = pathOut/f.with_suffix('').stem
        else:
            output = pathOut/f.stem
        results.append(
            runnerPool.apply_async(
                runAntismash,
                kwds={
                    "inputFilePath": f,
                    "title": f.name,
                    "taxon": args.taxon,
                    "completeness": args.completeness,
                    "cpu": args.threads,
                    "output": output,
                    "silent": True,
                    "dry": args.dry,
                    "geneFinding": args.geneFinding,
                    'overwrite': args.overwrite,
                }
            )
        )
    runnerPool.close()
    returns = []
    for res in tqdm(results):
        returns.append(res.get())
    logging.info('#' * 100)
    returns.sort(key=lambda x: x.name)
    for ret in returns:
        logging.info(ret)
    logging.info('#' * 100)


def getValidFiles(path: Path, allowGeneFinding: bool) -> list[Path]:
    newFiles: list[Path] = []
    fileStems = []
    allowedExts = (GBK_EXTENSIONS + FNA_EXTENSIONS if allowGeneFinding
                   else GBK_EXTENSIONS)
    for f in [f for f in path.iterdir() if f.is_file()]:
        if f.suffix in IMPLEMENTED_COMPRESSION_FORMATS:
            ext = f.with_suffix('').suffix
            stem = f.with_suffix('').stem
        else:
            ext = f.suffix
            stem = f.stem
        if ext in allowedExts and stem not in fileStems:
            newFiles.append(f)
            fileStems.append(stem)
    return newFiles


if __name__ == "__main__":
    main()
