import argparse
import logging
from typing import Literal
from pathlib import Path
from multiprocessing import Pool

from ..bioSequences.bio_seq_file_extensions \
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
        '--processes',
        type=int, default=8,
        help="Number of antismashes run in parallel"
    )
    parser.add_argument(
        '--threads', type=int, default=4,
        help="Number of threads for each antismash run"
    )
    parser.add_argument(
        '--taxon', type=str,
        default='bacteria'
    )
    parser.add_argument(
        '--completeness', type=Literal[1, 2, 10],
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
        '--geneFinding', type=Literal[
            'glimmerhmm', 'prodigal', 'prodigal-m', 'none', 'error'
        ],
        help=', '.join(['glimmerhmm', 'prodigal',
                       'prodigal-m', 'none', 'error']),
        default='none'
    )
    args = parser.parse_args()

    pathIn = Path(args.p.strip())
    pathOut = pathIn.parent / (pathIn.name + '-antismash')
    pathOut.mkdir(exist_ok=True)
    logging.basicConfig(filename=pathOut / 'antismash.log',
                        filemode='a', level='INFO')

    files = getValidFiles(pathIn, allowGeneFinding=args.allowGeneFinding)

    runnerPool = Pool(args.processes)
    results = []
    for f in files:
        results.append(
            runnerPool.apply_async(
                runAntismash,
                kwds={
                    "inputFilePath": f,
                    "title": f.name,
                    "taxon": args.taxon,
                    "completeness": args.completeness,
                    "cpu": args.threads,
                    "output": pathOut,
                    "silent": True,
                    "dry": args.dry,
                    "geneFindingTool": args.geneFinding,
                    'overwrite': args.overwrite,
                }
            )
        )
    runnerPool.close()
    returns = []
    if args.noProgress:
        for res in results:
            returns.append(res.get())
    else:
        from tqdm import tqdm
        for res in tqdm(results):
            returns.append(res.get())
    logging.info('#' * 100)
    returns.sort(key=lambda x: x[0])
    for ret in returns:
        logging.info(ret)


def getValidFiles(path: Path, allowGeneFinding: bool) -> list[Path]:
    newFiles: list[Path] = []
    allowedExts = (GBK_EXTENSIONS + FNA_EXTENSIONS if allowGeneFinding
                   else GBK_EXTENSIONS)
    for f in [f for f in path.iterdir() if f.is_file()]:
        if f.suffix in IMPLEMENTED_COMPRESSION_FORMATS:
            ext = f.with_suffix('').suffix
        else:
            ext = f.suffix
        if ext in allowedExts:
            newFiles.append(f)
    return newFiles


if __name__ == "__main__":
    main()
