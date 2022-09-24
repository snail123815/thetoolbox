import argparse
import shutil
from pathlib import Path
from pyBioinfo_modules.wrappers.antismash import clusterGbkGlobTxt

parser = argparse.ArgumentParser(prog='''
Gather antismash results from run_antismash.py
Name of each folder in the ROOT path will be used to name
the .gbk files of antismash result.
''')
parser.add_argument('p', type=Path, help="ROOT path to antismash results")
parser.add_argument('--out', type=Path,
                    help="Output root path", default=Path(''))
parser.add_argument('--link', action='store_true')
args = parser.parse_args()

outputRoot = (args.p.parent / (args.p.name + '_collection_gbks')
              if args.out == Path('') else args.out)

dirs = [d for d in args.p.iterdir() if d.is_dir()]
for d in dirs:
    outputDir = outputRoot / d.name
    outputDir.mkdir(parents=True, exist_ok=True)
    for gbk in d.glob(clusterGbkGlobTxt):
        if d.name not in gbk.name:
            target = outputDir / '_'.join([d.name, gbk.name])
        else:
            target = outputDir / gbk.name
        if args.link:
            target.symlink_to(gbk.resolve())
        else:
            shutil.copyfile(gbk, target, follow_symlinks=True)