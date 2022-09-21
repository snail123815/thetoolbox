from pyBioinfo.wrappers.busco import runBusco
from pathlib import Path
from multiprocessing import Pool
from tqdm import tqdm
from typing import Any

targetRoot = Path('~/gdata/MBT-collection/MBT-noRawData-collection')
targetDirs = [d for d in targetRoot.iterdir() if d.is_dir()
              and len(d.name.split('_')) == 1]

runnerPool = Pool(1)
# Do not run more than 1 process. TODO Check what is wrong, probably
# temp dir of busco causing problem.

results: list[tuple[Path, Any]] = []
for d in targetDirs:
    prokkaDir = sorted((d / 'annotation').glob('prokka*'),
                       key=lambda d: d.name, reverse=True)[0]
    targetProteome = list(prokkaDir.glob('*.faa'))[0]
    kwargs = {
        'targetProteome': targetProteome,
        'outName': prokkaDir.name + 'busco',
        'outPath': prokkaDir,
        'condaEnv': '~/genvs/quasan',
        'shell': 'zsh',
        'cpu': 20,
        'silent': True
    }
    results.append((
        d,
        runnerPool.apply_async(
            runBusco,
            kwds=kwargs
        )
    ))

runnerPool.close()
with tqdm(results) as tqdmBar:
    for d, res in tqdmBar:
        tqdmBar.set_description(str(d))
        res.get()
