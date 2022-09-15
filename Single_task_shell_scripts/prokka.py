import subprocess
from pathlib import Path
from typing import Literal
from datetime import datetime
import os


def getArgs():
    import argparse
    parser = argparse.ArgumentParser(description='Run prokka for fasta files')
    parser.add_argument('--ncpu', type=int, default=4)
    parser.add_argument('fastas', type=str, nargs="+", required=True)
    return parser.parse_args()


def runProkka(
    fastaPath: Path,
    gcode: int = 11,
    gram: Literal['pos', 'neg'] = 'pos',
    center: str | None = None,
    genus: str | None = None,
    species: str | None = None,
    strain: str | None = None,
    locustag: str | None = None,
    condaEnv: Path | None = None,
    cpu: int = 4,
    shell: Literal['bash', 'zsh'] = 'zsh',
    output: Path = Path("."),
    prefix: str = 'prokka'
) -> Path:

    print(f'Running prokka for {fastaPath}')

    unzip = False
    if fastaPath.suffix == '.gz':
        gzip = subprocess.run(f'gzip -dkf {fastaPath}'.split(' '),
                              capture_output=True)
        assert gzip.returncode == 0
        fastaPath = fastaPath.with_suffix('')
        assert fastaPath.exists(), '\n'.join([
            f'Unzip file {fastaPath} failed with error message:\n',
            gzip.stderr.decode(), gzip.stdout.decode()
        ])
        unzip = True

    timeStr = datetime.now().strftime(r'%Y%m%d%H%M')
    outdir = output/("_".join(item for item in
                              [prefix, genus, species, strain, timeStr]
                              if item is not None))
    cmd = 'prokka --compliant --addgenes --mincontiglen 200 --rfam' + \
        f' --gcode {gcode}' + \
        f' --gram {gram}' + \
        f' --cpu {cpu}' + \
        f' --outdir {outdir}' + \
        (f' --centre {center}' if center is not None else "") + \
        (f' --genus {genus}' if genus is not None else "") + \
        (f' --strain {strain}' if strain is not None else "") + \
        (f' --species {species}' if species is not None else "") + \
        (f' --locustag {locustag}' if locustag is not None else "") + \
        f' {fastaPath}'

    activateEnvCmd = f'eval "$(micromamba shell hook --shell={shell})"' + \
        f' && micromamba activate {condaEnv}' + \
        f' && {cmd}'
    commandResult = subprocess.run(
        activateEnvCmd, capture_output=True, shell=True,
        executable=shell
    )
    if commandResult.returncode != 0:
        print('Failed prokka:')
        print(cmd)
        print(commandResult.stdout.decode())
        print(commandResult.stderr.decode())

    if unzip:
        os.remove(str(fastaPath))

    return outdir.resolve()


def main():
    runProkka(Path('abc'), condaEnv='~/genvs/quasan')
    pass


if __name__ == "__main__":
    main()
