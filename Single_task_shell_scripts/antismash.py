import subprocess
from pathlib import Path
from typing import Literal
from datetime import datetime
import os


def getArgs():
    import argparse
    parser = argparse.ArgumentParser(
        description='Run antismash for genbank files'
    )
    parser.add_argument('--ncpu', type=int, default=4)
    parser.add_argument('genbankfiles', type=str, nargs="+", required=True)
    return parser.parse_args()


def runAntismash(
    genbankFilePath: Path,
    title: str | None = None,
    description: str | None = None,
    taxon: Literal['bacteria', 'fungi'] = 'bacteria',
    completeness: Literal[1, 2, 3] = 2,
    condaEnv: Path | None = None,
    cpu: int = 4,
    output: Path | None = None,
    shell: Literal['bash', 'zsh'] = 'zsh',
    prefix: str = 'antismash',
    silent: bool = False
) -> Path:

    if not silent:
        print(f'Running antiSMASH for {genbankFilePath}')

    unzip = False
    if genbankFilePath.suffix == '.gz':
        gzip = subprocess.run(f'gzip -dkf {genbankFilePath}'.split(' '),
                              capture_output=True)
        assert gzip.returncode == 0
        genbankFilePath = genbankFilePath.with_suffix('')
        assert genbankFilePath.exists(), '\n'.join([
            f'Unzip file {genbankFilePath} failed with error message:\n',
            gzip.stderr.decode(), gzip.stdout.decode()
        ])
        unzip = True

    try:
        timeStr = datetime.now().strftime(r'%Y%m%d%H%M')
        prefix = ("_".join(item for item in
                           [prefix, title, f'level{completeness}', timeStr]
                           if item is not None))
        if output is None:
            output = genbankFilePath.parent
        outdir = output/prefix
        cmd = (f'antismash --cpus {cpu} --genefinding-tool none'
               + f' --taxon {taxon}'
               + f' --html-title {prefix}'
               + f' --output-dir {outdir}')
        if description is not None:
            cmd += f' --html-description {description}'

        if completeness > 1:
            cmd += ' --cb-general --cb-knownclusters --cb-subclusters'
            cmd += ' --clusterhmmer'
            if taxon == 'fungi':
                cmd += ' --cassis'
        if completeness > 2:
            cmd += ' --asf'
            cmd += ' --rre'  # needs fimo
            cmd += ' --fullhmmer'
            cmd += ' --tigrfam'
            cmd += ' --smcog-trees'

        cmd += f' {genbankFilePath}'

        activateEnvCmd = (
            f'eval "$(micromamba shell hook --shell={shell})"'
            + f' && micromamba activate {condaEnv}'
            + f' && {cmd}'
        )
        commandResult = subprocess.run(
            activateEnvCmd, capture_output=True, shell=True,
            executable=shell
        )
        if commandResult.returncode != 0:
            print('Failed antismash:')
            print(cmd)
            print(commandResult.stdout.decode())
            print(commandResult.stderr.decode())
    finally:
        if unzip:
            os.remove(str(genbankFilePath))

    return outdir.resolve()


def main():
    runAntismash(Path('abc'), condaEnv='~/genvs/quasan')
    pass


if __name__ == "__main__":
    main()
