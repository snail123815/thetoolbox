import subprocess
from pathlib import Path
from typing import Literal
from datetime import datetime
import os
from _environment_settings import \
    CONDAEXE, ANTISMASH_ENV, SHELL, getActivateEnvCmd
from decompress import decompFileIfCompressed, IMPLEMENTED_COMPRESSION_FORMATS
from bioSequences.bio_seq_file_extensions import FNA_EXTENSIONS


def getArgs():
    import argparse
    parser = argparse.ArgumentParser(
        description='Run antismash for genbank files'
    )
    parser.add_argument('--ncpu', type=int, default=4)
    parser.add_argument('--completeness',
                        type=int, default=2)
    parser.add_argument('--dry', action='store_true')
    parser.add_argument('--geneFinding', type=str, default='none')
    parser.add_argument('files', type=Path, nargs="+")
    return parser.parse_args()


def runAntismash(
    inputFilePath: Path,
    title: str | None = None,
    description: str | None = None,
    taxon: Literal['bacteria', 'fungi'] = 'bacteria',
    completeness: Literal[1, 2, 10] = 2,
    condaExe: Literal['conda', 'mamba', 'micromamba'] = CONDAEXE,
    condaEnv: Path | None = ANTISMASH_ENV,
    cpu: int = 4,
    output: Path | None = None,
    shell: Literal['bash', 'zsh'] = SHELL,
    prefix: str = 'antismash',
    geneFinding: Literal[
        'glimmerhmm', 'prodigal', 'prodigal-m', 'none', 'error'
    ] = 'error',
    defaultGeneFinding: str = 'prodigal',
    silent: bool = False,
    dry: bool = False,
    overwrite: bool = False
) -> Path:

    if not silent:
        print(f'Running antiSMASH for {inputFilePath}')

    inputFilePath, unzip = decompFileIfCompressed(inputFilePath)

    try:
        if output is None:
            # timeStr = datetime.now().strftime(r'%Y%m%d%H%M')
            prefix = ("_".join(item for item in
                               [
                                prefix,
                                title,
                                f'level{completeness}',
                                # timeStr
                               ]
                               if item is not None))
            outdir = inputFilePath.parent / prefix
        else:
            outdir = output
        if (outdir / 'index.html').exists():
            if overwrite:
                outdir.rmdir()
            else:
                raise FileExistsError(str(outdir))
        elif outdir.exists():
            outdir.rmdir()

        cmd = (f'antismash --cpus {cpu}'
               + ' --minimal'
               + ' --skip-zip-file'
               + f' --taxon {taxon}'
               + f' --html-title {prefix}'
               + f' --output-dir {outdir}')
        if inputFilePath.suffix in FNA_EXTENSIONS and geneFinding == 'none':
            cmd += f' --genefinding-tool {defaultGeneFinding}'
        else:
            cmd += f' --genefinding-tool {geneFinding}'
        if description is not None:
            cmd += f' --html-description {description}'

        if completeness >= 2:
            cmd = cmd.replace(' --minimal', '')
            cmd += ' --cb-knownclusters'
            cmd += ' --cb-subclusters'
            cmd += ' --asf'
        if completeness >= 3:
            cmd += ' --cb-general'
            cmd += ' --cc-mibig'
            cmd += ' --clusterhmmer'
            cmd += " --pfam2go"
            if taxon == 'fungi':
                cmd += ' --cassis'
        if completeness >= 4:
            cmd += ' --rre'
            cmd += ' --fullhmmer'
            cmd += ' --tigrfam'
            cmd += ' --smcog-trees'

        cmd += f' {inputFilePath}'

        if not silent:
            print(cmd)

        cmd = ' && '.join([
            getActivateEnvCmd(condaEnv, condaExe, shell),
            cmd
        ])

        if dry:
            print(cmd)
        else:
            commandResult = subprocess.run(
                cmd, capture_output=True, shell=True,
                executable=shell
            )
            if commandResult.returncode != 0:
                print('Failed antismash:')
                print(cmd)
                print(commandResult.stdout.decode())
                print(commandResult.stderr.decode())
    finally:
        if unzip:
            os.remove(str(inputFilePath))

    return outdir.resolve()


def main():
    args = getArgs()
    for f in args.files:
        if f.suffix in IMPLEMENTED_COMPRESSION_FORMATS:
            prefix = f.with_suffix('').stem + '_antismash'
        else:
            prefix = f.stem + '_antismash'
        runAntismash(f, cpu=args.ncpu,
                     dry=args.dry,
                     geneFinding=args.geneFinding,
                     prefix=prefix,
                     completeness=args.completeness)


if __name__ == "__main__":
    main()
