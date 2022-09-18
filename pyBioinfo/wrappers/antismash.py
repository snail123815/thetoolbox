import subprocess
import logging
from pathlib import Path
from typing import Literal
import shutil
from datetime import datetime
import os
import sys
from _environment_settings import \
    CONDAEXE, ANTISMASH_ENV, SHELL, getActivateEnvCmd
from decompress import decompFileIfCompressed, IMPLEMENTED_COMPRESSION_FORMATS
from bioSequences.bio_seq_file_extensions import FNA_EXTENSIONS


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Run antismash for genbank files'
    )
    parser.add_argument('--ncpu', type=int, default=4)
    parser.add_argument(
        '--taxon', type=str,
        choices=['bacteria', 'fungi'],
        default='bacteria'
    )
    parser.add_argument('--completeness',
                        type=int,
                        choices=[1, 2, 3, 10],
                        default=2)
    parser.add_argument('--dry', action='store_true')
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
        default='error'
    )
    parser.add_argument('files', type=Path, nargs="+")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format='%(message)s',
                        handlers=[logging.StreamHandler(sys.stdout)])
    for f in args.files:
        if f.suffix in IMPLEMENTED_COMPRESSION_FORMATS:
            prefix = f.with_suffix('').stem + '_antismash'
        else:
            prefix = f.stem + '_antismash'
        runAntismash(f, cpu=args.ncpu,
                     dry=args.dry,
                     taxon=args.taxon,
                     geneFinding=args.geneFinding,
                     prefix=prefix,
                     completeness=args.completeness)


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
    addDateTimeToPrefix: bool = False,
    geneFinding: Literal[
        'glimmerhmm', 'prodigal', 'prodigal-m', 'auto', 'error'
    ] = 'error',
    defaultGeneFinding: str = 'prodigal',
    silent: bool = False,
    dry: bool = False,
    overwrite: bool = False
) -> Path:

    logging.info(f'Running antiSMASH for {inputFilePath}')

    inputFilePath, unzip = decompFileIfCompressed(inputFilePath)

    try:
        if output is None:
            prefix = ("_".join(item for item in
                               [
                                   prefix,
                                   title,
                                   f'level{completeness}',
                               ]
                               if item is not None))
            if addDateTimeToPrefix:
                timeStr = datetime.now().strftime(r'%Y%m%d%H%M')
                prefix += "_" + timeStr
            outdir = inputFilePath.parent / prefix
        else:
            outdir = output
        if (outdir / 'index.html').exists():
            if overwrite:
                shutil.rmtree(outdir)
            else:
                raise FileExistsError(str(outdir))
        elif outdir.exists():
            shutil.rmtree(outdir)

        cmd = (f'antismash --cpus {cpu}'
               + ' --minimal'
               + ' --skip-zip-file'
               + f' --taxon {taxon}'
               + f' --html-title {prefix}'
               + f' --output-dir {outdir}')
        if inputFilePath.suffix in FNA_EXTENSIONS and geneFinding == 'auto':
            cmd += f' --genefinding-tool {defaultGeneFinding}'
        elif geneFinding == 'auto':
            cmd += f' --genefinding-tool none'
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
            logging.info(cmd)

        cmd = ' && '.join([
            getActivateEnvCmd(condaEnv, condaExe, shell),
            cmd
        ])

        if dry:
            logging.info(cmd)
        else:
            commandResult = subprocess.run(
                cmd, capture_output=True, shell=True,
                executable=shell
            )
            if commandResult.returncode != 0:
                logging.error('Failed antismash:')
                logging.error(cmd)
                logging.error(commandResult.stdout.decode())
                logging.error(commandResult.stderr.decode())
    finally:
        if unzip:
            os.remove(str(inputFilePath))
    logging.info(f'Done antiSMASH for {inputFilePath}')

    return outdir.resolve()


if __name__ == "__main__":
    main()
