# For RNA-Seq alignment, with each file or pair of files per sample. Out
# put aligned file per sample.


import subprocess
from tempfile import NamedTemporaryFile
from Bio import SeqIO
import os
import time
import argparse
import logging
from pathlib import Path
from pyBioinfo_modules.basic.decompress import splitStemSuffixIfCompressed
from pyBioinfo_modules.wrappers._environment_settings \
    import SHORTREADS_ENV, SHELL, withActivateEnvCmd
from typing import IO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--raw',
        required=True,
        type=Path,
        nargs='+',
        help='path(s) to raw data')
    parser.add_argument(
        '--out',
        required=True,
        type=Path,
        help='path to output alignment')
    parser.add_argument(
        '--genome',
        nargs="+",
        required=True,
        type=Path,
        help='path to genome file(s), also supports indexed genome (*.bt2)')
    parser.add_argument(
        '--isPe',
        action='store_true',
        help='set if you have pairend')

    # TODO this should be removed
    parser.add_argument('--rawExt', help='extension ".fastq.gz" or ".fq.gz"')
    parser.add_argument(
        '--pesuffix',
        nargs=2,
        help=('suffix to pairend file name. Will impute by default.'
              ' eg. a_1_.fq.gz and a_2_.fq.gz, --pesuffix _1_ _2_'))

    parser.add_argument(
        '--ncpu',
        default=1,
        type=int,
        help='number of cpu to use')
    parser.add_argument(
        '--sampleNames',
        nargs="+",
        help=('List of sample names.'
              ' Set when sample names are random string'
              ' which will unpair the samples after removing pairend suffix'))

    args = parser.parse_args()

    # TODO check existance of files before creating output dir
    if not args.out.is_dir():
        args.out.mkdir(exist_ok=True)

    toBamNcpu = max(args.ncpu // 8, 1)

    logging.basicConfig(
        filename=os.path.join(
            args.out / 'align.log'),
        level=logging.DEBUG)
    logging.debug(args)
    logging.info('=' * 20 + getTime() + '=' * 20)

    genomeBowtie2Idx = buildBowtie2idx(
        [args.genome], out=args.out / 'genomeIdx')

    # Gether files for each sample, max depth 2
    filePaths = [f for f in args.raw.iterdir() if f.is_file()]
    [filePaths.extend(f for f in d.iterdir() if f.is_file())
     for d in args.raw.iterdir() if d.is_dir()]
    for f in filePaths.copy():
        ext = splitStemSuffixIfCompressed(f)[1]
        if ext.split('.')[0] not in ['fastq', 'fq']:
            filePaths.pop(f)
    assert len(filePaths) > 0, f'Files not found in {args.raw}'
    FILE_PATHS = sorted(filePaths)
    FILE_PARTS = [splitStemSuffixIfCompressed(f, fullSuffix=True) for
                  f in FILE_PATHS]
    FILE_NAMES = [fp[0] for fp in FILE_PARTS]
    FILE_FULLEXTS = [fp[1] for fp in FILE_PARTS]
    assert len(set(FILE_FULLEXTS)) == 1, (
        'All files should have same format, '
        f'yet multiple found {set(FILE_FULLEXTS)}.')
    assert len(set(FILE_NAMES)) == len(FILE_NAMES), (
        'File names has duplicates'
        ', maybe from different dirs?\n'
        f'{set([fn for fn in FILE_NAMES if FILE_NAMES.count(fn) > 1])}')

    if args.sampleNames is not None:
        samples = args.sampleNames
    else:
        samples = []

    sampleFileDict: dict[str, list[Path]] = {}

    peSfx: list[str] = []
    if args.isPe:
        peSfx = imputePeSuffix(FILE_PATHS)
        # TODO make it only accecpt file stems
        if samples == []:
            samples = (sorted(list(set(f[:-len(peSfx[0])] for f in FILE_NAMES)))
                       if len(samples) == 0 else samples)
            assert (len(FILE_PATHS) / len(samples)) % 2 == 0
        for s in samples:
            fns = [FILE_PATHS[i] for i, fn in enumerate(FILE_NAMES) if s in fn]
            assert len(fns) > 0, f'Did not find files for sample {s}'
            assert len(fns) % 2 == 0, (
                f'Number of files for sample {s} is not paired.\n'
                f'{fns}'
            )
            sampleFileDict[s] = fns

    else:  # single end reads
        if samples == []:
            for fn, fp in zip(FILE_NAMES, FILE_PATHS):
                sampleFileDict[fn] = [fp]
            samples = FILE_NAMES
        elif len(samples) == 1:
            sampleFileDict[samples[0]] = FILE_PATHS
        else:
            for s in samples:
                fps = [FILE_PATHS[i] for i, fn in enumerate(FILE_NAMES)
                       if s in fn]
                assert len(fps) > 0, f'Did not find files for sample {s}'
                sampleFileDict[s] = fps

    tInit = time.time()
    logging.info(f'Samples to process: {samples}')
    for i, (s, fps) in enumerate(sampleFileDict.items()):
        ts = time.time()
        logging.info(f'Processing {i+1}/{len(samples)}: {s}')

        # prepare align arguments
        cmdList = ['bowtie2', '-x', genomeBowtie2Idx, '-p', args.ncpu]
        if args.isPe:
            assert len(peSfx) == 2
            samples1: list[str] = []
            samples2: list[str] = []
            for fp in fps:
                assert not (peSfx[0] in fp.stem and peSfx[1] in fp.stem)
                if peSfx[0] in fp.stem:
                    samples1.append(str(fp))
                elif peSfx[1] in fp.stem:
                    samples2.append(str(fp))
                else:
                    raise ValueError(f'File {fp} not bound to PE {peSfx}')
            cmdList.extend([
                '-1', ','.join(samples1),
                '-2', ','.join(samples2)
            ])
        else:
            cmdList.extend([
                '-U', ','.join([str(fp) for fp in fps])
            ])

        # prepare convert to bam arguments
        target = args.out / f'{s.strip("_")}.bam'

        if target.is_file():
            logging.info(f'Found existing bam file, skip')
            pass
        cmdList.extend([
            '|', 'samtools', 'view', '-bS', '-@', toBamNcpu,
            "|", "samtools", "sort", '-@', toBamNcpu, "--write-index", '-o', target
        ])
        logging.info(' '.join(cmdList))

        # Start running both
        cmd = withActivateEnvCmd(' '.join(cmdList), SHORTREADS_ENV)
        result = subprocess.run(cmd, shell=True,
                                capture_output=True, executable=SHELL)
        if result.returncode != 0:
            logging.info('stderr: ' + result.stderr.decode())
            logging.info('stdout: ' + result.stdout.decode())
        # stderr has logging.info info from bowtie2
        logging.info(result.stderr.decode())
        logging.info(f'Finished in {diffTime(ts)}\n')

    logging.info(f'All done, time elapsed {diffTime(tInit)}')
    logging.info('=' * 20 + getTime() + '=' * 20 + '\n' * 2)


def buildBowtie2idx(fs: list[Path], out: Path, name=None) -> Path:
    fs = [f.resolve() for f in fs]
    out = out.resolve()

    # check if idx file exists
    genomePath = fs[0].parent
    genomeName = fs[0].stem
    logging.info(f'genome name: {genomeName}')

    out.mkdir(exist_ok=True)
    bt2_base = ('_'.join(f.stem for f in fs) if name is None else name)
    outIdxForUse = out / bt2_base

    if any((
        outIdxForUse.with_suffix('.1.bt2').is_file(),
        outIdxForUse.with_suffix('.1.bt21').is_file()
    )):
        logging.info(
            f'Found index file {out.glob(bt2_base)}, will not make new ones.')
        return outIdxForUse

    # convert gbk to fa
    tempFiles: list[IO] = []
    if os.path.splitext(fs[0])[1] not in ['.fa', '.fasta', '.fna', '.fsa']:
        # try gbk
        newFps: list[Path] = []
        try:
            for f in fs:
                for s in SeqIO.parse(f, 'genbank'):
                    newF = NamedTemporaryFile()
                    SeqIO.write(s, newF.name, 'fasta')
                    newFps.append(Path(newF.name))
                    tempFiles.append(newF)  # for closing these files later
        except Exception as err:
            logging.error(f'Unexpected error {err=}, {type(err)=}')
            raise
        fs = newFps

    logging.info('-' * 20 + 'Indexing genome ' + getTime() + '-' * 20)
    args = [
        'bowtie2-build',
        ','.join(str(f) for f in fs),
        os.path.join(out, bt2_base),
    ]
    logging.info(' '.join(args))
    cmd = withActivateEnvCmd(' '.join(args), SHORTREADS_ENV)
    result = subprocess.run(args=args, capture_output=True, executable=SHELL)
    (f.close() for f in tempFiles)
    if result.returncode != 0:
        logging.info('stderr: ' + result.stderr.decode())
        logging.info('stdout: ' + result.stdout.decode())
        logging.info('-' * 20 + 'Error Indexing genome' + getTime() + '-' * 20)
        raise Exception
    logging.info('-' * 20 + 'DONE Indexing genome' + getTime() + '-' * 20)
    logging.info('\n' * 2)

    return outIdxForUse


def imputePeSuffix(
    rawFiles: list[Path],
    peSfx: list[str] = []
):
    if len(peSfx) == 0:
        return peSfx
    assert len(rawFiles) % 2 == 0 and len(rawFiles) != 0, (
        'Pair end reads should be in pairs, however'
        f' only {len(rawFiles)} found.')
    extensions = []
    fileNames = []
    for f in rawFiles:
        fn, suffix = splitStemSuffixIfCompressed(f, ['.gz'], fullSuffix=True)
        fileNames.append(fn)
        extensions.append(suffix)
    extensions = sorted(list(set(extensions)))
    assert len(extensions) == 1, f'Multiple file formats found {extensions}'
    peFound = False
    for i in range(1, 12):
        s: set[str] = set(fn[-i:] for fn in fileNames)
        if len(s) == 2:
            peSfx = sorted(s)
            peFound = True
        if len(s) != 2 and peFound:
            break
    if peFound:
        logging.info(f'Found pairend suffix {peSfx}')
        return peSfx
    raise ValueError(f'pair end suffix not found in {rawFiles}')


def getTime():
    return time.strftime('%z, %a, %d %b %Y, %H:%M:%S', time.localtime())


def diffTime(a):
    d = abs(time.time() - a)
    h = int(d // 3600)
    return str(h).zfill(2) + time.strftime(':%M:%S', time.gmtime(d))


if __name__ == "__main__":
    main()
