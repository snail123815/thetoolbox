# For RNA-Seq alignment, with each file or pair of files per sample. Out
# put aligned file per sample.


from subprocess import run
from tempfile import NamedTemporaryFile
from Bio import SeqIO
import os
import time
import argparse
import logging
from pathlib import Path
from pyBioinfo_modules.basic.decompress import splitStemSuffixIfCompressed

env = '~/genvs/shortReads/'


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
        help='suffix to pairend file name. eg. a_1_.fq.gz and a_2_.fq.gz, then set this to --pesuffix _1_ _2_')
    parser.add_argument(
        '--ncpu',
        default=1,
        type=int,
        help='number of cpu to use')
    parser.add_argument(
        '--sampleNames',
        nargs="+",
        help='sample names if file names contain random generated string which will unpair the samples after removing pairend suffix')

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
    files = [f for f in args.raw.iterdir() if f.is_file()]
    [files.extend(f for f in d.iterdir() if f.is_file())
     for d in args.raw.iterdir() if d.is_dir()]
    for f in files.copy():
        if f.suffix == '.gz':
            ext = f.with_suffix('').suffix
        else:
            ext = f.suffix
        if ext not in ['.fastq', '.fq']:
            files.pop(f)
    assert len(files) > 0, f'Files not found in {args.raw}'
    files.sort()

    # peSfx defined here if not defined
    if isPe:
        peSfx = imputePeSuffix(files, rawExt, peSfx)

    # samples defined here if not defined
    if isinstance(samples, type(None)):
        fns = [os.path.split(f)[1][:-len(rawExt)] for f in files]
        if isPe:
            samples = list(set(f[:-len(peSfx[0])] for f in fns))
        else:
            samples = fns
        samples.sort()

    sampleFileDict = {}
    for s in samples:
        fs = sorted([f for f in files if s in f])
        if isPe:
            sampleFileDict[s] = {peSfx[0]: [f for f in fs if f[:-len(rawExt)].endswith(peSfx[0])],
                                 peSfx[1]: [f for f in fs if f[:-len(rawExt)].endswith(peSfx[1])]}
        else:
            sampleFileDict[s] = fs

    totalTs = time.time()
    assert len(samples) > 0
    logging.info(f'Samples to process: {samples}')
    for i, s in enumerate(samples):
        ts = time.time()
        logging.info(f'Processing {i+1}/{len(samples)}: {s}')

        # prepare align arguments
        args = ['bowtie2', '-x', genomeBowtie2Idx, '-p', ncpu]
        if isPe:
            args.extend([
                '-1', ','.join(sampleFileDict[s][peSfx[0]]),
                '-2', ','.join(sampleFileDict[s][peSfx[1]])
            ])
        else:
            args.extend([
                '-U', ','.join(sampleFileDict[s])
            ])

        # prepare convert to bam arguments
        target = os.path.join(outFolder, f'{s.strip("_")}.bam')

        if os.path.isfile(target):
            pass
        args.extend([
            '|', 'samtools', 'view', '-bS', '-@', toBamNcpu,
            "|", "samtools", "sort", '-@', toBamNcpu, "--write-index", '-o', target
        ])
        logging.info(' '.join(args))

        # Start running both
        result = run(args=' '.join(args), shell=True, capture_output=True)
        if result.returncode != 0:
            logging.info('stderr: ' + result.stderr.decode())
            logging.info('stdout: ' + result.stdout.decode())
        # stderr has logging.info info from bowtie2
        logging.info(result.stderr.decode())
        logging.info(f'Finished in {diffTime(ts)}\n')

    logging.info(f'All done, time elapsed {diffTime(totalTs)}')
    logging.info('=' * 20 + getTime() + '=' * 20 + '\n' * 2)


def buildBowtie2idx(fs: list[Path], out: Path, name=None):
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
    converted = []
    if os.path.splitext(fs[0])[1] not in ['.fa', '.fasta', '.fna', '.fsa']:
        # try gbk
        newFs = []
        try:
            for f in fs:
                for s in SeqIO.parse(f, 'genbank'):
                    newF = NamedTemporaryFile()
                    SeqIO.write(s, newF.name, 'fasta')
                    newFs.append(newF.name)
                    converted.append(newF)  # for closing these files later
        except BaseException as err:
            # python >= 3.8
            logging.error(f'Unexpected error {err=}, {type(err)=}')
        fs = newFs

    logging.info('-' * 20 + 'Indexing genome ' + getTime() + '-' * 20)
    args = [
        'bowtie2-build',
        ','.join(fs),
        os.path.join(out, bt2_base),
    ]
    logging.info(' '.join(args))
    result = run(args=args, capture_output=True)
    for f in converted:
        f.close()
    if result.returncode != 0:
        logging.info('stderr: ' + result.stderr.decode())
        logging.info('stdout: ' + result.stdout.decode())
        logging.info('-' * 20 + 'Error Indexing genome' + getTime() + '-' * 20)
        raise Exception
    logging.info('-' * 20 + 'DONE Indexing genome' + getTime() + '-' * 20)
    logging.info('\n' * 2)

    return os.path.join(out, bt2_base)


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
