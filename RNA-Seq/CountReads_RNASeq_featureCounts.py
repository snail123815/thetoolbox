from subprocess import run
from tempfile import NamedTemporaryFile
import os
import time
import argparse
import logging
from Bio import SeqIO
from BCBio import GFF


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='alignment folder (bam or sam files)')
    parser.add_argument('--output', required=True, help='gene counts outpuf folder')
    parser.add_argument('--gbk', required=True, help='genbank file with annotation')
    parser.add_argument('--ncpu', help='number of cpu to use')
    parser.add_argument('--isPe', action='store_true', help='set if pairend')
    parser.add_argument('-t', default='gene', help='target feature')
    parser.add_argument('-g', default='locus_tag', help='group factor')
    parser.add_argument('--fractionCounting', action='store_true', help='will add -M --fraction -O if True')
    parser.add_argument('--peLoose', action='store_true', help='will use loose configuration in pairend counting, set it if you donot want the features: -P, Check validity of paired-end distance when counting read pairs. Use -d and -D to set thresholds; -B, only count read pairs that have both ends aligned. (must set togeter with -P) ')


    args = parser.parse_args()
    alignFolder = args.input
    out = args.output
    gbk = args.gbk

    ncpu = args.ncpu
    isPe = args.isPe

    targetFeature = args.t
    groupFactor = args.g

    fractionCounting = args.fractionCounting
    peStrict = not args.peLoose

    logging.basicConfig(filename='!featureCounts.log', level=logging.DEBUG)

    # convert gbk to gff
    gffFile = NamedTemporaryFile('w+')
    GFF.write(SeqIO.parse(gbk,'genbank'), gffFile)
    gffFile.seek(0)
    gff = gffFile.name

    finalTs = time.time()
    if not os.path.isdir(out):
        os.makedirs(out)
    logging.info('='*20 + getTime() + '='*20)
    files = [f for f in os.listdir(alignFolder) if f.endswith('.bam')]
    if len(files) == 0:
        files = [f for f in os.listdir(alignFolder) if f.endswith('.sam')]
    assert len(files) != 0

    for i, f in enumerate(files):
        logging.info(f'Processing {i+1}/{len(files)}: ')
        ts = time.time()
        b = os.path.splitext(os.path.split(f)[-1])[0]
        args = [
            'featureCounts',
            '-T', ncpu,
            '-a', gff,
            '-t', targetFeature,
            '-g', groupFactor,
            '--minOverlap', '20', # Minimum number of overlapping bases in a read that is
                                  # required for read assignment. 1 by default. Number of
                                  # overlapping bases is counted from both reads if paired
                                  # end. If a negative value is provided, then a gap of up
                                  # to specified size will be allowed between read and the
                                  # feature that the read is assigned to.
            '--fracOverlap', '0.25', # Minimum fraction of overlapping bases in a read that is
                                     # required for read assignment. Value should be within range
                                     # [0,1]. 0 by default. Number of overlapping bases is
                                     # counted from both reads if paired end. Both this option
                                     # and '--minOverlap' option need to be satisfied for read
                                     # assignment.
            '-o', os.path.join(out, f'{b}.txt'),
            f'{alignFolder}/{f}'
        ]
        if fractionCounting:
            args.insert(3, '-M')  # Multi-mapping reads will also be counted. For a multi-
            # mapping read, all its reported alignments will be
            # counted. The 'NH' tag in BAM/SAM input is used to detect
            # multi-mapping reads.
            args.insert(4, '--fraction')  # Assign fractional counts to features. This option must
            # be used together with '-M' or '-O' or both. When '-M' is
            # specified, each reported alignment from a multi-mapping
            # read (identified via 'NH' tag) will carry a fractional
            # count of 1/x, instead of 1 (one), where x is the total
            # number of alignments reported for the same read. When '-O'
            # is specified, each overlapping feature will receive a
            # fractional count of 1/y, where y is the total number of
            # features overlapping with the read. When both '-M' and
            # '-O' are specified, each alignment will carry a fractional
            # count of 1/(x*y).
            args.insert(5, '-O')  # Assign reads to all their overlapping meta-features (or
            # features if -f is specified).

        if isPe:
            args.insert(3, '-p')  # If specified, fragments (or templates) will be counted
            # instead of reads. This option is only applicable for
            # paired-end reads; single-end reads are always counted as
            # reads.
            if peStrict:
                args.insert(4, '-P')  # Check validity of paired-end distance when counting read
                # pairs. Use -d and -D to set thresholds.
                args.insert(5, '-B')  # Only count read pairs that have both ends aligned. (must
                # set togeter with -P)
        logging.info(' '.join(args))
        res = run(args, capture_output=True)
        gffFile.close()
        if res.returncode != 0:
            logging.info(res.stdout.decode())
            logging.info(res.stderr.decode())
            raise Exception
        logging.info(res.stderr.decode())
        logging.info(f'Finished in {diffTime(ts)}\n')
    logging.info(f'All done, time elapsed {diffTime(finalTs)}')
    logging.info('='*20 + getTime() + '='*20)


def getTime():
    return time.strftime('%z, %a, %d %b %Y, %H:%M:%S', time.localtime())


def diffTime(a):
    d = abs(time.time()-a)
    h = int(d//3600)
    return str(h).zfill(2) + time.strftime(':%M:%S', time.gmtime(d))


if __name__ == "__main__":
    main()
