from pathlib import Path
import argparse
from Bio import SeqIO


from pyBioinfo_modules.bio_sequences.vcf_parser import (
    vcfParser, applyVariancesOnSeqRecords
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', type=Path)
    parser.add_argument('--genome', type=Path, help='genbank format')
    parser.add_argument('--output', type=Path)

    args = parser.parse_args()

    varianceDatas = vcfParser(args.vcf)

    seqRecordDict = {
        rec.id: rec for rec in SeqIO.parse(args.genome, 'genbank')
    }

    variantSeqRecordDict = applyVariancesOnSeqRecords(
        varianceDatas, seqRecordDict
    )

    SeqIO.write(variantSeqRecordDict.values(), args.output)


if __name__ == '__main__':
    main()
