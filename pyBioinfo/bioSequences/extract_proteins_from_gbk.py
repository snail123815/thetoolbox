import argparse
import subprocess
import os
from pathlib import Path


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError


def getProteins(seqObj, codonTable=11):
    """Extract proteins from a SeqRecord.""" \
        + """ If your input file have multiple contigs, do a loop"""
    proteins = []
    for feat in seqObj.features:
        if feat.type == "CDS":
            proteinLocustag = feat.qualifiers['locus_tag'][0]
            try:
                proteinGeneId = feat.qualifiers['gene'][0]
                proteinGeneId = proteinGeneId[0].upper() + proteinGeneId[1:]
            except KeyError:
                proteinGeneId = ''
            try:
                proteinTranslation = Seq(feat.qualifiers['translation'][0])
            except KeyError:
                proteinGene = seqObj.seq[feat.location.start:feat.location.end]
                if len(proteinGene) % 3 == 0:
                    proteinGene = (proteinGene if feat.location.strand == 1
                                   else proteinGene.reverse_complement())
                    try:
                        # Get a reliable tranlsation
                        proteinTranslation = proteinGene.translate(
                            to_stop=True,
                            cds=True,
                            table=codonTable
                        )
                    except TranslationError:
                        continue
                else:
                    continue
            proteinProduct = feat.qualifiers['product'][0]
            p = SeqRecord(proteinTranslation,
                          id=proteinLocustag,
                          name=proteinGeneId,
                          description=proteinProduct)
            proteins.append(p)
    return proteins


def getFaaFromGbk(gbkPath: Path, codonTable=11) -> Path:
    unzip = False
    if gbkPath.suffix == '.gz':
        gzipd = subprocess.run(
            f'gzip -dkf {gbkPath}'.split(),
            capture_output=True
        )
        assert gzipd.returncode == 0
        gbkPath = gbkPath.with_suffix('')
        assert gbkPath.exists(), '\n'.join([
            f'Unzip file {gbkPath} failed with error message:',
            gzipd.stderr.decode(), gzipd.stdout.decode()
        ])
        unzip = True
    try:
        faaPath = gbkPath.with_suffix('.faa')
        proteins = []
        for s in SeqIO.parse(str(gbkPath), 'genbank'):
            proteins.extend(getProteins(s, codonTable=codonTable))
        n = SeqIO.write(proteins, faaPath, 'fasta')
        print(f'Successfully wrote {n} proteins')
    except Exception as e:
        raise e
    finally:
        if unzip:
            os.remove(str(gbkPath))
    return faaPath


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('file', help='genbank file')

    args = argparser.parse_args()
    gbkPath = Path(args.file)
    faaPath = getFaaFromGbk(gbkPath)
    print(faaPath)
