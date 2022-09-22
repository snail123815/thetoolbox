import subprocess
import os
from pathlib import Path


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

from pyBioinfo_modules.basic.decompress import decompFileIfCompressed
from pyBioinfo_modules.bioSequences.bio_seq_file_extensions \
    import GBK_EXTENSIONS

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
    gbkPath, unzip = decompFileIfCompressed(gbkPath)
    faaPath = gbkPath.with_suffix('.faa')
    try:
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
