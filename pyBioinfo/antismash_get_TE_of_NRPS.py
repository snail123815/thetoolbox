import argparse
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pathlib import Path
from pyBioinfo_modules.wrappers.antismash import find_NRPS_TE_domain

parser = argparse.ArgumentParser()
parser.add_argument('-p', type=Path, help='Root path of all AntiSMASH folders')
args = parser.parse_args()

rootAllAntismashes = args.p

allAntismashes = [d for d in rootAllAntismashes.iterdir() if d.is_dir()]
asResultJsons = []
for d in allAntismashes:
    jsons = list(d.glob('*.json'))
    assert len(jsons) == 1, (
        f'AntiSMASH result in "{d.name}" has some problem:'
        f'    Found {len(jsons)} of ".json" files, but there should be one.'
    ) 
    asResultJsons.append(jsons[0])

teDomains: list[SeqRecord] = []

for asResultJson in asResultJsons:
    teDomains.extend(find_NRPS_TE_domain(asResultJson))

SeqIO.write(
    teDomains,
    (rootAllAntismashes / 'NRPS_TE_domains.fasta').open('w'),
    'fasta'
)