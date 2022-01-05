import argparse
import os
from Bio import SeqIO
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(usage=('Generate a plot showing the distribution of contig lengths of all sequences in a fasta file'))
parser.add_argument('fasta', help='path to fasta file')

args = parser.parse_args()
fa = args.fasta

title = os.path.basename(fa)

seqs = list(SeqIO.parse(fa,'fasta'))
seqs.sort(key=lambda x: len(x), reverse=True)
seqIds = [s.id for s in seqs]
plt.plot([len(s) for s in seqs])
plt.title(title)
i = 0
while i < 15:
    try:
        t = f"{seqIds[i]}: {len(seqs[i]):.2e}"
        plt.text(0.98, 0.98-0.05*i, t,
                 transform=plt.gca().transAxes,
                 va='top', ha='right'
                )
    except:
        pass
    i += 1


plt.show()


