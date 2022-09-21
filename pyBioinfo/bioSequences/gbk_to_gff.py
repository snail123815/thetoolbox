from pathlib import Path
from BCBio import GFF
from Bio import SeqIO
from ..wrappers.decompress \
    import decompFileIfCompressed, getSuffixAfterDecompression
from .bio_seq_file_extensions import GBK_EXTENSIONS

def gbkToGff(path: Path) -> Path:
    assert getSuffixAfterDecompression(path).lower() in GBK_EXTENSIONS
    gbkPath, didUnzip = decompFileIfCompressed(path)
    gffPath = gbkPath.with_suffix('.gff')
    with gffPath.open('w') as gff:
        GFF.write(SeqIO.parse(gbkPath, 'genbank'), gff)
    assert gffPath.is_file()
    if didUnzip:
        gbkPath.unlink()
    return gffPath