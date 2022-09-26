import subprocess
from pathlib import Path
from typing import Literal

from pyBioinfo_modules.wrappers._environment_settings \
    import MASH_ENV, getActivateEnvCmd


def mashSketchFiles(
    inputFiles: list[Path],
    output: Path,
    kmer: int,
    sketch: int,
    ncpu: int = 1,
    molecule: Literal['DNA', 'protein'] = "DNA"
) -> Path:
    """
    Calculates the distance between the query fasta files
    stored in the sketch file by using mash.
    """
    cmd = f"mash sketch -o {output} -k {kmer} -p {ncpu} -s {sketch}"
    cmd += (' -a' if molecule == 'protein' else '')
    cmd += ' '.join([str(f) for f in inputFiles])
    mashSketchRun = subprocess.run(
        getActivateEnvCmd(MASH_ENV) + cmd,
        shell=True, capture_output=True, check=True)
    outputMsh = Path(str(output) + '.msh')
    assert outputMsh.is_file()
    return outputMsh


def mashDistance(
    inputMsh: Path,
    outputFile: Path,
) -> Path:
    """
    Calculates the distance between the query fasta files
    stored in the sketch file by using mash.
    Parameters
    ----------
    outdir
        string, the path to output directory
    returns
    ----------
    """
    cmd = f"mash dist {inputMsh} {inputMsh} > {outputFile}"
    mashDistRun = subprocess.run(
        getActivateEnvCmd(MASH_ENV) + cmd,
        shell=True, check=True
        )
    assert outputFile.exists
    return outputFile
