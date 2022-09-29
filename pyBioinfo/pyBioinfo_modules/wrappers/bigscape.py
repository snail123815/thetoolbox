import subprocess
from pathlib import Path
from typing import Literal
import numpy as np
from tempfile import NamedTemporaryFile

from pyBioinfo_modules.wrappers._environment_settings \
    import CONDAEXE, SHELL, BIGSCAPE_ENV, PFAM_DB, withActivateEnvCmd


def runBigscape(
    inputPath: Path,
    outputPath: Path,
    cpus: int = 4,
    cutoffs: list[float] = [0.2,],
    bigscapeEnv = BIGSCAPE_ENV,
    pfamDb = PFAM_DB,
    condaExe = CONDAEXE,
    shell = SHELL
) -> Path:
    if not outputPath.is_dir():
        outputPath.mkdir(parents=True, exist_ok=False)
    cmd = f'bigscape.py --mode auto -c {cpus}'
    cmd += f' --cutoffs {" ".join([str(c) for c in cutoffs])}'
    cmd += f' --pfam_dir {pfamDb}'
    cmd += f' -i {inputPath}'
    cmd += f' -o {outputPath}'
    bigscapeRun = subprocess.run(
        withActivateEnvCmd(cmd, bigscapeEnv, condaExe, shell),
        shell=True,  capture_output=True
        )
    assert bigscapeRun.returncode == 0, bigscapeRun.stderr.decode()
    return outputPath