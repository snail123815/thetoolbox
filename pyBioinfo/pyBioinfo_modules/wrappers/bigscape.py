from shutil import which
import subprocess
from pathlib import Path
from typing import Literal
import numpy as np
from tempfile import NamedTemporaryFile

from pyBioinfo_modules.wrappers._environment_settings \
    import CONDAEXE, SHELL, BIGSCAPE_ENV, PFAM_DB, withActivateEnvCmd


def whichBigscape(shell=SHELL):
    try:
        bigscapeExe = 'bigscape.py'
        subprocess.run(withActivateEnvCmd('bigscape.py --version',BIGSCAPE_ENV,CONDAEXE, SHELL), shell=True,
                       check=True, executable=shell)
    except subprocess.CalledProcessError:
        bigscapeExe = 'bigscape'
    subprocess.run(withActivateEnvCmd(f'{bigscapeExe} --version',BIGSCAPE_ENV,CONDAEXE, SHELL), shell=True,
                   check=True, executable=shell)
    return bigscapeExe


bigscapeExe = whichBigscape(shell=SHELL)


def runBigscape(
    inputPath: Path,
    outputPath: Path,
    cpus: int = 4,
    cutoffs: list[float] = [0.2, ],
    bigscapeEnv=BIGSCAPE_ENV,
    pfamDb=PFAM_DB,
    condaExe=CONDAEXE,
    shell=SHELL
) -> Path:
    if not outputPath.is_dir():
        outputPath.mkdir(parents=True, exist_ok=False)
    cmd = f'{bigscapeExe} --mode auto -c {cpus}'
    cmd += f' --cutoffs {" ".join([str(c) for c in cutoffs])}'
    cmd += f' --pfam_dir {pfamDb}'
    cmd += f' -i {inputPath}'
    cmd += f' -o {outputPath}'
    # print(withActivateEnvCmd(cmd, bigscapeEnv, condaExe, shell))
    bigscapeRun = subprocess.run(
        withActivateEnvCmd(cmd, bigscapeEnv, condaExe, shell),
        shell=True, capture_output=True, executable=shell
    )
    # print(bigscapeRun.returncode)
    # print(bigscapeRun.stdout.decode())
    # print(bigscapeRun.stderr.decode())
    assert bigscapeRun.returncode == 0, bigscapeRun.stderr.decode()
    return outputPath
