
# TODO: write a script change this file
from typing import Literal
from pathlib import Path


SHELL: Literal['bash', 'zsh'] = 'zsh'
CONDAEXE: Literal['conda', 'mamba', 'micromamba'] = 'micromamba'

ANTISMASH_ENV: Path | None = Path.home()/'genvs/quasan'
BUSCO_ENV: Path | None = Path.home()/'genvs/quasan'
PROKKA_ENV: Path | None = Path.home()/'genvs/quasan'
MASH_ENV: Path | None = Path.home()/'genvs/phylophlan'


def getActivateEnvCmd(condaEnv, condaExe=CONDAEXE, shell=SHELL):
    if condaExe == 'micromamba':
        activateEnvCmd = f'eval "$(micromamba shell hook --shell={shell})"'
    else:
        activateEnvCmd = f'eval "$({condaExe} shell.{shell} hook)"'
    activateEnvCmd += f' && micromamba activate {condaEnv}'
    return activateEnvCmd
