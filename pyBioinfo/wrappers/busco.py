# Note that busco cannot run multiple instances from the same
# executable.
import subprocess
from pathlib import Path
from typing import Literal
from _environment_settings import \
    SHELL, CONDAEXE, BUSCO_ENV, getActivateEnvCmd


def runBusco(
    targetProteome: Path,
    outName: str,
    outPath: Path,
    condaEnv: Path | None = Path(BUSCO_ENV),
    condaExe: Literal['conda', 'micromamba', 'mamba'] = CONDAEXE,
    shell: Literal['bash', 'zsh'] = SHELL,
    silent: bool = False,
    cpu: int = 4
) -> Path:

    cmd = (
        f'busco --auto-lineage-prok -m prot -f -c {cpu}'
        + f' -i {targetProteome}'
        + f' --out_path {outPath}'
        + f' -o {outName}'
    )
    if not silent:
        print(cmd)

    cmd = ' && '.join([
        getActivateEnvCmd(condaEnv, condaExe, shell),
        cmd
    ])

    commandResult = subprocess.run(
        cmd, capture_output=True, shell=True,
        executable=shell
    )
    if commandResult.returncode != 0:
        print('Failed busco:')
        print(cmd)
        print(commandResult.stdout.decode())
        print(commandResult.stderr.decode())

    return outPath.resolve()
