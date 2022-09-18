
# TODO: write a script change this file
SHELL = 'zsh'
CONDAEXE = 'micromamba'

ANTISMASH_ENV = '~/genvs/quasan'
BUSCO_ENV = '~/genvs/quasan'
PROKKA_ENV = '~/genvs/quasan'

def getActivateEnvCmd(condaEnv, condaExe=CONDAEXE, shell=SHELL):
    if condaExe == 'micromamba':
        activateEnvCmd = f'eval "$(micromamba shell hook --shell={shell})"'
    else:
        activateEnvCmd = f'eval "$({condaExe} shell.{shell} hook)"'
    activateEnvCmd += f' && micromamba activate {condaEnv}'
    return activateEnvCmd