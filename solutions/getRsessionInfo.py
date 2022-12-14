import rpy2.robjects as robjects
import io
import os
import logging

def printToString(*args, **kwargs):
    with io.StringIO() as output:
        print(*args, file=output, **kwargs)
        contents = output.getvalue()
    return contents


def writeRSessionInfo(fileName, overwrite=True, logger=None):
    if isinstance(logger, type(None)):
        logger = logging.getLogger()
    sinfo = 'Attached packages in current R session:\n' + '=' * 80 + '\n'
    try:
        sinfo = 'Attached packages in current R session:\n' + '=' * 80 + '\n'
        for l in robjects.r('sessionInfo()["otherPkgs"]')[0]:
            x = printToString(l).split('\n')
            y = [a for a in x if len(a)>0 and not a.startswith('-- File:')]
            z = '\n'.join(y) + '\n\n' + "=" * 80 + '\n'
            sinfo += z
    except TypeError:
        logger.info('No R packages loaded in current session')
        return 0
    logger.info(sinfo)
    if not overwrite:
        if os.path.isfile(fileName):
            logger.info(f'File {fileName} exists, will not overwrite.')
            return 0
    with open(fileName, 'w') as fh:
        fh.write(sinfo)
    return 1
