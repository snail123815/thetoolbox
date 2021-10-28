import rpy2.robjects as robjects
import io
import logging

def printToString(*args, **kwargs):
    with io.StringIO() as output:
        print(*args, file=output, **kwargs)
        contents = output.getvalue()
    return contents


def writeRSession(fileName, logger=None):
    if isinstance(logger, type(None)):
        logger = logging.getLogger()
    with open(fileName, 'w') as fh:
        fh.write('Attached packages in current R session:\n' + '=' * 80 + '\n')
        sinfo = ''
        for l in robjects.r('sessionInfo()["otherPkgs"]')[0]:
            x = printToString(l).split('\n')
            y = [a for a in x if len(a)>0 and not a.startswith('-- File:')]
            z = '\n'.join(y) + '\n\n' + "=" * 80 + '\n'
            sinfo += z
        fh.write(sinfo)
        logger.info(sinfo)
        
