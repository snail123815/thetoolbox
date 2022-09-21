from pathlib import Path
import subprocess
from tempfile import NamedTemporaryFile
from tempfile import _TemporaryFileWrapper

IMPLEMENTED_COMPRESSION_FORMATS: list[str] = ['.gz', '.xz']


def getSuffixAfterDecompression(filePath: Path) -> str:
    if filePath.suffix in IMPLEMENTED_COMPRESSION_FORMATS:
        return filePath.with_suffix('').suffix
    else:
        return filePath.suffix


def decompressFile(filePath: Path) -> Path:
    match filePath.suffix:
        case '.gz':
            prog = 'gzip'
        case '.xz':
            prog = 'xz'
        case _:
            raise NotImplementedError(
                f'Decompress {filePath.suffix} file is not supported.')
    decompress = subprocess.run(
        f'{prog} -dkf {filePath.resolve()}'.split(' '),
        capture_output=True, check=True
    )
    resultFilePath = filePath.with_suffix('')
    if not resultFilePath.exists():
        raise FileNotFoundError('\n'.join([
            f'Unzip file {filePath} succeed but no output file found:',
            ' '.join(decompress.args),
            decompress.stderr.decode(), decompress.stdout.decode()
        ]))
    return resultFilePath


def decompFileIfCompressed(filePath: Path) -> tuple[Path, bool]:
    if filePath.suffix in IMPLEMENTED_COMPRESSION_FORMATS:
        return decompressFile(filePath), True
    else:
        return filePath, False


def decompressToTempTxt(filePath: Path) -> _TemporaryFileWrapper:
    # To access the file outside this function,
    # do NOT: a = decompressToTempTxt(filePath).name
    # it will kill the temporary file.
    # You can do with decompressToTempTxt(filePath) as t:
    match filePath.suffix:
        case '.gz':
            prog = 'gzip'
        case '.xz':
            prog = 'xz'
        case _:
            raise NotImplementedError(
                f'Decompress {filePath.suffix} file is not supported.')
    outTempFile = NamedTemporaryFile()
    with open(outTempFile.name, 'w') as out:
        subprocess.run([prog, '-dkc', filePath], stdout=out, check=True)
    return outTempFile
