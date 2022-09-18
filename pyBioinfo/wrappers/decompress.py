from pathlib import Path
import subprocess
from tempfile import NamedTemporaryFile
from tempfile import _TemporaryFileWrapper

IMPLEMENTED_COMPRESSION_FORMATS: list[str] = ['.gz', '.xz']


def decompressFile(filePath: Path) -> Path:
    match filePath.suffix:
        case '.gz':
            decompress = subprocess.run(
                f'gzip -dkf {filePath.resolve()}'.split(' '),
                capture_output=True
            )
            resultFilePath = filePath.with_suffix('')
        case '.xz':
            decompress = subprocess.run(
                f'xz -dkf {filePath.resolve()}'.split(' '),
                capture_output=True
            )
            resultFilePath = filePath.with_suffix('')
        case _:
            raise NotImplementedError(
                f'Decompress {filePath.suffix} file is not supported.')
    assert decompress.returncode == 0
    assert resultFilePath.exists(), '\n'.join([
        f'Unzip file {filePath} failed with error message:',
        decompress.stderr.decode(), decompress.stdout.decode()
    ])
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
    outTempFIle = NamedTemporaryFile()
    with open(outTempFIle.name, 'w') as out:
        unzip = subprocess.run([prog, '-dkc', filePath], stdout=out)
    if unzip.returncode != 0:
        print(f'unzip "{filePath}" failed')
        print(' '.join(unzip.args))
        print(unzip.stderr.decode())
        print(unzip.stdout.decode())
    return outTempFIle
