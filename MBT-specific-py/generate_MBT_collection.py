from pathlib import Path
from typing import TypedDict
from typing import NewType
from typing import cast
from pyBioinfo_modules.wrappers.prokka import runProkka
from pyBioinfo_modules.wrappers.antismash import runAntismash
from pyBioinfo_modules.basic.decompress import decompFileIfCompressed
from pyBioinfo_modules.bio_sequences.protein_from_gbk import getFaaFromGbk
from pyBioinfo_modules.bio_sequences.gbk_to_gff import gbkToGff
from tqdm import tqdm
from multiprocessing import Pool
import shutil

rootPath = (Path.home() / 'gdata/MBT-collection').resolve()

targetDirs = [
    rootPath / 'MBT-initial-96strains',
    rootPath / 'MBT-continue-collection',
    rootPath / 'MBT-noRawData-collection'
]
collectionGbkDir = (rootPath / 'collective-gbk').resolve()
collectionGffDir = (rootPath / 'collective-gff').resolve()
collectionFnaDir = (rootPath / 'collective-fna').resolve()
collectionFaaDir = (rootPath / 'collective-faa').resolve()
collectionAntismashDir = (rootPath / 'collective-antismash').resolve()
knownSpecies = [
    'Streptomyces',
    'Kitasatospora',
    'Paenibacillus',
    'Terrbacter',
    'Actinobacteria',
    'Streptosporangiales'
]
fnaExtensions = ['.fna', '.fasta', '.fa']
gbkExtensions = ['.gbk', '.gb', '.genbank', '.gbff']
gffExtensions = ['.gff']
faaExtensions = ['.faa']
highPriorityKeywords = ['RefSeq']
priorityExceptions = ['MBT56', 'MBT84']
defaultStrainName = 'Actinobacteria_bacterium_'
annotationPriority = ['ncbi', 'pgap', 'novogene', 'prokka']

allowedCompressions = ['.gz']
annotationUnderStrainRoot = Path('annotation')
assemblyUnderStrainRoot = Path('assembly')
antismashUnderStrainRoot = Path('antismash')

runAntismashParallel = 4


class StrainPaths(TypedDict):
    rootPath: Path
    fnaPath: Path
    gbkPath: Path | None
    gffPath: Path | None
    faaPath: Path | None
    antismashPath: Path | None


StrainsInfo = NewType('StrainsInfo', dict[str, StrainPaths])


def getSuffixWithCompression(p: Path) -> str:
    return (p.suffixes[-2] if p.suffix.lower()
            in allowedCompressions else p.suffix.lower())


def findAnnotation(
    strainRoot: Path,
    removePartial: bool = True
) -> tuple[Path, Path, Path] | tuple[None, None, None]:
    annotationRoot = strainRoot / annotationUnderStrainRoot
    try:
        annotationDirs = [d for d in annotationRoot.iterdir() if d.is_dir()]
        sourceDirs = []
        for source in annotationPriority:
            sourceDirs.extend([d for d in annotationDirs
                               if d.name.startswith(source)])
        sourceDirs.sort(key=lambda d: d.name, reverse=True)
        gbkfaaPairs = []
        for d in sourceDirs:
            gbks = [f for f in d.iterdir()
                    if getSuffixWithCompression(f) in gbkExtensions]
            for faaExt in faaExtensions:
                faas = [
                    gbk.with_suffix('').with_suffix(faaExt)
                    if gbk.suffix in allowedCompressions
                    else gbk.with_suffix(faaExt)
                    for gbk in gbks
                ]
            for gffExt in gffExtensions:
                gffs = [
                    gbk.with_suffix('').with_suffix(gffExt)
                    if gbk.suffix in allowedCompressions
                    else gbk.with_suffix(gffExt)
                    for gbk in gbks
                ]
            if len(gbks) == 0:
                if removePartial:
                    shutil.rmtree(str(d.resolve()))
            for gbk, gff, faa in zip(gbks, gffs, faas):
                if not faa.exists():
                    faaGenerated = getFaaFromGbk(gbk)
                    assert faaGenerated.exists()
                    faa = faaGenerated
                if not gff.exists():
                    gffGenerated = gbkToGff(gbk)
                    assert gffGenerated.exists()
                    gff = gffGenerated
                gbkfaaPairs.append((gbk, gff, faa))

        gbkPath, gffPath, faaPath = gbkfaaPairs[0]
    except (IndexError, FileNotFoundError):
        return None, None, None
    return gbkPath, gffPath, faaPath


def findAssembly(strainRoot: Path) -> Path | None:
    assemblyPath = strainRoot / assemblyUnderStrainRoot
    try:
        fnaFiles = [f for f in assemblyPath.iterdir()
                    if getSuffixWithCompression(f) in fnaExtensions]
        fnaPath = fnaFiles[0]
        if len(fnaFiles) > 1:
            print(f'Many assemblies found for strain {strainRoot}:' +
                  '\n'.join(f.name for f in fnaFiles) +
                  '\nThe first one will be used.')
    except (FileNotFoundError, IndexError):
        print(f'No valid assembly data found in dir {strainRoot}')
        return None
    return fnaPath.resolve()


def findAntismash(strainRoot: Path, removePartial: bool = True) -> Path | None:
    try:
        _antismashPaths = [p for p in strainRoot.iterdir()
                           if p.is_dir() and
                           p.name.startswith(antismashUnderStrainRoot.name)]
        antismashPaths: list[Path] = []
        for p in _antismashPaths:
            if 'index.html' not in [
                f.name for f in p.iterdir() if f.is_file()
            ]:
                if removePartial:
                    shutil.rmtree(str(p))
            else:
                antismashPaths.append(p)
        antismashPaths.sort(key=lambda p: p.name, reverse=True)
        return antismashPaths[0].resolve()
    except (FileNotFoundError, IndexError):
        return None


def generateStrainName(strain: str) -> str:
    nameSplit = strain.split('_')
    if nameSplit[0] in knownSpecies:
        strainName = "_".join(nameSplit)
    else:
        strainName = defaultStrainName + '_'.join(nameSplit)
    return strainName


strainsInfo = StrainsInfo({})
for root in [Path(d) for d in targetDirs]:
    singleStrainDirs = [d for d in root.iterdir() if d.is_dir()]
    for d in singleStrainDirs:
        gbkPath, gffPath, faaPath = findAnnotation(d)
        fnaPath = findAssembly(d)
        if fnaPath is None:
            continue
        antismashPath = findAntismash(d)
        strainName = generateStrainName(d.name)

        assert strainName not in strainsInfo, \
            f'Duplicated strain name {strainName} in dir {root}, ' \
            + f'duplicated item is {strainsInfo[strainName]}'

        strainsInfo[strainName] = {
            'rootPath': d.resolve(),
            'antismashPath': antismashPath,
            'gbkPath': gbkPath,
            'gffPath': gffPath,
            'faaPath': faaPath,
            'fnaPath': fnaPath
        }

strainNamesAll = list(strainsInfo.keys())
for strainName in strainNamesAll:
    for kw in highPriorityKeywords:
        if kw in strainName:
            # find the other strainName with the same prefix and remove that
            prefix = strainName.split(kw)[0]
            isException = any(exp in prefix for exp in priorityExceptions)
            if not isException:
                for sn in strainNamesAll:
                    if sn not in strainsInfo:
                        continue
                    if sn == strainName:
                        continue
                    if sn.startswith(prefix):
                        print(f'{sn} removed because it is a duplicate of '
                              + f'{strainName}')
                        del (strainsInfo[sn])

noAssemblyStrains = [strain for strain in strainsInfo if
                     strainsInfo[strain]['fnaPath'] is None]
assert len(noAssemblyStrains) == 0
noAnnotationStrains = [strain for strain in strainsInfo if
                       strainsInfo[strain]['gbkPath'] is None]
noAntismashStrains = [strain for strain in strainsInfo if
                      strainsInfo[strain]['antismashPath'] is None]
# assert all(s in noAntismashStrains for s in noAnnotationStrains)

print('No annotation strains:')
print(noAnnotationStrains)
print('No antismash strains:')
print(noAntismashStrains)

if len(noAnnotationStrains) > 0:
    with tqdm(noAnnotationStrains) as tqdmBar:
        for strain in tqdmBar:
            tqdmBar.set_description(f'Prokka for {strain}')
            strainNameSplits = strain.split('_')
            prokkaDir = runProkka(
                fastaPath=strainsInfo[strain]['fnaPath'],
                gcode=11,
                gram='pos',
                center='MBT',
                genus=strainNameSplits[0],
                species=strainNameSplits[1],
                strain='_'.join(strainNameSplits[2:]),
                locustag='_'.join(strainNameSplits[2:]),
                prokkaEnv=Path('~/genvs/quasan'),
                shell='zsh',
                cpu=16,
                output=strainsInfo[strain]['rootPath']
                / annotationUnderStrainRoot,
                silent=True
            )
            gbkPath, gffPath, faaPath = findAnnotation(
                strainsInfo[strain]['rootPath'])
            assert gbkPath is not None, \
                f'Annotation not found in {prokkaDir}'
            strainsInfo[strain]['gbkPath'] = gbkPath
            strainsInfo[strain]['gffPath'] = gffPath
            strainsInfo[strain]['faaPath'] = faaPath

if len(noAntismashStrains) > 0:
    runnerPool = Pool(runAntismashParallel)
    results = []
    for strain in noAntismashStrains:
        gbkPath = strainsInfo[strain]['gbkPath']
        assert gbkPath is not None, "You need to annotate" + \
            f" {strainsInfo[strain]['fnaPath']} first."
        kwargs = {
            "genbankFilePath": gbkPath,
            "title": strain,
            "taxon": 'bacteria',
            "completeness": 2,
            "condaEnv": Path('~/genvs/quasan'),
            "cpu": 16,
            "shell": 'zsh',
            "output": strainsInfo[strain]['rootPath'],
            "silent": True
        }
        results.append(
            (strain,
             runnerPool.apply_async(runAntismash, kwds=kwargs))
        )
    runnerPool.close()
    with tqdm(results) as tqdmBar:
        for strain, res in tqdmBar:
            strainsInfo[strain]['antismashPath'] = res.get()
            tqdmBar.set_description(f'Running antismash for {strain}')


noAnnotationStrains = [strain for strain in strainsInfo if
                       strainsInfo[strain]['gbkPath'] is None]
noAntismashStrains = [strain for strain in strainsInfo if
                      strainsInfo[strain]['antismashPath'] is None]
assert len(noAnnotationStrains) + len(noAnnotationStrains) == 0

for p in [collectionGbkDir, collectionFnaDir, collectionGffDir,
          collectionFaaDir, collectionAntismashDir]:
    p.mkdir(exist_ok=True)
for strain, strainPaths in strainsInfo.items():
    toGbk = collectionGbkDir / f'{strain}.gbk'
    toGff = collectionGffDir / f'{strain}.gff'
    toFaa = collectionFaaDir / f'{strain}.faa'
    toFna = collectionFnaDir / f'{strain}.fna'
    if not toGbk.exists():
        gbk, rmgbk = decompFileIfCompressed(cast(Path, strainPaths['gbkPath']))
        if rmgbk:
            shutil.move(gbk, toGbk)
        else:
            shutil.copyfile(gbk, toGbk)
    if not toGff.exists():
        gff, rmgff = decompFileIfCompressed(cast(Path, strainPaths['gffPath']))
        if rmgff:
            shutil.move(gff, toGff)
        else:
            shutil.copyfile(gff, toGff)
    if not toFaa.exists():
        faa, rmfaa = decompFileIfCompressed(cast(Path, strainPaths['faaPath']))
        if rmfaa:
            shutil.move(faa, toFaa)
        else:
            shutil.copyfile(faa, toFaa)
    if not toFna.exists():
        fna, rmfna = decompFileIfCompressed(strainPaths['fnaPath'])
        if rmfna:
            shutil.move(fna, toFna)
        else:
            shutil.copyfile(fna, toFna)
    antismashRoot = collectionAntismashDir / strain
    bgcGbks = list(cast(Path, strainPaths['antismashPath']).glob(
        '*region[0-9][0-9][0-9].gbk'))

    if len(bgcGbks) > 0:
        antismashRoot.mkdir(exist_ok=True)
        for f in bgcGbks:
            toBgcGbk = antismashRoot / f'{strain}_{f.name}'
            if not toBgcGbk.exists():
                shutil.copyfile(f, toBgcGbk)
