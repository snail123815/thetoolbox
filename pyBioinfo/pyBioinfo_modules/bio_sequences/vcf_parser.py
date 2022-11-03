from pathlib import Path
from typing import TypedDict, NotRequired, cast, Literal


class VarianceData(TypedDict):
    CHROM: str
    POS: int
    ID: str
    REF: str
    ALT: str
    QUAL: float
    FILTER: str
    INFO: str
    FORMAT: NotRequired[str]
    SAMPLES: NotRequired[dict[str, str]]

essentialVcfColumns = [
    'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'
]

def vcfParser(vcf: Path) -> list[VarianceData]:
    varianceDatas: list[VarianceData] = []
    with open(vcf, 'r') as vcfFh:
        getFormat: bool
        # Parse header lines
        for line in vcfFh:
            if line.startswith('##'):  # header lines
                continue
            elif line.startswith('#'):  # column names
                columns = line.strip()[1:].split('\t')
                break
        # Assert correct format
        assert 'columns' in locals(), f'Header line not found for {vcf}'
        assert len(columns) >= 8, \
            f'Header is not complete for {vcf}, {len(columns)} columns.'
        assert columns[:8] == essentialVcfColumns, \
            f'The file {vcf} is not a standard VCF file, \n{columns}'
        getFormat = len(columns) > 8 and columns[8] == 'FORMAT'
        if getFormat:
            samplesColumns = columns[9:]
            assert len(samplesColumns) > 0, (f'No sample is found in {vcf}'
                                             ' but FORMAT field present.')
        # Parse data
        for line in vcfFh:
            fields = line.strip().split('\t')
            varianceData = VarianceData(
                CHROM=fields[0],
                POS=int(fields[1]),
                ID=fields[2],
                REF=fields[3],
                ALT=fields[4],
                QUAL=(0 if fields[5] == '.' else float(fields[5])),
                FILTER=fields[6],
                INFO=fields[7],
            )
            if getFormat:
                varianceData['FORMAT'] = fields[8]
                varianceData['SAMPLES'] = {
                    s: fields[9+i] for i, s in enumerate(samplesColumns)
                }
            varianceDatas.append(varianceData)
    return varianceDatas
