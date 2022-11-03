from pathlib import Path
from typing import TypedDict, NotRequired, cast, Literal
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, ExactPosition


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
                    s: fields[9 + i] for i, s in enumerate(samplesColumns)
                }
            varianceDatas.append(varianceData)
    return varianceDatas


def filterVarianceData(
    varianceDatas: list[VarianceData],
    locationRange: tuple[int, int],
    quality: float,
    filter: bool,  # on FILTER column (PASS or not)
    # on INFO column, specific keys -> values
    info: dict[str, float | int | str],
    additional: dict[str, dict[str, float | int | str]],  # on SAMPLEs columns
    # {sampleName: {filter dict}}
) -> list[VarianceData]:
    '''
    TODO
    filter: bool, # on FILTER column (PASS or not)
    info: dict[str, float|int|str], # on INFO column, specific keys -> values
    additional: dict[str, dict[str, float|int|str]], # on SAMPLEs columns
    # {sampleName: {filter dict}}

    Is this function worth making? bcftools actually have filter tool.
    '''
    return varianceDatas


def applyVarianceDataOnSeqRecord(
    varianceDatas: list[VarianceData], seqRecord: SeqRecord, searchFrom: int = 0
) -> SeqRecord:
    assert all(varianceData['CHROM'] == seqRecord.id
               for varianceData in varianceDatas)
    assert all(seqRecord.seq[varianceData['POS'] - 1].lower() ==
               varianceData['REF'][0].lower()
               for varianceData in varianceDatas), (
        'Reference nucleotide do not match variance data')

    newSeqRecord = seqRecord.copy()
    features = seqRecord.features
    # https://www.insdc.org/submitting-standards/feature-table/#7.3
    excludeQualifiers = [
        'codon_start', 'altitude', 'anticodon', 'artificial_location',
        'bio_material', 'bound_moiety', 'cell_line', 'cell_type', 'chromosome',
        'circular_RNA', 'citation', 'clone_lib', 'collected_by',
        'collection_date', 'compare', 'country', 'cultivar',
        'culture_collection', 'db_xref', 'dev_stage', 'direction', 'EC_number',
        'ecotype', 'environmental_sample', 'estimated_length', 'exception',
        'experiment', 'focus', 'frequency', 'gap_type', 'germline',
        'haplogroup', 'haplotype', 'host', 'identified_by', 'inference',
        'isolate', 'isolation_source', 'lab_host', 'lat_lon',
        'linkage_evidence', 'macronuclear', 'map', 'mating_type',
        'metagenome_source', 'mobile_element_type', 'mod_base', 'mol_type',
        'number', 'operon', 'organelle', 'organism', 'pop_variant',
        'ribosomal_slippage', 'rpt_family', 'rpt_type', 'rpt_unit_range',
        'satellite', 'serotype', 'serovar', 'sex', 'specimen_voucher',
        'strain', 'sub_clone', 'sub_species', 'sub_strain', 'tag_peptide',
        'tissue_lib', 'tissue_type', 'transgenic', 'translation',
        'transl_except', 'transl_table', 'trans_splicing' 'type_material',
        'variety'
    ]
    for varianceData in varianceDatas:
        start = varianceData['POS'] - 1
        length = len(varianceData['REF'])
        upStreamSeq = newSeqRecord[:start]
        downStreamSeq = newSeqRecord[start + length:]
        newSeqRecord = upStreamSeq + Seq(varianceData['ALT']) + downStreamSeq
        affectedFeatures = set(
            newSeqRecord.features
        ).symmetric_difference(set(features))
        mutFeatures = []
        diffLen = len(varianceData['ALT']) - len(varianceData['REF'])

        # add affected features back
        for feat in affectedFeatures:
            newLoc = FeatureLocation(
                feat.location.start,
                feat.location.end + diffLen,
                feat.location.strand
            )
            newFeat = feat
            newFeat.location = newLoc

            for qualifier, contents in newFeat.qualifiers.items():
                # white list
                if qualifier in excludeQualifiers:
                    continue
                newContents = []
                for content in contents:
                    newContents.append('MUTATED_' + content)
                newFeat.qualifiers[qualifier] = newContents
                if qualifier == 'translation' and newFeat.type == 'CDS':
                    newDna = newSeqRecord[
                        newFeat.location.start: newFeat.location.end
                    ]
                    if newFeat.location.strand == -1:
                        newDna = newDna.reverse_complement()
                    try:
                        translateTable: int | Literal['Standard'] = int(
                            newFeat.qualifiers['transl_table'][0]
                        )
                    except IndexError:
                        translateTable = 'Standard'
                    newFeat.qualifiers[qualifier] = [
                        newDna.translate(
                            translateTable, to_stop=True, cds=False
                        )
                    ]

            mutFeatures.append(newFeat)
        newSeqRecord.features.append(mutFeatures)

    # check if there is features in this region
    # for i, feat in enumerate(seqRecord.features[searchFrom:]):
    #     if len(feat) == len(seqRecord):
    #         continue
    #     if varianceData['POS'] < feat.location.start:
    #         continue
    #     elif feat.location.start <= varianceData['POS'] <= feat.location.end:
    #         newSearchFrom = searchFrom + i
    #     else:
    #         newSearchFrom = searchFrom + i


def applyVariancesOnSeqRecords(
    varianceDatas: list[VarianceData],
    seqRecordDict: dict[str, SeqRecord]
) -> dict[str, SeqRecord]:
    variantSeqRecordDict: dict[str, SeqRecord] = {}
    totalN = 0
    for seqid, seqRec in seqRecordDict.items():
        varianceDatas_specific = [
            varianceData for varianceData in varianceDatas
            if varianceData['CHROM'] == seqid
        ]
        for varianceData in varianceDatas_specific:

            pass
    for varianceData in varianceDatas:
        pass
    return variantSeqRecordDict
