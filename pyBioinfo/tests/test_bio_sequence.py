import unittest
from pathlib import Path

from pyBioinfo_modules.bio_sequences.vcf_parser import (
    vcfParser, VarianceData, essentialVcfColumns
)


class Test_vcfParser(unittest.TestCase):

    def setUp(self) -> None:
        self.standardMinColumns = [
            'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'
        ]

    def test_essentialColumns(self):
        self.assertListEqual(self.standardMinColumns, essentialVcfColumns)

    def test_minimal_vcf(self):
        minvcf_1 = Path('tests/test_data/vcf/minimal_1.vcf')
        minvcf_2 = Path('tests/test_data/vcf/minimal_2.vcf')
        minvcfData_1 = vcfParser(minvcf_1)
        minvcfData_2 = vcfParser(minvcf_2)
        oneVd = VarianceData(
            CHROM='CP009124',
            POS=2981667,
            ID='.',
            REF='T',
            ALT='C',
            QUAL=656.7,
            FILTER='PASS',
            INFO='AF=1.0000;AD=180;DP=180'
        )
        anotherVd = VarianceData(
            CHROM='CP009124',
            POS=4637796,
            ID='.',
            REF='CGGC',
            ALT='C',
            QUAL=0,
            FILTER='PASS',
            INFO='AF=1.0000;AD=170;DP=176'
        )
        self.assertEqual(len(minvcfData_1), 8)
        self.assertEqual(len(minvcfData_2), 17)
        self.assertIn(oneVd, minvcfData_1)
        self.assertIn(anotherVd, minvcfData_2)

    def test_complex_vcf(self):
        compVcf_1 = Path('tests/test_data/vcf/complex_1.vcf')
        compVcf_2 = Path('tests/test_data/vcf/complex_2.vcf')
        compvcfData_1 = vcfParser(compVcf_1)
        compvcfData_2 = vcfParser(compVcf_2)
        onecompVd = VarianceData(
            CHROM='CP009124.1',
            POS=2981667,
            ID='.',
            REF='T',
            ALT='C',
            QUAL=5858.77,
            FILTER='PASS',
            INFO=('AC=2;AF=1.00;AN=2;DP=143;ExcessHet=3.0103;FS=0.000;MLEAC=2;'
                  'MLEAF=1.00;MQ=60.00;QD=31.26;SOR=0.889'),
            FORMAT='GT:AD:DP:GQ:PL',
            SAMPLES={'mut': '1/1:0,143:143:99:5887,430,0'}
        )
        anotherVd = VarianceData(
            CHROM='CP009124.1',
            POS=3637029,
            ID='.',
            REF=('CCCCCGCATGCGGCCGCGCGGGGCGGCGGAACCTTTCCGCGGTCCGGGTGGGCCGGACCC'
                 'GCGGTCGGTTCCGCCG'),
            ALT='C',
            QUAL=1258.73,
            FILTER='PASS',
            INFO=('AC=1;AF=0.500;AN=2;BaseQRankSum=1.958;ClippingRankSum=0.000;'
                  'DP=117;ExcessHet=3.0103;FS=11.928;MLEAC=1;MLEAF=0.500;'
                  'MQ=60.00;MQRankSum=0.000;QD=11.34;ReadPosRankSum=-0.899;'
                  'SOR=0.490'),
            FORMAT='GT:AD:DP:GQ:PL',
            SAMPLES={'wt': '0/1:74,37:111:99:1296,0,2954'}
        )
        self.assertEqual(len(compvcfData_1), 18)
        self.assertEqual(len(compvcfData_2), 12)
        self.assertIn(onecompVd, compvcfData_1)
        self.assertIn(anotherVd, compvcfData_2)
        self.assertTrue(all('FORMAT' in vd for vd in compvcfData_1))
        self.assertTrue(all('SAMPLES' in vd for vd in compvcfData_1))
