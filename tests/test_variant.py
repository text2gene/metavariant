
import unittest

from metavariant import Variant
from metavariant.exceptions import CriticalHgvsError

class TestVariant(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_allow_gene_name_in_hgvs_string(self):
        seqvar = Variant('NM_001165963.1(SCN1A):c.4338_6030del')
        assert 'NM_001165963.1:c.4338_6030del' == '%s' % seqvar    

    def test_bad_hgvs_string_returns_None(self):
        assert Variant('boogers') is None

    def test_crazy_long_hgvs_strings_just_to_see_what_happens(self):

        # https://www.ncbi.nlm.nih.gov/clinvar/38771295/    RCV000008537
        seqvar = Variant('NP_001121636.1:p.Gln197_Gln208delinsGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGlnGln')        # https://www.ncbi.nlm.nih.gov/clinvar/36484151/

        # https://www.ncbi.nlm.nih.gov/clinvar/variation/68075/
        seqvar = Variant('LRG_702p1:p.Glu129_Asp145delinsGluIleLysValIleSerGlyIleLeuThrGlnGlyArgCysAspIleGluIleLysValIleSerGlyIleLeuThrGlnGlyArgCysAspIleAsp')

        # RCV000087437
        seqvar = Variant('LRG_3p1:p.Arg1139_Gly1140insValSerSerThrGluArgTyrTyrArgSerThrCysPheArgCysLeuHisPheArgLysIlePheTrpHisCysAspValMetIleLeuSerLeu')

        # also RCV000087437
        seqvar = Variant('NP_000081.1:p.Arg1139_Gly1140insValSerSerThrGluArgTyrTyrArgSerThrCysPheArgCysLeuHisPheArgLysIlePheTrpHisCysAspValMetIleLeuSerLeu')

