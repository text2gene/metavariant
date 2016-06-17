
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

