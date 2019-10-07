import unittest

from metavariant.ncbi import NCBIVariationService, SPDI, HGVS 
from metavariant.hgvs_samples import hgvs_c, hgvs_g, hgvs_n, hgvs_p, intronic


# Tests for SPDI, HGVS, and any others we implement over time.

class TestNCBIVariationServices(unittest.TestCase):

    def test_NCBIVariationService_class(self):

        vs = NCBIVariationService()

    def test_SPDI_class(self):
        spdi = SPDI()
        #assert spdi.endpoint 

    def test_HGVS_class(self):
        hgvs = HGVS()
        #assert hgvs.endpoint

