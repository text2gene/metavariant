import unittest

from metavariant import NCBIEnrichedLVG

test_hgvs_c = 'NM_001232.3:c.919G>C'
test_hgvs_n = 'NM_194248.2:n.285C>T'


class TestNCBIEnrichedLVG(unittest.TestCase):

    def test_NCBIEnrichedLVG(self):
        # Stub of a test to flesh out over time as we find edge cases and whatnot.
        lvg = NCBIEnrichedLVG(test_hgvs_c)
        assert lvg.hgvs_c[0] == test_hgvs_c
        assert len(lvg.hgvs_g) > 0

        lvg = NCBIEnrichedLVG(test_hgvs_n)
        assert len(lvg.seqvars) == 1

