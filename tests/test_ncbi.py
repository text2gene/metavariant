import unittest

from metavariant.ncbi import *

test_hgvs_c = 'NM_001232.3:c.919G>C'
test_hgvs_n = 'NM_194248.2:n.285C>T'


class TestNCBIEnrichedLVG(unittest.TestCase):
    # Stubs of tests to hold places for the various NCBI endpoints.

    def test_SPDI(self):
        pt = SPDI()


    def test_HGVS(self):
        pt = HGVS()

