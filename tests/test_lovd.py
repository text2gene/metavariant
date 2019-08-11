import unittest

from metavariant import LOVDVariantsForGene

class TestLOVD(unittest.TestCase):

    def test_workingness_of_LOVDVariantsForGene(self):
        # too experimental to trust results yet.
        # just a "test" to help remember to keep it maintained.
        variants = LOVDVariantsForGene('FOXP2')
        assert len(variants) > 1

