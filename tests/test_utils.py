from __future__ import absolute_import, print_function, unicode_literals

import unittest

from metavariant.utils import strip_gene_name_from_hgvs_text


hgvs_text_with_gene_name = 'NM_006516.2(SLC2A1):c.1109T>C'
hgvs_text_no_gene_name = 'NM_006516.2:c.1109T>C'

class TestUtils(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_strip_gene_name_from_hgvs_text(self):
        result = strip_gene_name_from_hgvs_text(hgvs_text_no_gene_name)
        assert result == hgvs_text_no_gene_name

    def test_no_gene_name_in_hgvs_text_returns_same_string(self):
        result = strip_gene_name_from_hgvs_text(hgvs_text_no_gene_name)
        assert result == hgvs_text_no_gene_name

