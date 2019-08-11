
import unittest

from metavariant import VariantLVG
from metavariant.exceptions import CriticalHgvsError

SEQVAR_MAX_LEN_SAMPLE = 'NM_001204.6:c.51_814del'

class TestVariantLVG(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_del_numeric_or_chars(self):
        hgvs_text = 'NM_005228.3:c.2240_2257del18'
        lex = VariantLVG(hgvs_text)

        expected_c = [
                "NM_005228.3:c.2240_2257del",
                ]
        expected_g = 'NC_000007.14:g.55174777_55174794del'      # same: 'NC_000007.14:g.55174777_55174794delTAAGAGAAGCAACATCTC'
        expected_p = 'NP_005219.2:p.(Leu747_Pro753delinsSer)'
        expected_n = 'NM_005228.3:n.2486_2503del'               # same: 'NM_005228.3:n.2486_2503delTAAGAGAAGCAACATCTC'

        for c_hgvs_text in expected_c:
            assert c_hgvs_text in lex.hgvs_c

        assert expected_g in lex.hgvs_g
        assert expected_p in lex.hgvs_p
        assert expected_n in lex.hgvs_n
    
    def test_plain_lvgno_variant_mappings(self):
        #TODO
        hgvs_text = 'NM_194248.1:c.158C>T'

    def test_bad_hgvs_text_raises_CriticalHgvsError(self):
        bad_hgvs_text = 'NM_004628.4:c.621_622ins83'
    
        with self.assertRaises(CriticalHgvsError):
            VariantLVG(bad_hgvs_text)

    def test_gene_name(self):
        hgvs_text = 'NM_000548.3:c.826_827del'
        lex = VariantLVG(hgvs_text)
        assert lex.gene_name == 'TSC2'

    def test_seqvar_max_len_100(self):
        lex = VariantLVG(SEQVAR_MAX_LEN_SAMPLE, seqvar_max_len=100)
        for seqvar in lex.seqvars:
            assert len('%s' % seqvar) <= 100

    def test_seqvar_max_len_255(self):
        lex = VariantLVG(SEQVAR_MAX_LEN_SAMPLE, seqvar_max_len=255)
        for seqvar in lex.seqvars:
            assert len('%s' % seqvar) <= 255

    def test_to_json(self):
        lex = VariantLVG('NM_000020.2(ACVRL1):c.1031G>A')
        json_result = lex.to_json()
        assert 'NM_000020.2:c.1031G>A' in json_result
        assert 'NM_001077401.1' in json_result
        assert 'ACVRL1' in json_result

    def test_from_json(self):
        json_str = '{"hgvs_text": "NM_005228.3:c.2240_2257del18", "hgvs_n": ["NM_005228.3:n.2486_2503delTAAGAGAAGCAACATCTC"], "gene_name": "EGFR", "transcripts": ["NM_005228.3"], "hgvs_g": ["NC_000007.14:g.55174777_55174794delTAAGAGAAGCAACATCTC"], "hgvs_c": ["NM_005228.3:c.2240_2257del18"]}'
        lex = VariantLVG.from_json(json_str)
        assert lex.hgvs_text == 'NM_005228.3:c.2240_2257del18'
        assert lex.gene_name == 'EGFR'
        assert 'NC_000007.14:g.55174777_55174794delTAAGAGAAGCAACATCTC' in lex.hgvs_g
        assert 'NM_005228.3:n.2486_2503delTAAGAGAAGCAACATCTC' in lex.hgvs_n



lex = VariantLVG("NM_005228.3:c.2240_2257del18")

json_str = lex.to_json()
print(json_str)

lex = VariantLVG.from_json(json_str)


