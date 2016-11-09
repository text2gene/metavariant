from __future__ import absolute_import, print_function

import unittest

from metavariant.components import VariantComponents, findall_aminochanges_in_text, parse_components_from_aminochange


AA_synonyms = {'Leu653Arg': ['L653R', 'Leu653Arg', '(Leu653Arg)'],
               'Cys344Tyr': ['(Cys344Tyr)', 'C344Y', 'Cys344Tyr'],
              }

good = 0
bad = 0
total = 6


class TestAminoChangeSearch(unittest.TestCase):

    def test_findall_aminochanges_in_text(self):
        for a_chg in AA_synonyms.keys():
            found = findall_aminochanges_in_text('%r' % AA_synonyms[a_chg])
            for synonym in AA_synonyms[a_chg]:
                assert synonym in found

    def test_parse_components_from_aminochange(self):
        for a_chg_list in AA_synonyms.values():
            for a_chg in a_chg_list:
                assert parse_components_from_aminochange(a_chg) is not None

    def test_posedits_from_aminochange(self):
        achg = 'Leu653Arg'
        comp = VariantComponents(aminochange=achg)
        assert comp.posedit == achg
        assert 'L653R' in comp.posedit_slang


