from __future__ import absolute_import, print_function, unicode_literals

# utils.py contains utility functions that should require few, if any, 3rd party python libraries.

import re

# regular expressions
re_hgvs_with_gene_name = re.compile('^\w+\.\d+(?P<gene_name>\(\w+\)):')


def strip_gene_name_from_hgvs_text(hgvs_text):
    """ If gene name is embedded in the hgvs_text string, strip it out and
    return the HGVS string without the gene name part.

    :param: hgvs_text (str)
    :return: hgvs_text with gene_name removed (if needed)
    :rtype: str
    """
    # see if the gene name is embedded in the string, e.g. NM_003331.4(TYK2):c.3318_3319insC
    match = re_hgvs_with_gene_name.match(hgvs_text)
    if match:
        gene_str = match.groupdict()['gene_name']
        return hgvs_text.replace(gene_str, '')
    else:
        return hgvs_text
