from __future__ import absolute_import, print_function, unicode_literals

import sys, re

import requests
import xmltodict

from .exceptions import LOVDRemoteError

re_accession = re.compile('\/variants\/\w+\/(?P<acc_id>NM_\d+.\d+)\?')

LOVD_GENE_VARIANTS_URL = 'http://databases.lovd.nl/shared/api/rest.php/variants/%s'


def _parse_lovd_variants_by_gene_name_response(xml_blob, symbol):
    """ Parses LOVD response as provided by query_lovd_for_variants_by_gene_name

    :param xml_blob: (str) XML response from LOVD
    :param symbol: (str) gene name
    :return: dictionary of converted XML
    """
    variant_dicts = xmltodict.parse(xml_blob)['feed']['entry']
    out = []
    symbol = symbol.upper()

    for vari in variant_dicts:
        alt_link = (vari['link'][0]['@href'])
        accession_num = re_accession.findall(alt_link)[0]
        hgvs_text = vari['title'].replace(symbol, re_accession.findall(alt_link)[0])
        out.append(hgvs_text)

    return out


def _query_lovd_for_variants_by_gene_name(symbol):
    """ Submits query to lovd.nl for specified gene name (Symbol).

    :param symbol: (str)
    :return: unparsed response content (str)
    :raises: LOVDRemoteError if not response.ok (message contains response.status_code)
    """

    response = requests.get(LOVD_GENE_VARIANTS_URL % symbol)
    if response.ok:
        return response.content
    else:
        raise LOVDRemoteError('lovd.nl returned HTTP %r' % response.status)


def get_variants_for_gene_name(symbol):
    """ Given HUGO gene name symbol (e.g. "ACVRL1", "FANCA"), return list of variants found on LOVD.
    
    WARNING: Not all genes can be found on lovd.nl -- see http://databases.lovd.nl/shared/genes

    :param symbol: (str)
    :return: variants (list) 
    """
    api_response = _query_lovd_for_variants_by_gene_name(symbol)
    return _parse_lovd_variants_by_gene_name_response(api_response, symbol)


LOVDVariantsForGene = get_variants_for_gene_name


