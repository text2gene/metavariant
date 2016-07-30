from __future__ import absolute_import, print_function, unicode_literals

import sys, re

import requests
import xmltodict

from lxml import etree
from lxml.html import HTMLParser

from .exceptions import LOVDRemoteError

re_transcript = re.compile('\/variants\/\w+\/(?P<transcript>NM_\d+.\d+)\?')

LOVD_DEFAULT_DOMAIN = 'databases.lovd.nl'

LOVD_API_GENE_VARIANTS_URL = 'http://{domain}/shared/api/rest.php/variants/{symbol}'
LOVD_PAGE_GENE_VARIANTS_URL = 'http://databases.lovd.nl/shared/view/{symbol}#page_size=1000&page=1'

re_pubmed_id = re.compile('pubmed\/(?P<pmid>\d+)')
re_doi = re.compile(r'(10[.][0-9]{2,}(?:[.][0-9]+)*/(?:(?!["&\'])\S)+)')


def construct_hgvs_name(transcript, coding_sequence):
    return '%s:%s' % (transcript, coding_sequence)


class LOVDVariant(object):

    def __init__(self, vdict=None, table_row=None, transcript=None):

        self.gene_name = None
        self.title = None
        self.transcript = transcript
        self.hgvs_text = None
        self.id = None
        self.DBID = None
        self.cDNA = None
        self.gDNA = None
        self.mRNA = None
        self.times_reported = None
        self.link = None
        self.alt_link = None
        self.references = {'pmid': [], 'doi': []}

        if vdict is not None:
            self.__dict__.update(vdict)
        elif table_row is not None:
            self.__dict__.update(_get_vdict_from_tr(table_row))
            if transcript:
                self.hgvs_text = construct_hgvs_name(transcript, self.cDNA)

    def to_dict(self):
        return self.__dict__    


def _get_references_from_td(td_elem):
    refs = {'doi': [], 'pmid': []}
    for span in td_elem.findall('span'):
        if not span.text:
            continue

        if span.text.startswith('Journal'):
            stuff = span.get('onmouseover')
            doi = re_doi.findall(stuff)[0].strip('\\')
            if doi:
                refs['doi'].append(doi)
            
        if span.text.startswith('PubMed'):
            if span.text.startswith('PubMed'):
                stuff = span.get('onmouseover')
                pmid = re_pubmed_id.findall(stuff)[0]
                if pmid:
                    refs['pmid'].append(pmid)
    return refs


def _get_text_from_link_in_td(td_elem):
    try:
        return td_elem.find('a').text
    except AttributeError:
        return ''


def _get_vdict_from_tr(tr_elem):
    tds = tr_elem.getchildren()
   
    return { 'allele': tds[6].text,
             'cDNA': _get_text_from_link_in_td(tds[2]), 
             'exon': tds[1].text,
             'effect': tds[0].text,
             'frequency': tds[16].text,
             'gDNA': _get_text_from_link_in_td(tds[7]),
             'haplotype': tds[5].text,
             'DBID': tds[10].text,
             'protein': tds[4].text,
             'references': _get_references_from_td(tds[12]),
             'remarks': tds[25].text,
             'mRNA': tds[3].text
            }


def _parse_content_text(text):
    """ <content type="text">
          symbol:FANCA
          id:0000116087
          position_mRNA:c.1626+1_1627-1
          position_genomic:chr16:?
          Variant/DNA:c.(1626+1_1627-1)_(2151+1_2152-1)del
          Variant/DBID:FANCA_000723
          Times_reported:1
        </content>
    """
    out = {}
    for part in text.split('\n'):
        part = part.strip()
        try:
            key, val = part.split(':', 1)
        except IndexError:
            continue

        if part.startswith('symbol'):
            out['gene_name'] = val
        if part.startswith('id'):
            out['id'] = val
        if part.startswith('Variant/DBID'):
            out['DBID'] = val
        if part.startswith('Variant/DNA'):
            out['cDNA'] = val
        if part.startswith('position_mRNA'):
            out['mRNA'] = val
        if part.startswith('position_genomic'):
            out['gDNA'] = val
        if part.startswith('Times_reported'):
            out['times_reported'] = val
    return out


def _parse_entry(entry):
    vdict = _parse_content_text(entry['content']['#text'])
    vdict['link'] = entry['link'][1]['@href']
    vdict['alt_link'] = (entry['link'][0]['@href'])
    vdict['title'] = entry['title']
    vdict['transcript'] = re_transcript.findall(vdict['alt_link'])[0]
    vdict['hgvs_text'] = construct_hgvs_name(vdict['transcript'], entry['title'].replace(vdict['gene_name']+ ':', ''))
    return vdict


def _parse_lovd_variants_by_gene_name_response(xml_blob):
    """ Parses LOVD response as provided by query_lovd_for_variants_by_gene_name

    Returns list of dictionaries containing pared-down and interpreted variant info.

    :param xml_blob: (str) XML response from LOVD
    :return: list of dictionaries (one per variant entry)
    """
    entries = xmltodict.parse(xml_blob)['feed']['entry']
    out = []
    for entry in entries:
        out.append(_parse_entry(entry))
    return out


def _query_lovd_for_variants_by_gene_name(symbol, domain=LOVD_DEFAULT_DOMAIN):
    """ Submits query to lovd.nl for specified gene name (Symbol), retrieving the API
    content for the given gene.

    :param symbol: (str)
    :param domain: (str) domain name of target LOVD database (default: databases.lovd.nl)
    :return: unparsed response content (str)
    :raises: LOVDRemoteError if not response.ok (message contains response.status_code)
    """

    response = requests.get(LOVD_API_GENE_VARIANTS_URL.format(symbol=symbol, domain=domain))
    if response.ok:
        return response.content
    else:
        raise LOVDRemoteError('lovd.nl returned HTTP %r' % response.status)


def _load_lovd_variant_page_by_gene_name(symbol, domain=LOVD_DEFAULT_DOMAIN):
    """ Submits query to lovd.nl for specified gene name (Symbol), retrieving the WEBSITE
    content for the given gene.

    :param symbol: (str)
    :param domain: (str) domain name of target LOVD database (default: databases.lovd.nl)
    :return: unparsed response content (str)
    :raises: LOVDRemoteError if not response.ok (message contains response.status_code)
    """
    response = requests.get(LOVD_PAGE_GENE_VARIANTS_URL.format(symbol=symbol, domain=domain))
    if response.ok:
        return response.content
    else:
        raise LOVDRemoteError('lovd.nl returned HTTP %r' % response.status)


def get_variants_for_gene_name(symbol):
    """ Given HUGO gene name symbol (e.g. "ACVRL1", "FANCA"), return the set of unique variants 
        found on LOVD.
    
    NOTE: Not all genes can be found on lovd.nl -- see http://databases.lovd.nl/shared/genes

    :param symbol: (str)
    :return: variants (list) 
    """
    api_response = _query_lovd_for_variants_by_gene_name(symbol)
    vdicts = _parse_lovd_variants_by_gene_name_response(api_response)
    lovd_variants = set([vdict['hgvs_text'] for vdict in vdicts])
    return lovd_variants


def get_variants_with_annotations_for_gene_name(symbol):
    """ Given HUGO gene name symbol (e.g. "ACVRL1", "FANCA"), return the set of unique variants 
        found on LOVD paired with references (citations, annotations, etc).
    
    NOTE: Not all genes can be found on lovd.nl -- see http://databases.lovd.nl/shared/genes

    :param symbol: (str)
    :return: variants (list) 
    """
    content = _load_lovd_variant_page_by_gene_name(symbol)
    htm = etree.fromstring(content, parser=HTMLParser()).find('body')
    try:
        table = htm.cssselect('#viewlistTable_CustomVL_VIEW_%s' % symbol)[0]
    except:
        table = htm.cssselect('#viewlistTable_CustomVL_VOTunique_VOG_%s' % symbol)[0]
    transcripts = re.findall('using the (?P<acc_id>NM_\d+\.\d+) transcript reference sequence.', content)
    if transcripts:
        tran = transcripts[0]
    else:
        tran = None 

    return [LOVDVariant(table_row=tr, transcript=tran) for tr in table.findall('tr')]


LOVDVariantsForGene = get_variants_for_gene_name

LOVDAnnotatedVariants = get_variants_with_annotations_for_gene_name

