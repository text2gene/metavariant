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
LOVD_PAGE_GENE_VARIANTS_URL = 'http://{domain}/shared/view/{symbol}#page_size=1000&page=1'

re_pubmed_id = re.compile('pubmed\/(?P<pmid>\d+)')
re_doi = re.compile(r'(10[.][0-9]{2,}(?:[.][0-9]+)*/(?:(?!["&\'])\S)+)')


def construct_hgvs_name(transcript, coding_sequence):
    return '%s:%s' % (transcript, coding_sequence)


class LOVDVariant(object):

    def __init__(self, vdict=None, table_row=None, transcript=None, has_haplotype=False):

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
            self.__dict__.update(_get_vdict_from_tr(table_row, has_haplotype))
            if transcript:
                self.hgvs_text = construct_hgvs_name(transcript, self.cDNA)

    def load_from_api(self):
        """ Loads information from LOVD variant API, which does not return certain items like
        dbSNP number or references."""
        pass

    def load_from_page(self):
        """ Loads information from LOVD variant page, which is slower than an API call (and will
        probably break more readily over time), but includes references and other useful info."""
        pass

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


def _get_lovd_link_from_td(td_elem):
    return td_elem.find('a').get('href')


def _get_vdict_from_tr(tr_elem, has_haplotype=False):
    tds = tr_elem.getchildren()

    # there are 2 different styles of page. one has Haplotype at td-5, the other has no Haplotype column.
    # we have to parse these out differently.

    if has_haplotype:
        return { 'effect': tds[0].text,
                 'exon': tds[1].text,
                 'cDNA': _get_text_from_link_in_td(tds[2]), 
                 'lovd_link': _get_lovd_link_from_td(tds[2]),
                 'mRNA': tds[3].text,
                 'protein': tds[4].text,
                 'haplotype': tds[5].text,
                 'allele': tds[6].text,
                 'gDNA': _get_text_from_link_in_td(tds[7]),
                 'DBID': tds[10].text,
                 'references': _get_references_from_td(tds[12]),
                 'frequency': tds[16].text,
                 'remarks': tds[25].text,
                }
    else:
        return { 'effect': tds[0].text,
                 'exon': tds[1].text,
                 'cDNA': _get_text_from_link_in_td(tds[3]), 
                 'lovd_link': _get_lovd_link_from_td(tds[2]),
                 'mRNA': tds[4].text,
                 'protein': tds[5].text,
                 'gDNA': _get_text_from_link_in_td(tds[6]),
                 'DBID': tds[9].text,
                 'remarks': tds[10].text,
                 'references': _get_references_from_td(tds[11]),
                 'dbSNP': tds[11].text,
                 'genetic_origin': tds[12].text,
                 'frequency': tds[14].text,
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


def _query_lovd_api_for_variants_by_gene_name(symbol, domain):
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
        raise LOVDRemoteError('lovd.nl returned HTTP %r' % response.status_code)


def _load_lovd_variant_page_by_gene_name(symbol, domain):
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
        raise LOVDRemoteError('lovd.nl returned HTTP %r' % response.status_code)


def get_variants_for_gene_name(symbol, domain=LOVD_DEFAULT_DOMAIN):
    """ Given HUGO gene name symbol (e.g. "ACVRL1", "FANCA"), return the set of unique variants 
        found on LOVD.
    
    NOTE: Not all genes can be found on lovd.nl -- see http://databases.lovd.nl/shared/genes

    :param symbol: (str)
    :param domain: (str) [default: lovd.nl]
    :return: variants (list) 
    """
    api_response = _query_lovd_api_for_variants_by_gene_name(symbol, domain)
    vdicts = _parse_lovd_variants_by_gene_name_response(api_response)
    lovd_variants = set([vdict['hgvs_text'] for vdict in vdicts])
    return lovd_variants


def get_variants_with_annotations_for_gene_name(symbol, domain=LOVD_DEFAULT_DOMAIN):
    """ Given HUGO gene name symbol (e.g. "ACVRL1", "FANCA"), return the set of unique variants 
        found on LOVD paired with references (citations, annotations, etc).
    
    NOTE: Not all genes can be found on lovd.nl -- see http://databases.lovd.nl/shared/genes

    :param symbol: (str)
    :param domain: (str) [default: lovd.nl]
    :return: variants (list) 
    """
    content = _load_lovd_variant_page_by_gene_name(symbol, domain)
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

    # determine which type of table this is by what the headers look like
    has_haplotype = False
    ths = table.find('thead').getchildren()[0].findall('th')
    for th in ths:
        title = th.get('title')
        if title and 'haplotype' in th.get('title').lower():
            has_haplotype = True
    return [LOVDVariant(table_row=tr, transcript=tran, has_haplotype=has_haplotype) for tr in table.findall('tr')]


LOVDVariantsForGene = get_variants_for_gene_name

LOVDAnnotatedVariants = get_variants_with_annotations_for_gene_name

