""" Provides NCBIEnrichedLVG object NCBI Variant Services client. """

import requests
import urllib
import time

from .lvg import VariantLVG, Variant
from .exceptions import CriticalHgvsError, NCBIRemoteError  #, MetaVariantException
from .utils import strip_gene_name_from_hgvs_text

# ACKNOWLEDGEMENTS:
# 
# This module was built from example code provided by the NCBI, written by SPDI format paper 
# co-author Lon Phan (lonphan@ncbi.nlm.nih.gov).
# 
#   https://github.com/ncbi/dbsnp/blob/master/tutorials/Variation%20Services/spdi_batch.py
#
# spdi_batch Jupyter notebook example:
#   https://hub.gke.mybinder.org/user/ncbi-dbsnp-wogm4vv0/notebooks/tutorials/Variation%20Services/Jupyter_Notebook/spdi_batch.ipynb
# 
# SPDI format specification (paper alluded to above):
# https://www.biorxiv.org/content/10.1101/537449v3.full


def api_request(url, retry=0, wait=.1):
    """ Generalized function for request via API expecting a JSON return.

    If retry is set to a positive integer, this function recurses for each retry
    until retry decrements to 0.

    :param url: preconstructed URL to request via HTTP
    :param retry: (int) how many times to retry in case of a timeout.
    :param wait: (float) how long in seconds to wait between retries
    :return: response object (from requests lib) if status_code was 200.
    :raises: NCBIRemoteError
    """
    try:
        res = requests.get(url)
    except requests.exceptions.Timeout:
        # Maybe set up for a retry, or continue in a retry loop
        if retry > 0:
            return api_request(url, retry-1)
        else:
            raise NCBIRemoteError('Timeout during api_request to %s' % url)
    except requests.exceptions.TooManyRedirects:
        # Tell the user their URL was bad and try a different one
        raise NCBIRemoteError('Bad URL / too many redirects: %s' % url)
    except requests.exceptions.RequestException as err:
        # catastrophic error. bail.
        raise NCBIRemoteError('Oh noes! Error while loading URL {}: {}'.format(url, error))
    if (res.status_code == 200):
        return res
    else:
        raise NCBIRemoteError('ERROR {} while loading URL {}'.format(res.status_code, url))


class NCBIVariantService(object):
    """ Base class for NCBI Variant Services API endpoints. """

    api_rootURL = 'https://api.ncbi.nlm.nih.gov/variation/v0/'

    def compose_url(self, suffix):
        return self.api_rootURL + suffix

    def fetch(self, url, retry=0):
        return api_request(url, retry)

    def parse(self, data):
        print(data)

    def batch(self, long_list):
        print(long_list)


class SPDIgonzales(NCBIVariantService):
    """ Handles the SPDI endpoint of the NCBI Variant Services API. """
    # /spdi/{spdi}/*  [GET]
    #           contextual, vcf_fields, canonical_representative, all_equivalent_contextual, rsids

    endpoint = 'spdi'

    def contextual(self, spdi, retry=0):
        suffix = 'contextual'
        url = self.compose_url(spdi, suffix)
        return self.fetch(url, retry)

    def vcf_fields(self, spdi, retry=0):
        suffix = 'vcf_fields'
        url = self.compose_url(spdi, suffix)
        return self.fetch(url, retry)

    def 



# endpoints:
#
# /spdi/{spdi}/*  [GET]
#           contextual, vcf_fields, canonical_representative, all_equivalent_contextual, rsids
#
#
# /hgvs/{hgvs}/contextuals  [GET]
# /hgvs/batch/contextuals   [POST]
#
# /vcf/{chrom}/{pos}/{ref}/{alts}/contextuals   [GET]]
# /vcf/file/set_rsids       [POST]
#   q                                                                                                                                                       
# /refsnp/{rsid}            [GET]
# /



class NCBIEnrichedLVG(VariantLVG):

    """ Creates a true LVG object by subclassing from VariantLVG and using data drawn from an NCBIReport
    lookup.  See VariantLVG (from metavariant) for additional documentation.

    Examples:
        lex = NCBIEnrichedLVG('NM_000249.3:c.1958T>G')

        seqvar = Variant('NM_000249.3:c.1958T>G')
        lex = NCBIEnrichedLVG(seqvar)
    """

    VERSION = 2
    LVG_MODE = 'ncbi_enriched'

    def __init__(self, hgvs_text_or_seqvar, **kwargs):
        self.hgvs_text = strip_gene_name_from_hgvs_text('%s' % hgvs_text_or_seqvar)
        self.seqvar = Variant(hgvs_text_or_seqvar)
        self.ncbierror = None
        if self.seqvar is None:
            raise CriticalHgvsError('Cannot create SequenceVariant from input %s' % self.hgvs_text)
        try:
            report = get_ncbi_variant_report(self.hgvs_text)
            self.variants = ncbi_report_to_variants(report)
        except NCBIRemoteError as error:
            log.debug('Skipping NCBI enrichment; %r' % error)
            self.ncbierror = '%r' % error
            self.error = ''
            self.variants = {'c': {}, 'g': {}, 'p': {}, 'n': {}}
            self.variants[self.seqvar.type][self.hgvs_text] = self.seqvar

        super(NCBIEnrichedLVG, self).__init__(self.hgvs_text,
                                              hgvs_c=self.hgvs_c,
                                              hgvs_g=self.hgvs_g,
                                              hgvs_p=self.hgvs_p,
                                              hgvs_n=self.hgvs_n,
                                              **kwargs)

# ---- RETURN FORMATS ----
#
# spdi{
# description:	
# A single contextual allele in SPDI notation. Contextual allele means that applying the Blossom Precision Correction Algorithm would leave the fields unchanged.
#
# seq_id*	string($ascii)
# example: NC_000001.23
# The RefSeq/Genbank Accession.Version for the reference sequence
#
# position*	integer
# The 0-based boundary position where the deletion starts. That is, position 0 starts the deletion immediately before the first nucleotide and position 1 starts the deletion between the first and second nucleotides.
#
# deleted_sequence*	string($ascii)
# The IUPAC sequence of nucleotides/amino-acids to delete from the reference. This can be empty, which is how a pure insertion is represented.
#
# inserted_sequence*	string($ascii)
# The IUPAC sequence of nucleotides/amino-acids to insert after perforing the deletion. Amino-acids use the single character encoding. All nucleotide codes including the ones for ambiguity are allowed.
#
# warnings	warnings[...]
# }

# 


