""" Provides NCBIEnrichedLVG object NCBI Variant Services client. """

import requests
#import urllib
import time
import json

from .config import log
from .lvg import VariantLVG, Variant
from .exceptions import CriticalHgvsError, NCBIRemoteError  #, MetaVariantException
from .utils import strip_gene_name_from_hgvs_text

# TODO: allow configuration through .config via env variables.
NCBI_RETRIES = 1
NCBI_WAIT_TIME = 1.0

# this URL will probably stay hard-coded here.
NCBI_VS_API_ROOT = 'https://api.ncbi.nlm.nih.gov/variation/v0/'

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


def api_request(url, retry=0, wait=1.0, data=None):
    """ Generalized function for shipping a GET or POST request via HTTP.

    If retry is set to a positive integer, this function recurses for each retry
    until retry decrements to 0.

    POST requests are automatically set when data is set.  You can set data to
    an empty dictionary if you have no data to send.

    :param url: preconstructed URL to request via HTTP
    :param retry: (int) how many times to retry in case of a timeout.
    :param wait: (float) how long in seconds to wait between retries
    :param data: [default: None] if set, switches to POST method.
    :return: response object (from requests lib) if status_code was 200.
    :raises: NCBIRemoteError
    """
    try:
        if data is not None:
            res = requests.post(url, data=data)
        else:
            res = requests.get(url)
    except requests.exceptions.Timeout:
        # Maybe set up for a retry, or continue in a retry loop
        if retry > 0:
            time.sleep(wait)
            return api_request(url, retry-1)
        else:
            raise NCBIRemoteError('Timeout during api_request to %s' % url)
    except requests.exceptions.TooManyRedirects:
        # Tell the user their URL was bad and try a different one
        raise NCBIRemoteError('Bad URL / too many redirects: %s' % url)
    except requests.exceptions.RequestException as err:
        # catastrophic error. bail.
        raise NCBIRemoteError('Oh noes! Error while loading URL {}: {}'.format(url, err))
    if (res.status_code == 200):
        return res
    else:
        raise NCBIRemoteError('ERROR {} while loading URL {}'.format(res.status_code, url))


class NCBIVariationService(object):
    """ Base class for NCBI Variant Services API endpoints. """

    _endpoint = None

    def _compose_url(self, text, suffix=''):
        """ Composes an API url based on the global NCBI_VS_API_ROOT
            plus the subclass's class variable _endpoint
            plus the supplied data query (text)
            plus any relevant suffix (e.g. 'vcf_fields').

        Raises TypeError if the subclass's _endpoint class variable has not been set.

        :param text: typically the data point being submitted to the _endpoint for query.
        :param suffix: a qualifier for the _endpoint  [default: '']
        :return url: composed url suitable for submission via GET or POST
        :rtype: str
        :raises: TypeError
        """
        #TODO: HTTP encode text.  Use urljoin instead of concats.
        return NCBI_VS_API_ROOT + self._endpoint + '/' + text + '/' + suffix

    def _fetch(self, url, retry=NCBI_RETRIES, wait=NCBI_WAIT_TIME):
        return api_request(url, retry, wait)

    def _fetch_POST(self, url, data, retry=NCBI_RETRIES, wait=NCBI_WAIT_TIME):
        return api_request(url, retry, wait, data)


# aka SPDI_Gonzales :-D 
class SPDI(NCBIVariationService):
    """ Handles the SPDI endpoint of the NCBI Variant Services API. 

    /spdi/{spdi}/*  [GET]
    * contextual
    * vcf_fields
    * canonical_representative
    * all_equivalent_contextual
    * rsids
    """
    _endpoint = 'spdi'

    def contextual(self, spdi, retry=0):
        suffix = 'contextual'
        url = self._compose_url(spdi, suffix)
        return self._fetch(url, retry)

    def vcf_fields(self, spdi, retry=0):
        suffix = 'vcf_fields'
        url = self._compose_url(spdi, suffix)
        return self._fetch(url, retry)

    def canonical_representative(self, spdi, retry=0):
        suffix = 'canonical_representative'
        url = self._compose_url(spdi, suffix)
        return self._fetch(url, retry)

    def all_equivalent_contextual(self, spdi, retry=0):
        suffix = 'all_equivalent_contextual'
        url = self._compose_url(spdi, suffix)
        return self._fetch(url, retry) 

    def rsids(self, spdi, retry=0):
        suffix = 'rsids'
        url = self._compose_url(spdi, suffix)
        try:
            res = self._fetch(url, retry)
            return json.loads(req.text)['data']['rsids']
        except NCBIRemoteError as err:
            log.debug(err)
            return None


class HGVS(NCBIVariationService):
    # /hgvs/{hgvs}/contextuals  [GET]
    # /hgvs/batch/contextuals   [POST]

    _endpoint = 'hgvs'

    def contextuals(self, hgvs_text, retry=0):
        suffix = 'contextuals'
        url = self._compose_url(hgvs_text, suffix)
        return self._fetch(url, retry)

    def batch_contextuals(self, data, retry=0):
        url = self._compose_url('batch', 'contextuals')
        return self._fetch_POST(url, data, retry)

# endpoints:
#
# /vcf/{chrom}/{pos}/{ref}/{alts}/contextuals   [GET]]
# /vcf/file/set_rsids       [POST]
#   q                                                                                                                                                       
# /refsnp/{rsid}            [GET]
# /


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


