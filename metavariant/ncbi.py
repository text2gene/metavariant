""" Provides NCBIEnrichedLVG object and NCBI Variant report functions. """

from .lvg import VariantLVG, Variant
from .exceptions import CriticalHgvsError, NCBIRemoteError  #, MetaVariantException
from .utils import strip_gene_name_from_hgvs_text




def ncbi_report_to_variants(report):
    """ Parses Hgvs_* strings from NCBI report and creates a "variants" dictionary
    like the following (mimicking the VariantLVG.variants attribute):

    {seqtype: { 'hgvs_string': SequenceVariant object }

    :param report: list of strings representing NCBI Variation Reporter output
    :return: dict as per structure above
    """
    variants = {'p': {}, 'c': {}, 'g': {}, 'n': {}, 'm': {}, 'r': {}}
    for rep_part in report:
        for seqtype in variants.keys():
            hgvs_text = rep_part.get('Hgvs_%s' % seqtype, '').strip()
            if hgvs_text:
                # set up data structure just like VariantLVG object, i.e.:
                # {seqtype: { 'hgvs_string': SequenceVariant object }
                seqvar = Variant(hgvs_text)
                if seqvar:
                    # Sometimes NCBI has variant strings that do not parse. Variant() function returns None in these cases.
                    variants[seqvar.type][str(seqvar)] = seqvar
    return variants


def ncbi_report_to_pubmeds(report):
    """ Parses PMIDs from NCBI report and returns as list of strings.

    :param report: list of strings representing NCBI Variation Reporter output
    :return: list of pubmeds found in report
    """
    return [int(item) for item in report[0]['PMIDs']]



def get_ncbi_variant_report(hgvs_text):
    """
    Return results from API query to the NCBI Variant Reporter Service
    See documentation at:
    https://www.ncbi.nlm.nih.gov/variation/tools/reporter

    :param hgvs_text: ( c.DNA | r.RNA | p.Protein | g.Genomic )
    :return: list containing each dict of parsed results
    """
    response = requests.get("https://www.ncbi.nlm.nih.gov/projects/SNP/VariantAnalyzer/var_rep.cgi?annot1={}".format(urllib.quote(hgvs_text)))

    if 'Error' in response.text:
        error_str = 'The NCBI Variant Report Service returned an error: "{}"\n'.format(response.text)
        error_str += 'To reproduce, visit: https://www.ncbi.nlm.nih.gov/projects/SNP/VariantAnalyzer/var_rep.cgi?annot1={}'.format(hgvs_text)
        raise NCBIRemoteError(error_str)

    keys = []
    report = []

    for line in response.text.split('\n'):
        if not line.strip() or line.startswith('.') or line.startswith('##') or line.startswith('Submitted'):
            continue

        if line.startswith('# '):
            keys = line.strip('# ').split('\t')
        else:
            values = line.split('\t')
            outd = dict(zip(keys, values))

            # convert PMIDs from semicolon- or comma-separated string into python list
            if len(outd.get('PMIDs', '')) > 0:
                outd['PMIDs'] = outd['PMIDs'].replace(', ', ';').split(';')
            report.append(outd)

    if report == []:
        raise NCBIRemoteError('The NCBI Variant Report Service returned an empty report.\nTo reproduce, visit: https://www.ncbi.nlm.nih.gov/projects/SNP/VariantAnalyzer/var_rep.cgi?annot1={}'.format(hgvs_text))

    return report



class NCBIEnrichedLVG(VariantLVG):

    """ Creates a true LVG object by subclassing from VariantLVG and using data drawn from an NCBIReport
    lookup.  See VariantLVG (from metavariant) for additional documentation.

    *** ! To use the cached version of this object, use LVGEnriched ! ***

    Examples:
        lex = NCBIEnrichedLVG('NM_000249.3:c.1958T>G')

        seqvar = Variant('NM_000249.3:c.1958T>G')
        lex = NCBIEnrichedLVG(seqvar)
    """

    VERSION = 1
    LVG_MODE = 'ncbi_enriched'

    def __init__(self, hgvs_text_or_seqvar, **kwargs):
        self.hgvs_text = strip_gene_name_from_hgvs_text('%s' % hgvs_text_or_seqvar)
        self.seqvar = Variant(hgvs_text_or_seqvar)
        self.ncbierror = None
        if self.seqvar is None:
            raise CriticalHgvsError('Cannot create SequenceVariant from input %s' % self.hgvs_text)
        try:
            report = NCBIReport(self.hgvs_text)
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


