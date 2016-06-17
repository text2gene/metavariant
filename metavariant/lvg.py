from __future__ import absolute_import, print_function, unicode_literals

import re
import logging

import hgvs.parser
import hgvs.variantmapper
from hgvs.exceptions import HGVSDataNotAvailableError, HGVSParseError

from .components import VariantComponents
from .config import get_uta_connection, PKGNAME
from .exceptions import CriticalHgvsError

log = logging.getLogger(PKGNAME)

# === UTA Connection setup. === #
uta = get_uta_connection()
mapper = hgvs.variantmapper.EasyVariantMapper(uta)
hgvs_parser = hgvs.parser.Parser()

# regular expressions
re_hgvs_with_gene_name = re.compile('^\w+\.\d+(?P<gene_name>\(\w+\)):')


def _seqvar_map_func(in_type, out_type):
    func_name = '%s_to_%s' % (in_type, out_type)
    return getattr(mapper, func_name)


def variant_to_gene_name(seqvar):
    """
    Get HUGO Gene Name (Symbol) for given sequence variant object.

    Input seqvar must be of type 'n', 'c', or 'p'.

    :param variant: hgvs.SequenceVariant
    :return: string gene name (or None if not available).
    """
    if seqvar.type in ['n', 'c', 'p']:
        try:
            tx_identity = uta.get_tx_identity_info(seqvar.ac)
        except HGVSDataNotAvailableError:
            return None

        if tx_identity is not None:
            return tx_identity[-1]
        else:
            return None
    else:
        return None


def _seqvar_to_seqvar(seqvar, base_type, new_type, transcript=None):
    if base_type == new_type:
        return None

    if base_type == 'p':
        return None

    if new_type == 'p' and base_type == 'g':
        return None

    map_seqvar = _seqvar_map_func(base_type, new_type)
    try:
        if base_type == 'g':
            if transcript:
                return map_seqvar(seqvar, transcript)
            else:
                return None
        return map_seqvar(seqvar)
    except NotImplementedError:
        log.debug('Cannot map %s to %s: hgvs raised NotImplementedError', seqvar, new_type)
        return None
    except HGVSDataNotAvailableError as error:
        log.debug('Cannot map %s to %s: hgvs raised HGVSDataNotAvailableError (%r)', seqvar, new_type, error)
        return None
    except Exception as error:
        # catch the general case to be robust around the hgvs library's occasional volatility.
        log.debug('Cannot map %s to %s: unexpected Exception (%r)', seqvar, new_type, error)
        return None


class VariantLVG(object):

    VERSION = '0.0.2'
    LVG_MODE = 'lvg'

    def __init__(self, hgvs_text_or_seqvar, **kwargs):
        self.hgvs_text = str(hgvs_text_or_seqvar)

        # use the hgvs library to get us some info about this HGVS string.
        self.seqvar = self.parse(hgvs_text_or_seqvar)

        if self.seqvar is None:
            raise CriticalHgvsError('Cannot create SequenceVariant from input %s (see hgvs_lexicon log)' % hgvs_text_or_seqvar)

        # initialize transcripts list
        self.transcripts = set(kwargs.get('transcripts', []))

        # fill in all the different ways to talk about this variant in each sequence type.
        self.variants = {'g': dict(), 'c': dict(), 'n': dict(), 'p': dict()}

        # collect any variants that were supplied at instantiation ("enrichment")
        for input_hgvs_c in kwargs.get('hgvs_c', []):
            self.variants['c'][str(input_hgvs_c)] = Variant(input_hgvs_c)

        for input_hgvs_g in kwargs.get('hgvs_g', []):
            self.variants['g'][str(input_hgvs_g)] = Variant(input_hgvs_g)

        for input_hgvs_n in kwargs.get('hgvs_n', []):
            self.variants['n'][str(input_hgvs_n)] = Variant(input_hgvs_n)

        for input_hgvs_p in kwargs.get('hgvs_p', []):
            self.variants['p'][str(input_hgvs_p)] = Variant(input_hgvs_p)

        try:
            self.variants[self.seqvar.type][str(self.seqvar)] = self.seqvar
        except KeyError:
            log.warn('Ignoring supplied SequenceVariant of type "%s" (not supported) -- (input was %s).' % (self.seqvar.type, self.seqvar))

        if self.variants['c']:
            # attempt to derive all 4 types of SequenceVariants from all available 'c'.
            for var_c in list(self.variants['c'].values()):
                for this_type, value in list(self.variants.items()):
                    new_seqvar = _seqvar_to_seqvar(var_c, 'c', this_type)
                    if new_seqvar:
                        self.variants[this_type][str(new_seqvar)] = new_seqvar

        # Now that we have a 'g', collect all available transcripts.
        if self.variants['g']:
            for var_g in list(self.variants['g'].values()):
                for trans in self.get_transcripts(var_g):
                    self.transcripts.add(trans)

        # In the case of starting with a 'g' type...
        if self.seqvar.type == 'g' and self.transcripts:
            # we still need to collect 'c' and 'n' variants
            for trans in self.transcripts:
                var_c = _seqvar_to_seqvar(self.seqvar, 'g', 'c', trans)
                if var_c:
                    self.variants['c'][str(var_c)] = var_c
                var_n = _seqvar_to_seqvar(self.seqvar, 'g', 'n', trans)
                if var_n:
                    self.variants['n'][str(var_n)] = var_n

        # With a list of transcripts, we can do g_to_c and g_to_n
        if self.transcripts:
            for trans in self.transcripts:
                for var_g in list(self.variants['g'].values()):
                    # Find all available 'c'
                    if not trans.startswith('NR'):
                        var_c = _seqvar_to_seqvar(var_g, 'g', 'c', trans)
                        if var_c:
                            self.variants['c'][str(var_c)] = var_c

                    # Find all available 'n'
                    var_n = _seqvar_to_seqvar(var_g, 'g', 'n', trans)
                    if var_n:
                        self.variants['n'][str(var_n)] = var_n

        # map all newly found 'c' to 'p'
        for var_c in list(self.variants['c'].values()):
            new_seqvar = _seqvar_to_seqvar(var_c, 'c', 'p')
            if new_seqvar:
                self.variants['p'][str(new_seqvar)] = new_seqvar

        # find out the gene name of this variant.
        self._gene_name = None

    @property
    def gene_name(self):
        """ Lazy-loaded gene name based on hgvs lookup. """
        if self._gene_name is None:
            for seqvar in self.variants['c'].values() + self.variants['n'].values() + self.variants['p'].values():
                name = variant_to_gene_name(seqvar)
                if name:
                    self._gene_name = name
                    break
        return self._gene_name

    @staticmethod
    def get_transcripts(var_g):
        return mapper.relevant_transcripts(var_g)

    @staticmethod
    def parse(hgvs_text_or_seqvar):
        """ Parse input through hgvs_parser if text, do nothing if SequenceVariant.
        Return SequenceVariant object.

        Allow all potential hgvs parsing errors and type errors to flow upwards.

        :param hgvs_text_or_seqvar: string or SequenceVariant object
        :return: SequenceVariant object
        """
        if type(hgvs_text_or_seqvar) == hgvs.variant.SequenceVariant:
            return hgvs_text_or_seqvar

        # see if the gene name is embedded in the string, e.g. NM_003331.4(TYK2):c.3318_3319insC
        match = re_hgvs_with_gene_name.match(hgvs_text_or_seqvar) 
        if match:
            gene_str = match.groupdict()['gene_name']
            hgvs_text = hgvs_text_or_seqvar.replace(gene_str, '')
        else:
            hgvs_text = hgvs_text_or_seqvar

        try:
            return hgvs_parser.parse_hgvs_variant(str(hgvs_text))
        except HGVSParseError as error:
            log.info('Cannot create SequenceVariant from hgvs_text "%s": %r', hgvs_text, error)
            # Examples:
            #  HGVSParseError(u'NP_068780.2:p.Tyr?His: char 17: expected a digit',)
            #  HGVSParseError(u'NM_004628.4:c.621_622ins83: char 24: Syntax error',)
            return None

    @property
    def hgvs_c(self):
        return self.variants['c'].keys()

    @property
    def hgvs_g(self):
        return self.variants['g'].keys()

    @property
    def hgvs_p(self):
        return self.variants['p'].keys()

    @property
    def hgvs_n(self):
        return self.variants['n'].keys()

    @property
    def seqvars(self):
        """ returns a flat list of all SequenceVariant objects contained in
        the self.variants dictionary """
        out = []
        for seqvar_dict in (self.variants.values()):
            out = out + list(seqvar_dict.values())
        return out

    def to_dict(self, with_gene_name=True):
        """Returns contents of object as a 2-level dictionary.

        Supply with_gene_name = True [default: True] to return gene_name as well.

        (gene_name is a lazy-loaded magic attribute, and may take a second or two).
        """
        outd = {'variants': {},
                'transcripts': list(self.transcripts),
                'seqvar': self.seqvar,
                'hgvs_text': self.hgvs_text,
               }

        if with_gene_name:
            outd['gene_name'] = self.gene_name

        # turn each seqvar dictionary into just a list of its values (the seqvar objects).
        for seqtype, seqvar_dict in (self.variants.items()):
            outd['variants'][seqtype] = list(seqvar_dict.values())

        return outd

    def __str__(self):
        out = 'HGVS input: %s\n' % self.hgvs_text
        out += '%r' % self.seqvar
        return out


### API Convenience Functions

Variant = VariantLVG.parse

