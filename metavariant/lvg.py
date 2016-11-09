from __future__ import absolute_import, print_function, unicode_literals

import re
import json
import logging

import hgvs.parser
import hgvs.variantmapper
from hgvs.exceptions import HGVSDataNotAvailableError, HGVSParseError

from .components import VariantComponents
from .config import get_uta_connection, PKGNAME
from .exceptions import CriticalHgvsError, RejectedSeqVar
from .utils import strip_gene_name_from_hgvs_text

log = logging.getLogger(PKGNAME)

# === UTA Connection setup. === #
uta = get_uta_connection()
mapper = hgvs.variantmapper.EasyVariantMapper(uta)
hgvs_parser = hgvs.parser.Parser()


def _seqvar_map_func(in_type, out_type):
    func_name = '%s_to_%s' % (in_type, out_type)
    return getattr(mapper, func_name)


def seqvar_length(seqvar):
    comp = VariantComponents(seqvar)
    return len(comp.posedit)


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


def _seqvar_to_seqvar(seqvar, base_type, new_type, transcript=None, maxlen=None):
    """ Using UTA, translate the input seqvar (SequenceVariant object) into 
    the desired new_type of sequence variant.  If the base_type is 'g', a transcript
    label will be required.

    (Optional) Specify `maxlen` (int) characters to restrict resultant SequenceVariant 
    to a reasonable size. Returns None if str(SequenceVariant) > maxlen.

    :param seqvar: (SequenceVariant)
    :param base_type: (str) single-letter abbrev for variant type to map FROM
    :param new_type: (str) single-letter abbrev for variant type to map TO
    :param transcript: (str) [default: None]
    :param maxlen: (int) max length of resultant str(SequenceVariant) to allow
    :return: SequenceVariant or None
    """

    if base_type == new_type:
        return None

    if base_type == 'p':
        return None

    if new_type == 'p' and base_type == 'g':
        return None

    map_seqvar = _seqvar_map_func(base_type, new_type)

    result_seqvar = None

    try:
        if base_type == 'g':
            if transcript:
                result_seqvar = map_seqvar(seqvar, transcript)
            else:
                return None
        result_seqvar = map_seqvar(seqvar)
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

    # if sequence variant maps out to longer than maxlen chars, return None
    try:
        if maxlen and seqvar_length(result_seqvar) > maxlen:
            return None
    except RejectedSeqVar:
        # generally, a "rejected" seqvar is a protein with an empty edit, which is OK to keep.
        pass

    return result_seqvar


class VariantLVG(object):

    #TODO: remove these two class variables (when ready)...
    VERSION = '0.0.2'
    LVG_MODE = 'lvg'

    def __init__(self, hgvs_text_or_seqvar, **kwargs):
        """ Creates a VariantLVG object given an HGVS string or SequenceVariant object.

        Enrichment:

            VariantLVG provides for "enrichment" of lexical variant generation by allowing
            more transcripts and variations to be supplied at instantiation. Just use the 
            appropriate keyword for the type of information, remembering that the "enrichment"
            keyword arguments are all lists.

        Keywords:
            hgvs_c (list): see Enrichment above
            hgvs_g (list): ""
            hgvs_n (list): ""
            hgvs_p (list): ""
            gene_name: accepts any string as gene name (should be HGNC)
            transcripts (list): list of strings describing valid alternative transcripts for seqvar
            seqvar_max_len (int): restrict posedit lengths to this number of characters (or fewer).    
        """

        self.hgvs_text = strip_gene_name_from_hgvs_text('%s' % hgvs_text_or_seqvar)
        self.seqvar = self.parse(hgvs_text_or_seqvar)

        self._gene_name = kwargs.get('gene_name', None)

        seqvar_max_len = kwargs.get('seqvar_max_len', None)


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
                    new_seqvar = _seqvar_to_seqvar(var_c, 'c', this_type, maxlen=seqvar_max_len)
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
                var_c = _seqvar_to_seqvar(self.seqvar, 'g', 'c', trans, maxlen=seqvar_max_len)
                if var_c:
                    self.variants['c'][str(var_c)] = var_c
                var_n = _seqvar_to_seqvar(self.seqvar, 'g', 'n', trans, maxlen=seqvar_max_len)
                if var_n:
                    self.variants['n'][str(var_n)] = var_n

        # With a list of transcripts, we can do g_to_c and g_to_n
        if self.transcripts:
            for trans in self.transcripts:
                for var_g in list(self.variants['g'].values()):
                    # Find all available 'c'
                    if not trans.startswith('NR'):
                        var_c = _seqvar_to_seqvar(var_g, 'g', 'c', trans, maxlen=seqvar_max_len)
                        if var_c:
                            self.variants['c'][str(var_c)] = var_c

                    # Find all available 'n'
                    var_n = _seqvar_to_seqvar(var_g, 'g', 'n', trans, maxlen=seqvar_max_len)
                    if var_n:
                        self.variants['n'][str(var_n)] = var_n

        # map all newly found 'c' to 'p'
        for var_c in list(self.variants['c'].values()):
            new_seqvar = _seqvar_to_seqvar(var_c, 'c', 'p', maxlen=seqvar_max_len)
            if new_seqvar:
                self.variants['p'][str(new_seqvar)] = new_seqvar

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

        hgvs_text = strip_gene_name_from_hgvs_text(hgvs_text_or_seqvar)

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

    def _simple_dict(self):
        return {'gene_name': self.gene_name,
                'hgvs_c': self.hgvs_c,
                'hgvs_g': self.hgvs_g,
                'hgvs_n': self.hgvs_n,
                'hgvs_p': self.hgvs_p,
                'hgvs_text': self.hgvs_text,
                'transcripts': list(self.transcripts),
               }

    def to_json(self):
        """ Returns a JSON representation of this object, one that can be used to 
        instantiate this object again using the "from_json" class method.

        :return: (str) json repr
        """
        return json.dumps(self._simple_dict())

    @classmethod
    def from_json(cls, json_str):
        """ Allows instantiation of VariantLVG object from pre-established values 
        set in a JSON data structure.  The structure is as follows:
                
                {
                'hgvs_text': <hgvs_text>,
                'gene_name': <gene_name>,
                'hgvs_c': [hgvs_c1,hgvs_c2...],
                'hgvs_g': [hgvs_g1,hgvs_g2...],
                'hgvs_n': [hgvs_n1,hgvs_n2...],
                'hgvs_p': [hgvs_p1,hgvs_p2...],
                'transcripts': [<transcript1>,<transcript2>...]
               }

        At least hgvs_text must be set.
        """
        inpd = json.loads(json_str)
        hgvs_text = inpd.pop('hgvs_text')

        return cls(hgvs_text, **inpd)

    def __str__(self):
        out = 'HGVS input: %s\n' % self.hgvs_text
        out += '%r' % self.seqvar
        return out


### API Convenience Functions

Variant = VariantLVG.parse

