from __future__ import print_function, unicode_literals

import logging

from .config import PKGNAME
from .exceptions import RejectedSeqVar

log = logging.getLogger(PKGNAME)

def lowercase_all_the_keys(some_dict):
    """ lowercases all the keys in the supplied dictionary.

    :param: dict
    :return: dict
    """
    return dict((key.lower(), val) for key, val in some_dict.iteritems())

amino_acid_map = { 'Ala': 'A',
                   'Arg': 'R',
                   'Asn': 'N',
                   'Asp': 'D',
                   'Cys': 'C',
                   'Glu': 'E',
                   'Gln': 'Q',
                   'Gly': 'G',
                   'His': 'H',
                   'Ile': 'I',
                   'Leu': 'L',
                   'Lys': 'K',
                   'Met': 'M',
                   'Phe': 'F',
                   'Pro': 'P',
                   'Ser': 'S',
                   'Thr': 'T',
                   'Trp': 'W',
                   'Tyr': 'Y',
                   'Val': 'V',
                 }

dna_nucleotides = ['A','C','T','G']
#rna_nucleotides = ['A','C','U','G']

official_to_slang_map = {'>': ['->', '-->', '/'],
                         'Ter': ['*', 'X'],
                        }

class VariantComponents(object):
    """
    VariantComponents

    Takes an hgvs.SequenceVariant object and parses it into components (position, reference,
    alternate, edit type, and sequence type), allowing further lexical manipulation for the
    purposes of generating pubmed evidence source "slang".

    Can also be instantiated from keyword components, helping organize this information into 
    an object, which will also perform seqtype inference (if necessary) and official posedit
    syntax construction.

    Note: FrameShift (FS) and Duplication (DUP) are currently not well supported.

    Attributes:

        seqtype: the sequence type of this seqvar (one of 'c', 'g', 'g', 'n')
        edittype: the type of mutation represented by this variant ('SUB', 'DEL', 'FS', etc)
        pos: position of the edit
        ref: reference sequence at given position
        alt: alternate (or "wildtype") at given position

    Properties:

        posedit: returns the HGVS "official" construction of this seqvar's position + edit information.
        posedit_slang: returns a list of algorithmically generated "slang" for given seqvar's posedit.

    Usage:

        # from SequenceVariant object:
        comp = VariantComponents(seqvar)

        # from individual components:
        comp = VariantComponents(seqtype='c', edittype='SUB', pos='128', ref='C', alt='T')
    """

    def __init__(self, seqvar=None, **kwargs):
        kwargs = lowercase_all_the_keys(kwargs)

        self.seqvar = seqvar
        if self.seqvar:
            self.seqtype, self.edittype, self.ref, self.pos, self.alt = self.parse(seqvar)
            #TODO: get FS_Pos and DupX out of seqvar when applicable.
        else:
            # names of keywords match capitalization used in MySQL m2p_* tables in pubtatordb
            self.seqtype = kwargs.get('seqtype', '').strip()
            self.edittype = kwargs.get('edittype', '').strip()
            self.ref = kwargs.get('ref', '').strip()
            self.pos = kwargs.get('pos', '').strip()
            self.alt = kwargs.get('alt', '').strip()
            self.fs_pos = kwargs.get('fs_pos', '').strip()
            self.dupx = kwargs.get('dupx', '').strip()

        if self.edittype.upper() == 'DELINS':
            # normalize DELINS to INDEL (synonymous)
            self.edittype = 'INDEL'

        if not self.seqtype:
            self.seqtype = self._infer_seqtype()

        #TODO: (attempt to) construct a seqvar if instantiated from kwargs.

    def _infer_seqtype(self):
        # If SeqType is none and REF in [u] or ALT in [u] --> then RNA
        # If SeqType is none and REF in [AminoAcidsList] and ALT in [AminoAcidsList] --> then Protein
        # If SeqType is none and REF in [a,c,t,g] and ALT in [a,c,t,g] --> then DNA or RNA
        
        refalt = self.ref.upper() + self.alt.upper()

        if 'u' in refalt:
            # Definitely RNA: there's no "U" amino acid and no "U" in the DNA nucleotides.
            return 'n'

        for char in refalt:
            if char not in dna_nucleotides:
                if char in list(amino_acid_map.values()):
                    return 'p'

        # it's "probably" DNA, but we can't know for sure.
        return ''

    @staticmethod
    def parse(seqvar):
        """ return tuple of sequence variant components as
        (seqtype, edittype, ref, pos, alt)
        """
        if seqvar.type.strip() == '':
            log.warn('SequenceVariant has empty seqtype. (%r)' % seqvar)
            seqtype = ''
        else:
            seqtype = seqvar.type

        ref = alt = edittype = pos = ''

        #if seqvar.posedit.edit in ['?', '=']:
        if hasattr(seqvar.posedit.edit, 'lower'):
            # if the edit attribute is a str or unicode type thing:
            raise RejectedSeqVar('SequenceVariant missing edit information. (%r)' % seqvar)

        edittype = seqvar.posedit.edit.type.upper()

        try:
            ref = seqvar.posedit.edit.ref
            alt = seqvar.posedit.edit.alt
        except AttributeError as error:
            # We expect a 'dup' variant not to have an alt. warn us if this is a different situation.
            if not seqvar.posedit.edit.type == 'dup':
                log.warn('SequenceVariant %s: %s', seqvar, error)

        # Termination (STOP) codon: normalize '*' to 'X'
        if alt == '*':
            alt = 'X'

        if seqtype == 'p':
            try:
                pos = '%s' % seqvar.posedit.pos.start.pos
                ref = '%s' % seqvar.posedit.pos.start.aa
            except AttributeError:
                raise RejectedSeqVar('Protein entry incomplete (unusable).')

        else:
            if seqvar.posedit.pos.end != seqvar.posedit.pos.start:
                # compose a "range"
                pos = '%s_%s' % (seqvar.posedit.pos.start, seqvar.posedit.pos.end)
            else:
                pos = '%s' % seqvar.posedit.pos.start

        return seqtype, edittype, ref, pos, alt

    def to_mysql_dict(self):
        outd = { 'Ref': self.ref,
                 'Alt': self.alt,
                 'SeqType': self.seqtype,
                 'EditType': self.edittype,
                 'Pos': self.pos,
                 }

        if self.edittype == 'FS':
            outd['FS_Pos'] = self.fs_pos

        if self.edittype == 'DUP':
            outd['DupX'] = self.dupx

        return outd

    @property
    def posedit(self):
        """ Returns the official lexeme representing this variant's position and edit information. """
        return '%s' % self.seqvar.posedit

    def _posedit_slang_protein(self):
        out = set()
        posedit = self.posedit.replace('(', '').replace(')', '')
        for item in official_to_slang_map['Ter']:
            out.add(posedit.replace('Ter', item))

        fs_pos = posedit.find('fs')
        if fs_pos > -1:
            out.add(posedit[:fs_pos + len('fs')])
        else:
            # e.g. Lys2569Gly produces "K2569G"
            out.add('%s%s%s' % (self.ref, self.pos, self.alt))
        return list(out)

    def _posedit_slang_SUB(self):
        """ Handles the Substitution case for generating posedit slang from Components. """
        # Examples based on input hgvs_text 'NM_014874.3:c.891C>T'
        out = set()

        # E.g. ["891C>T", "891C->T", "891C-->T", "891C/T"]
        for slang_symbol in official_to_slang_map['>']:
            out.add(self.posedit.replace('>', slang_symbol))

        # E.g. "C891T"
        out.add(self.ref + self.pos + self.alt)
        return list(out)

    def _posedit_slang_DEL(self):
        """ Handles the Deletion case for generating posedit slang from Components. """
        # Examples based on input hgvs_text 'NM_007294.3:c.4964_4982delCTGGCCTGACCCCAGAAGA'
        #
        # E.g. '4964_4982del'
        return [self.pos + 'del']

    def _posedit_slang_DUP(self):
        """ Handles the Duplication case for generating posedit slang from Components. """
        # Examples based on input hgvs_text 'NM_025114.3:c.6869dupA'
        #
        # E.g. '6869dup'
        return [self.pos + 'dup']

    def _posedit_slang_INDEL(self):
        """ Handles the Insertion case for generating posedit slang from Components.

        Returns an empty list (no slang terms known).
        """
        return []

    def _posedit_slang_INS(self):
        """ Handles the Insertion case for generating posedit slang from Components.

        Returns an empty list (no slang terms known).
        """
        return []

    @property
    def posedit_slang(self):
        """ If supported, returns alternative lexeme that may represent this variant's position and edit info in the wild. """
        if self.seqtype == 'p':
            return self._posedit_slang_protein()

        try:
            slang_method = getattr(self, '_posedit_slang_%s' % self.edittype)
            return slang_method()
        except AttributeError:
            raise NotImplementedError('Cannot currently handle EditType %s' % self.edittype)

    def to_dict(self):
        return self.__dict__

    def __str__(self):
        return '%r' % self.__dict__

    def __repr__(self):
        return '%r' % self.__dict__

