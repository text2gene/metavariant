"""Provides metavariant specific exceptions."""

class MetaVariantException(Exception):
    """ Generic exception class from which all metavariant exceptions are subclassed. """


class RejectedSeqVar(MetaVariantException):
    """ Indicates a sequence variant that cannot be considered functionally
    useful for the purposes of lookup in databases (pubtator, etc).

    Message should contain information about why the SequenceVariant fails.
    """
    pass


class CriticalHgvsError(MetaVariantException):
    """ Raised when a SequenceVariant cannot be parsed from the input hgvs_text; should be used only in cases
        where failure to create a SequenceVariant means a complete (critical) failure to proceed with LVG.

        Message should contain information passed along from the HgvsParseError.

        e.g. "HGVSParseError(u'NM_004628.4:c.621_622ins83: char 24: Syntax error',)"
    """
    pass


class NCBIRemoteError(MetaVariantException):
    """ Raised when NCBI fails to report on a variant for any reason (often mysterious/arbitrary ones). """
    pass


class LOVDRemoteError(MetaVariantException):
    """ Raised when response from LOVD contains a not-ok status code.
    Error message should contain status code. 
    """
    pass




