******************************************************
MetaVariant
******************************************************

Python toolkit for parsing HGVS strings describing genetic variants (mutations)
into components and generating lexical synonyms of these descriptors.

Requirements
============

You'll need a UNIX-like environment to use this package. Both OS X and Linux have been confirmed to work.

Additionally, your system must contain the C development version of the client library
for Postgres 9.3 or later (needed to connect to UTA services). 

On OS X using brew:

.. code-block:: bash

  $ brew install postgresql

On Linux using Ubuntu:

.. code-block:: bash

  $ apt-get install libpq-dev

Python library dependencies for variant_lexicon are the following::

  numpy
  biopython
  uta
  hgvs

The above "should" be successfully installed by using setup.py and/or by installing
this package via pip.


Optional but Recommended: IPython
---------------------------------

Install ipython via pypi:

.. code-block:: bash

  $ pip install ipython
  
More installation options and instructions are available on `the iPython installation page <http://ipython.org/ipython-doc/stable/install/install.html>`_.

Variant: shortcut to SequenceVariant
====================================

`Variant` is an API shortcut to instantiating a SequenceVariant object from an HGVS string. It accepts
both gene-name-bearing and non-gene-name-bearing HGVS strings. 

See `Human Genome Variation Society <http://www.hgvs.org/>`_ for explanation of HGVS strings.

Examples of valid HGVS strings::

  NM_198056.2:c.4786T>A
  NC_000023.10:g.32841489delT
  NP_009225.1:p.Glu1000_Glu1001del
  NM_001165963.1(SCN1A):c.384_1662dup
  

Note: some valid HGVS strings currently cannot be parsed by the hgvs library (as of 6/17/2016). Examples::

  NC_000001.11:g.(?_215888426)_(215900905_?)del
  NM_000251.2:c.212-?_366+?dup

When an illogical or invalid HGVS string is provided to `Variant`, this function returns None. (Tune into 
the "metavariant" logger at INFO level for information.)

`Variant` will also accept a SequenceVariant object as its argument, returning the same variant (identity).


VariantLVG: Lexical Variant Generation
======================================

`VariantLVG` takes an hgvs.SequenceVariant object or a plain HGVS string and uses the Universal Transcript
Archive (UTA) to find as many related transcripts and sequence variants as possible that could be used as
"synonyms" of the provided variant.

This class raises `CriticalHgvsError` upon instantiation if the supplied HGVS string fails to create a 
valid SequenceVariant object.

Usage example:

.. code-block:: python

  hgvs_text = 'NM_198056.2:c.4786T>A'
  lex = VariantLVG(hgvs_text)

  print(lex.variants)
  print(lex.transcripts)
  print(lex.gene_name)

Enrichment
----------

`VariantLVG` provides for "enrichment" of lexical variant generation by allowing
more transcripts and variations to be supplied at instantiation. Just use the 
appropriate keyword for the type of information, remembering that the "enrichment"
keyword arguments are all lists.

Keyword Arguments
-----------------

- transcripts (list): list of strings describing valid alternative transcripts for seqvar
- seqvar_max_len (int): restrict posedit lengths to this number of characters (or fewer).  
- hgvs_c (list): see Enrichment above
- hgvs_g (list): see Enrichment above
- hgvs_n (list): see Enrichment above
- hgvs_p (list): see Enrichment above

Attributes
----------

- hgvs_text: original hgvs string from instantiation
- seqvar: original SequenceVariant from instantiation
- transcripts: list of strings indicating related transcripts
- variants: 2-level dictionary of shape { seqtype: { hgvs_text: seqvar } }

Properties
----------

+ gene_name: returns HUGO gene-name if it can be ascertained using UTA. (Lazy-loaded attribute _gene_name.)
+ hgvs_c: returns flat list of c.DNA hgvs strings from variants
+ hgvs_g: returns flat list of g.DNA hgvs strings from variants
+ hgvs_n: returns flat list of n.RNA hgvs strings from variants 
+ hgvs_p: returns flat list of protein hgvs strings from variants
+ seqvars: returns flat list of SequenceVariant objects from variants

Methods
-------

+ to_dict(): returns non-underscored attributes (seqvar, hgvs_text, transcripts, seqvars) as dictionary
+ to_json(): returns a serialized JSON string representation of the object which can be used to instantiate this LVG again.
+ from_json(json_str): takes serialized JSON representation of this object and rebuilds LVG from its details.

VariantComponents: Parsing and "Slang"
======================================

`VariantComponents` allows instantiation in two different ways: using an hgvs.SequenceVariant object
and using a set of components as keyword arguments.

Usage starting from a SequenceVariant object::

    comp = VariantComponents(seqvar)

Usage starting from individual components::

    comp = VariantComponents(seqtype='c', edittype='SUB', pos='322', ref='C', alt='T')

Usage starting from "aminochange" string::

    comp = VariantComponents(aminochange='V777A')
    comp = VariantComponents(aminochange='Leu653Gly')

(If starting with "aminochange" string, the `seqvar` and `edittype` attributes will be None.)

If no seqtype is supplied, VariantComponents tries to infer the sequence type heuristically (e.g. the presence
of a "U" in the ref or the alt implies this is an RNA string).

VariantComponents may raise a `RejectedSeqVar` exception during instantiation (see "Exceptions" below).

A VariantComponents object provides access to the following attributes and properties::
    
   seqtype: the sequence type of this seqvar (one of 'c', 'g', 'g', 'n')
   edittype: the type of mutation represented by this variant ('SUB', 'DEL', 'FS', etc)
   pos: position of the edit
   ref: reference sequence at given position (aka "wildtype")
   alt: alternate (or "mutant") at given position

   posedit: returns the HGVS "official" construction of this seqvar's position + edit information.
   posedit_slang: returns a list of algorithmically generated "slang" for given seqvar's posedit.


Exceptions
==========

All exceptions can be found and imported from metavariant.exceptions.

`CriticalHgvsError`: raised when input HGVS string fails to instantiate a SequenceVariant object within the VariantLVG object.

`RejectedSeqVar`: raised inside VariantComponents when input sequence components fail certain tests of completeness. For example, a protein seqvar will throw this Exception if the protein effect string is only a "?" (i.e. unknown protein effect).  A "SUB" (substitution) will fail the completeness test if an "alt" is not provided in the instantiated components.

`MetaVariantException`: base exception class from which all metavariant package exceptions are subclassed.

`NCBIRemoteError`: raised when NCBI doesn't respond well to a request. Special characters or improper data will do this.

Setting UTA Server
==================

When you find yourself outgrowing the public UTA server, you may want to install your own UTA server. 
The instructions can be found on the [biocommons/uta repository README](https://bitbucket.org/biocommons/uta). Both "installing from database dumps" 
and installing from docker have been tested working with metavariant.

*How to Set UTA Host Variables*

metavariant connects to a UTA server as soon as it is imported, so reconfiguring which UTA server is used
happens at the environment variable level.  The relevant variables::

    UTA_HOST (default: 'default')
    UTA_PORT (default: 5432)
    UTA_USER (default: 'uta_admin')
    UTA_PASS (default: 'uta_admin')
    UTA_SCHEMA (default: 'uta_20150903')
    UTA_TIMEOUT (default: 3) -- how long to wait before giving up on a connection

When you set up your own postgres server for the UTA database and you connect on the same server, the only
environment variable you probably need to change is `UTA_HOST` (set it to "localhost").

You may have to do more postgres administration to get your preferred configuration going, which is outside the scope of this README.


Support and Maintenance
=======================

This library was developed in-house at Text2Gene, LLC.

It is provided to the community free of charge by way of the Apache 2.0 License.

You are free to modify it for commercial and non-commercial uses; just don't try to sell it as-is.

Contributions, extensions, pull requests, issues, suggestions, and swear words all happily accepted, 
in that order.

naomi@text2gene.com
2016 and onwards

