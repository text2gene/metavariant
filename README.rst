******************************************************
Variant-Lexicon
******************************************************

Python toolkit for parsing HGVS strings describing genetic variants (mutations)
into components and generating lexical synonyms of these descriptors.

Requirements
------------

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


HgvsLVG: Lexical Variant Generation
-----------------------------------

(Discussion)

Let's say you've been handed an HGVS string describing a seen variant. Examples::

  NM_
  NC_
  NP_

*VariantLVG*


HgvsComponents: Parsing and "Slang"
-----------------------------------

(Discussion)

Let's say you have discovered in clinical literature a some text referring to mutations,
described in a few different ways (aka "variant slang"), for example::

  410 T>C
  Arg459Pro

*VariantComponents*


Support and Maintenance
-----------------------

This library was developed in-house at Text2Gene, LLC.

It is provided to the community free of charge by way of the Apache 2.0 License.

You are free to modify it for commercial and non-commercial uses; just don't try to sell it as-is.

Contributions, extensions, bug reports, suggestions, and swear words all happily accepted, 
in that order.

naomi@nthmost.com
2016 and onwards

