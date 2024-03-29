import glob, os

from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext as _build_ext


class build_ext(_build_ext):
    # http://stackoverflow.com/questions/21605927/why-doesnt-setup-requires-work-properly-for-numpy/21621493
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

setup (
    name = 'metavariant',
    version = '0.3.0',
    description = 'a lexical manipulation toolkit for genetic variant descriptors (hgvs, etc)',
    long_description = open('README.rst').read(),
    long_description_content_type = 'text/x-rst',
    url = 'https://github.com/text2gene/metavariant',
    author = 'Naomi Most',
    maintainer = 'Naomi Most',
    author_email = 'naomi@nthmost.com',
    maintainer_email = 'naomi@text2gene.com',
    license = 'Apache 2.0',
    packages = find_packages(),
    cmdclass = {'build_ext': build_ext},
    setup_requires = ['setuptools', 'numpy'],
    install_requires = [
        'numpy',
        #'uta',     # let this be installed by hgvs
        'hgvs>=1.3.0.post0',
        'xmltodict',
        'pytest',
        ],
    )

