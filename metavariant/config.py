from __future__ import absolute_import, unicode_literals

import os, socket
from configparser import ConfigParser

import hgvs.dataproviders.uta

PKGNAME = 'metavariant'

DEBUG = bool(os.getenv('%s_DEBUG' % PKGNAME, False))
ENV = os.getenv('%s_ENV' % PKGNAME, 'dev')

UTA_HOST = os.getenv('UTA_HOST', 'default')
UTA_PORT = os.getenv('UTA_PORT', 5432)
UTA_SCHEMA = os.getenv('UTA_SCHEMA', 'uta_20150903')
UTA_TIMEOUT = os.getenv('UTA_TIMEOUT', 3)

####
import logging
log = logging.getLogger(PKGNAME)
if DEBUG:
    log.setLevel(logging.DEBUG)
else:
    log.setLevel(logging.INFO)
####

#log.debug('%s config dir: %s' % (PKGNAME, CFGDIR))
#log.debug('%s env: %s' % (PKGNAME, ENV))

def get_uta_connection(host=UTA_HOST, port=UTA_PORT, timeout=UTA_TIMEOUT, schema=UTA_SCHEMA, 
                        username='uta_admin', password='anonymous'):
    """ Returns an open connection to a UTA host at given coordinates, if possible.
    
    If host=='default', returns connection to default UTA server (uta.biocommons.org).

    :return: open UTA connection
    :raises: socket, uta, hgvs Exceptions
    """

    timeout = int(timeout)
    port = int(port)

    uta_cnxn_tmpl = 'postgresql://{user}:{pwd}@{host}:{port}/uta/{schema}/'

    if host == 'default':
        return hgvs.dataproviders.uta.connect()

    socket.create_connection((host, port), timeout=timeout)
    return hgvs.dataproviders.uta.connect(uta_cnxn_tmpl.format(host=host, port=port, schema=schema, user=username, pwd=password), pooling=True)

