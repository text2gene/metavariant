"""Handles env variable configuration usage and sets up UTA connection."""

import os, socket
from configparser import ConfigParser

import hgvs.dataproviders.uta

PKGNAME = 'metavariant'

DEBUG = bool(os.getenv('%s_DEBUG' % PKGNAME, False))
ENV = os.getenv('%s_ENV' % PKGNAME, 'dev')

UTA_HOST = os.getenv('UTA_HOST', 'default')
UTA_PORT = os.getenv('UTA_PORT', 5432)
UTA_SCHEMA = os.getenv('UTA_SCHEMA', 'uta_20171026')
UTA_TIMEOUT = os.getenv('UTA_TIMEOUT', 3)
UTA_USER = os.getenv('UTA_USER', 'uta_admin')
UTA_PASS = os.getenv('UTA_PASS', 'uta_admin')

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
                        username=UTA_USER, password=UTA_PASS):
    """ Returns an open connection to a UTA host at given coordinates, if possible.
    
    If host=='default', returns connection to default UTA server (uta.biocommons.org).

    :return: open UTA connection
    :raises: socket, uta, hgvs Exceptions
    """

    timeout = int(timeout)
    port = int(port)

    uta_cnxn_tmpl = 'postgresql://{user}:{pwd}@{host}:{port}/uta/{schema}/'
    cnxn_desc = uta_cnxn_tmpl.format(host=host, port=port, schema=schema, user=username, pwd=password)

    if host == 'default':
        return hgvs.dataproviders.uta.connect()

    socket.create_connection((host, port), timeout=timeout)
    try:
        return hgvs.dataproviders.uta.connect(cnxn_desc, pooling=True)
    except:
        raise Exception('Cannot connect to UTA at {}'.format(cnxn_desc))

