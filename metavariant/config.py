from __future__ import absolute_import, unicode_literals

import os, socket
from configparser import ConfigParser

import hgvs.dataproviders.uta

PKGNAME = 'metavariant'


default_config_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config')

CFGDIR = os.getenv('METAVARIANT_CONFIG_DIR', default_config_dir)

DEBUG = bool(os.getenv('%s_DEBUG' % PKGNAME, False))
ENV = os.getenv('%s_ENV' % PKGNAME, 'dev')

UTA_HOST = os.getenv('UTA_HOST', 'default')
UTA_PORT = os.getenv('UTA_PORT', 'default')

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

configs = [os.path.join(CFGDIR, x) for x in os.listdir(CFGDIR) if x.find(ENV+'.ini') > -1]

CONFIG = ConfigParser()
CONFIG.read(configs)

def get_uta_connection(host='default', port='default'):
    """ Returns an open UTA connection to the first available UTA host, as prioritized
    in the config file (see uta_* items under the [hgvs] section).

    If no connection can be made, returns None.

    Also sets UTA_HOST and UTA_PORT globals in this file.

    :return: open UTA connection or None if all connections fail
    """

    global UTA_HOST
    global UTA_PORT

    UTA_HOST = host
    UTA_PORT = port

    timeout = int(CONFIG.get('hgvs', 'uta_connection_timeout'))
    uta_cnxn_tmpl = CONFIG.get('hgvs', 'uta_connection_tmpl')

    for pair in host_port_pairs:
        host = pair['host']
        port = int(pair['port'])
        if host == 'default':
            return hgvs.dataproviders.uta.connect()
        else:
            try:
                socket.create_connection((host, port), timeout=timeout)
                log.info('Connected to UTA host %s on port %i', host, port)
                UTA_HOST = host
                UTA_PORT = port
                return hgvs.dataproviders.uta.connect(uta_cnxn_tmpl.format(host=host), pooling=True)
            except Exception as error:
                log.info('Could not connect to UTA host %s on port %i: %r', host, port, error)

    # if we get this far and nothing can be reached, return the biocommons default UTA host
    return hgvs.dataproviders.uta.connect()

