
try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources

import json    
    
from .mollusc import *
from . import mollusc_data


sp_pars = {par['ref']: par for par in         
    json.loads(pkg_resources.read_text(mollusc_data, 'sp_pars.json'))}


class MolluscExample(Mollusc, object):
    
    def __init__(self, ref):
        super().__init__(**sp_pars[ref])
        
