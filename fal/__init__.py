try:
    from ._version import __version__
except(ImportError):
    pass

__abspath__ = '/Users/pcargile/Astro/gitrepos/FALv2/'

from . import fitting
from . import prep
from . import utils
