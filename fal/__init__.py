try:
    from ._version import __version__
except(ImportError):
    pass

from . import fitting
from . import prep
from . import utils
