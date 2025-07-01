__version__ = '0.11'

# The engine to be used for calculations
# Accepted values are 'python' or 'fortran'
# LIBEQ_ENGINE = 'python'
LIBEQ_ENGINE = 'fortran'

from .consol import initial_guess
# from .consol import consol as solve_concentration
from .jacobian import amatrix, dlogcdlogbeta
