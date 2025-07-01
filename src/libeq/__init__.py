from .jacobian import amatrix, dlogcdlogbeta


__version__ = '0.11'

# The engine to be used for calculations
# Accepted values are 'python' or 'fortran'
# LIBEQ_ENGINE = 'python'
LIBEQ_ENGINE = 'fortran'
