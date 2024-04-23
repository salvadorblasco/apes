"""
Common parameters and constants used applicationwide.

.. module:: consts.py
.. moduleauthor:: Salvador Blasco <salvador.blasco@gmail.com>
"""

# from collections import namedtuple
# Model = namedtuple('Model', ('name', 'const', 'stoich', 'const_flags',
#                              'const_error', 'enthalpy', 'enthalpy_error',
#                              'enthalpy_flags'))

LOGK = 2.3025851      # ln(10) = 1/log(e)
R = 8.314472          # gas constant(J ∕ (K⋅mol))
fRTnF = 25.6926       # mV
NERNST = 25.6926      # mV
RoverF = 0.086173424  # mV/K

FM_NONE, FM_SIM, FM_EMF, FM_CALOR, FM_SPEC, FM_NMR = range(6)

# Refinement methods
METHOD_LM = 0   # method Levenberg-Marquardt, least squares
METHOD_NM = 1   # method Nelder-Mead, simplex

# Refinement flags
RF_IGNORE = -1
RF_CONSTANT = 0
RF_REFINE = 1
RF_CONSTRAINT1 = 2
RF_CONSTRAINT2 = 3
RF_CONSTRAINT3 = 4
RF_CONSTRAINT4 = 5
RF_CONSTRAINT5 = 6
RF_CONSTRAINT6 = 7

REFINE_BLABELS = ("ignore", "constant", "refine", "constraint 1",
                  "constraint 2", "constraint 3", "constraint 4",
                  "constraint 5", "constraint 6")

REFINE_LABELS = REFINE_BLABELS[1:]

WEIGHT_UNIT = 0
WEIGHT_AUTO = 1

OUTPUT_MINIMAL = 0
OUTPUT_REGULAR = 1
OUTPUT_VERBOSE = 2
