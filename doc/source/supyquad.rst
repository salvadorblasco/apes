SUPYQUAD script
===============

SUPYQUAD is a python implementation of SUPERQUAD.[1]_ It aims to use the same
data format to fit equilibrium constants from potentiometric data.

Usage
-----

``$ supyquad.py [-h] [--max-iterations N] [--verbosity VERBOSITY] [--version]
              [--leastsquares | --simplex] < input > output``

optional arguments:
  -h, --help            show this help message and exit
  --max-iterations N    Sets the maximum number of iterations
  --verbosity VERBOSITY, -v VERBOSITY
                        increase output verbosity
  --version             show program's version number and exit
  --leastsquares, -lm, -ls
                        sets the Levenberg-Marquardt algorithm (least squares)
                        as the algorithm to be used for fitting.
  --simplex, -nm        sets the Nelder-Mead algorithm (simplex) as the
                        algorithm to be used for fitting.

Where ``input`` is an input with complies with the input format for SUPERQUAD. 
The input is processed and the output is sent to the standard output. 

.. [1] Peter Gans, Antonio Sabatini, Alberto Vacca, *J. Chem. Soc. Dalton Trans.* 1985, 1195-1200
