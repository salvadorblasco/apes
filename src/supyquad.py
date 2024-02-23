#!/usr/bin/python3
# -*- encoding: utf-8 -*-

"""
This script reads an input from the stardard input and writes results
to standard output. Standard input data must match a superquad file
format.

usage: SUPYQUAD [-h] [--max-iterations N] [--verbosity VERBOSITY] [--version]
                [--leastsquares | --simplex]

A python implementation of the SUPERQUAD algoritm

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
"""

import sys
import platform
import argparse

import numpy as np

import consts
import libio
import libemf
import libeq
import libaux
import libmath
import report

VERSION = "0.3"
description = """A python implementation of the SUPERQUAD algoritm"""

parser = argparse.ArgumentParser(prog="SUPYQUAD",
                                 description=description)
parser.add_argument("--max-iterations", metavar='N', type=int,
                    help="Sets the maximum number of iterations")
parser.add_argument("--threshold", metavar='T', type=float,
                    help="Sets the convergence criteria")
parser.add_argument("--verbosity", "-v", type=int, choices=(0, 1, 2),
                    help="Sets the level of verbosity for the output")
parser.add_argument('--version', action='version', version='%(prog)s 2.0')

method = parser.add_mutually_exclusive_group()
method.add_argument("--leastsquares", "-lm", "-ls",
                    help="""sets the Levenberg-Marquardt algorithm (least
                         squares) as the algorithm to be used for fitting.""",
                    action="store_true", default=True)
method.add_argument("--simplex", "-nm",
                    help="""sets the Nelder-Mead algorithm (simplex) as the
                         algorithm to be used for fitting.""",
                    action="store_true")
parser.add_argument("-f", metavar='FILE',
                    help="""The data file to read. If not provided the
                         standard input is used""", default=False)

args = parser.parse_args()

print("SUPYQUAD ", VERSION, '\n')
print("Environment:")
print("\tplatform ", platform.platform(), platform.processor())
print("\tpython ", sys.version)
print("\tnumpy ", np.version.version)
print("\tlibeq ", libeq.__version__, " engine: ", libeq.LIBEQ_ENGINE)
print("\tlibemf ", libemf.__version__)
print("\tlibaux ", libaux.__version__)
print("\tlibio ", libio.__version__)

# if args.f:
#     d = libio.importSuperquad(args.f)
# else:
#     d = libio.importSuperquad(sys.stdin)

source = args.f if args.f else sys.stdin
d = libio.importSuperquad(source)

method = 1 if args.simplex else 0
freport = report.spyq_nm_iteration if args.simplex else report.spyq_lm_iteration

title = next(d)
control_numbers = next(d)
labels = next(d)
temperature = next(d)
log10B = next(d)
P = next(d)
Bflags = next(d)

titrations = []
electrodes = []
for amounts, electr, titdata in d:
    titrations.append({
        'T0': amounts[2],
        'buret': amounts[3],
        'Tflags': amounts[4],
        'V': titdata[0],
        'V0': electr[0],
        'error_V': electr[1],
        'emf': np.array(titdata[1])})

    electrodes.append({
        'E0': electr[4],
        'hindex': electr[3] - 1,
        'fRTnF': consts.RoverF*(temperature+273.15)/electr[2],
        'error_emf': electr[5]})

E = len(P)
S = len(P[0])
n_curves = len(titrations)

print("\n#### SUPYQUAD ####")
print("\nINPUT:")
print("TITLE: ", title)

weights = [np.fromiter(libmath.weighting_slope(d1['V'], d1['emf'],
                                               d1['error_V'],
                                               d2['error_emf']),
                       dtype=np.float)
           for d1, d2 in zip(titrations, electrodes)]

# for d1, d2 in zip(titrations, electrodes):
#     print(d1['V'], d1['emf'], d1['error_V'], d2['error_emf'])
# sys.exit(0)

for curve, weight in zip(titrations, weights):
    curve['weights'] = weight

report.spyq_start(log10B, Bflags, P, titrations, weights, verbosity=1)

newB, newC, alt = libemf.emffit(np.array(log10B), Bflags, np.array(P),
                                titrations, electrodes,
                                report=freport,
                                method=method)

report.spyq_finalconvg(newB, errors=alt['error_beta'], titrations=titrations,
                       residuals=alt['residuals'],
                       correlation=alt['correlation'])

print("PROGRAM FINISHED")

# _emf_ = [np.array(t['emf']) for t in titrations]
# _v_ = [np.array(t['V']) for t in titrations]
# np.savez('plot_emf.npz', free_concentration=newC,  emf=_emf_, titre=_v_, **alt)
