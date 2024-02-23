SH=/bin/bash
.PHONY=fortran


signature:
	f2py3 -m mfeq -h mfeq.pyf --overwrite-signature eq.f90


fortran:
	f2py3 -llapack -c -m mfeq mfeq.pyf eq.f90

# clibeq: gcc -Wall -shared eq.c
