#!/usr/bin/python
# -*- encoding: utf-8 -*-

import sys

from pyss import simulation, make_labels
from matplotlib import pyplot
import numpy as np

from libaux import make_T, make_ranges


def find_zeros(T):
    assert isinstance(T, np.ndarray)
    return [ n for n in range(T.shape[1]) if np.all(T[:,n] == 0) ]

labels = sys.stdin.readline().split()
flags = sys.stdin.readline().split()
formats = sys.stdin.readline().split()
assert len(labels) == len(flags)
t = np.array( 
    [ 
        map(float, sys.stdin.readline().split()),
        map(float, sys.stdin.readline().split())])
assert t.shape == (2, len(labels))

pT = map(lambda x: True if ('p' in x) else False, flags)

K_raw = list()
P_raw = list()
for line in sys.stdin.readlines():
    if line.strip() == '':          # consider newline the end
        break
    a = line.split()
    assert len(a) == len(labels) + 2
    K_raw.append(float(a[0]))
    P_raw.append(map(int,a[1:-1]))
    #toplot.append(a[-1] != '0')
    formats.append(a[-1])

toplot = map(
    lambda x: False if x=='0' or 'x' in x else True,
    formats  )

# assertions section
# I. only one 'x'
#for i in flags1:
#    if 'x' in i
#    n += 1

assert sum(1 if 'x' in i else 0 for i in flags) == 1, "There can only be one 'x'"
assert sum(1 if '%' in i else 0 for i in flags) == 1, "There can only be one '%'"

N = 50
P = np.array(P_raw)
E = P.shape[0]
S = P.shape[1]
logB = np.array(K_raw).reshape((E,1))
B = np.power(10, logB)
T = zip(t[0], t[1])
fullT = make_T(T, N)

for n, i in enumerate(flags):
    if '%' in i:
        relative = n
    if 'x' in i:
        x = n

full_labels = make_labels(labels, P)
#print "B = ", B
#print "P = ", P
#print "T = ", T
#print "pT = ", pT
#print "x = ", x

to_remove = find_zeros(t)

H, finalC = simulation(B, P, T, pT, N, x)

if 'p' in flags[x]:
    pH = -np.log10(H)
else:
    pH = H

print "ylabel=\\%% formaci\\'on relativa al %s" % full_labels[relative]
print "ignore_lower_than=5"
# styles
for i, s in enumerate(formats):
    if s != '0':
        print "override_datastyle[%d]=%s" % (i + 3 + S, s)
print
print "x=%d from %f to %f" % (x, T[n][0],T[n][1]) 
print "y=%s" % make_ranges(toplot, S+3)
print "y/$[%d]*n[%d]*100" % (relative+2,relative)
for i in range(S):
    if i == x:
        continue
    print "0\t" * (S+2), 
    for j in range(S):
        #if j == x:
        #    continue
        if j == i:
            print "1\t",
        else:
            print "0\t",
    
    for j in range(E):
        print "%d\t" % P[j,i],
    print 

print
print "#",
for i in range(E+2*S-1):
    print str(i+1) + "\t",
print
print "N\t",
if 'p' in flags[x]:
    print "pH",
else:
    print "H",

for i in range(S):
    print "\tT" + labels[i],
for i in range(S+E):
    print "\t" + full_labels[i],
print
for i, h in enumerate(pH):
    print i, "\t",
    print "%.4f" % h + '\t',
    for j in range(S):
        print "%.4f\t" % fullT[i,j],
    for j in range(S+E):
        print "%.4e\t" % finalC[i,j],
    print
#
#print finalC

for i in range(E+S):
    if toplot[i]:
        pyplot.plot(pH, finalC[:,i], formats[i],  label=full_labels[i])

pyplot.xlabel("pH")
pyplot.legend()
#pyplot.show()


