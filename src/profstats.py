#!/usr/bin/python3

import pstats
p = pstats.Stats('profile.dat')
p.strip_dirs().sort_stats('time').print_stats()
