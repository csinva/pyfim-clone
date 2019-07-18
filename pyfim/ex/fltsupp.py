#!/usr/bin/python
#-----------------------------------------------------------------------
# File    : fltsupp.py
# Contents: test fim with floating point support
#           (needs compilation with "make ADDFLAGS=-DSUPP=double")
# Author  : Christian Borgelt
# History : 2016.10.06 file created
#-----------------------------------------------------------------------
from sys import argv
from fim import apriori, eclat, fpgrowth, fim

#-----------------------------------------------------------------------

tid = int(argv[1])
if   tid < -2:
    print(fpgrowth.__doc__)
elif tid < -1:
    print(eclat.__doc__)
elif tid <  0:
    print(apriori.__doc__)
else:
    tracts = {(1, 2, 3): 0.5,
              (1, 4, 5): 1.2,
              (2, 3, 4): 0.8,
              (1, 2, 3, 4): 0.3,
              (2, 3): 1.5,
              (1, 2, 4): 0.9,
              (4, 5): 0.6,
              (1, 2, 3, 4): 1.0,
              (3, 4, 5): 0.7}
    print('transactions:')
    for t in tracts: print(t, tracts[t])
    if   tid < 1:
        print  ('apriori(tracts, supp=-1.6, zmin=2):')
        for r in apriori(tracts, supp=-1.6, zmin=2): print(r)
    elif tid < 2:
        print  ('eclat(tracts, supp=-1.6, zmin=2):')
        for r in eclat(tracts, supp=-1.6, zmin=2): print(r)
    elif tid < 3:
        print  ('fpgrowth(tracts, supp=-1.6, zmin=2):')
        for r in fpgrowth(tracts, supp=-1.6, zmin=2): print(r)
    else:
        print  ('fim(tracts, supp=-1.6, zmin=2, report=\'#\'):')
        for r in fim(tracts, supp=-1.6, zmin=2, report='#'): print(r)
