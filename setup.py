#!/usr/bin/python

from distutils.core import setup, Extension
from os             import remove

dirs   =['util/src',
         'math/src',
         'tract/src',
         'apriori/src',
         'eclat/src',
         'fpgrowth/src',
         'sam/src',
         'relim/src',
         'carpenter/src',
         'ista/src',
         'accretion/src',
         'pyfim/src']
headers=['util/src/arrays.h',
         'util/src/memsys.h',
         'util/src/symtab.h',
         'util/src/random.h',
         'util/src/sigint.h',
         'util/src/fntypes.h',
         'math/src/gamma.h',
         'math/src/chi2.h',
         'math/src/ruleval.h',
         'tract/src/tract.h',
         'tract/src/fim16.h',
         'tract/src/patspec.h',
         'tract/src/clomax.h',
         'tract/src/report.h',
         'tract/src/patred.h',
         'apriori/src/istree.h',
         'apriori/src/apriori.h',
         'eclat/src/eclat.h',
         'fpgrowth/src/fpgrowth.h',
         'fpgrowth/src/fpgpsp.h',
         'sam/src/sam.h',
         'relim/src/relim.h',
         'carpenter/src/repotree.h',
         'carpenter/src/carpenter.h',
         'ista/src/pfxtree.h',
         'ista/src/pattree.h',
         'ista/src/ista.h',
         'accretion/src/accretion.h']
sources=['util/src/arrays.c',
         'util/src/memsys.c',
         'util/src/symtab.c',
         'util/src/random.c',
         'util/src/sigint.c',
         'math/src/chi2.c',
         'math/src/gamma.c',
         'math/src/ruleval.c',
         'tract/src/tract.c',
         'tract/src/fim16.c',
         'tract/src/patspec.c',
         'tract/src/clomax.c',
         'tract/src/report.c',
         'tract/src/patred.c',
         'apriori/src/istree.c',
         'apriori/src/apriori.c',
         'eclat/src/eclat.c',
         'fpgrowth/src/fpgrowth.c',
         'fpgrowth/src/fpgpsp.c',
         'sam/src/sam.c',
         'relim/src/relim.c',
         'carpenter/src/repotree.c',
         'carpenter/src/carpenter.c',
         'ista/src/pfxtree.c',
         'ista/src/pattree.c',
         'ista/src/ista.c',
         'accretion/src/accretion.c',
         'pyfim/src/pyfim.c']
macros =[('NDEBUG',      None),
         ('IDMAPFN',     None),
         ('TATREEFN',    None),
         ('TA_SURR',     None),
         ('PSP_ESTIM',   None),
         ('ISR_PATSPEC', None),
         ('ISR_CLOMAX',  None),
         ('ISR_NONAMES', None),
         ('ACC_ABORT',   None),
         ('APR_ABORT',   None),
         ('ECL_ABORT',   None),
         ('FPG_ABORT',   None),
         ('SAM_ABORT',   None),
         ('RELIM_ABORT', None),
         ('CARP_ABORT',  None),
         ('ISTA_ABORT',  None)]

with open('MANIFEST.in', 'wt') as out:
    for h in headers: out.write('include '+h+'\n')
with open('README', 'wt') as out:
    out.write('Call help() on the functions\n')
    out.write('  fim\n')
    out.write('  arules\n')
    out.write('  apriori\n')
    out.write('  eclat\n')
    out.write('  fpgrowth\n')
    out.write('  sam\n')
    out.write('  relim\n')
    out.write('  carpenter\n')
    out.write('  ista\n')
    out.write('  accretion\n')
    out.write('  apriacc\n')
    out.write('  patspec\n')
    out.write('  estpsp\n')
    out.write('for explanations about their parameters.\n')

setup(name='fim',
      version='6.28 (2017.03.24)',
      description='Frequent Item Set Mining and Association Rule Induction for Python',
      long_description='Frequent Item Set Mining and Association Rule Induction for Python',
      author='Christian Borgelt',
      author_email='christian@borgelt.net',
      maintainer='Christian Borgelt',
      maintainer_email='christian@borgelt.net',
      url='http://www.borgelt.net/pyfim.html',
      download_url='http://www.borgelt.net/pyfim.html',
      platforms='Linux/Unix,Microsoft Windows',
      license='MIT License (Expat License)',
      ext_modules=[Extension('fim', sources, include_dirs=dirs,
                             define_macros=macros)])

remove('MANIFEST.in')
remove('README')
