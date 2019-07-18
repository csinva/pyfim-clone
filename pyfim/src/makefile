#-----------------------------------------------------------------------
# File    : makefile
# Contents: build shared object for frequent item set mining in Python
# Author  : Christian Borgelt
# History : 2011.07.08 file created
#           2011.07.22 module ruleval added
#           2013.10.19 module patspec added
#           2013.11.07 modules tabread and scanner removed
#           2013.08.08 modules gamma, chi2 and ruleval added
#           2014.08.25 modules ista and carpenter added
#           2015.03.04 interrupt signal handling added
#           2015.03.05 preprocessor definitions XXX_ABORT added 
#           2015.08.19 module patred added
#           2016.04.20 creation of dependency files added
#-----------------------------------------------------------------------
SHELL    = /bin/bash
THISDIR  = ../../pyfim/src
UTILDIR  = ../../util/src
MATHDIR  = ../../math/src
TRACTDIR = ../../tract/src
APRIDIR  = ../../apriori/src
ECLATDIR = ../../eclat/src
FPGDIR   = ../../fpgrowth/src
SAMDIR   = ../../sam/src
RELIMDIR = ../../relim/src
CARPDIR  = ../../carpenter/src
ISTADIR  = ../../ista/src
ACCDIR   = ../../accretion/src

CC       = gcc -std=c99
PY2INC   = -I/usr/include/python2.7
PY3INC   = -I/usr/include/python3.6
INCS     = -I$(UTILDIR) -I$(MATHDIR)  -I$(TRACTDIR) \
           -I$(APRIDIR) -I$(ECLATDIR) -I$(FPGDIR)   \
           -I$(SAMDIR)  -I$(RELIMDIR) -I$(CARPDIR)  \
           -I$(ISTADIR) -I$(ACCDIR)
CFBASE   = -Wall -Wextra -Wconversion -Wno-unused-parameter \
           -Wimplicit-fallthrough=0 -Wno-overlength-strings \
           -pedantic -c -fPIC $(ADDFLAGS)
CFLAGS   = $(CFBASE) -DNDEBUG -DQUIET -O3 -funroll-loops
# CFLAGS   = $(CFBASE) -DNDEBUG -O3 -DALIGN8
# CFLAGS   = $(CFBASE) -g
# CFLAGS   = $(CFBASE) -g -DSTORAGE
# CFLAGS   = $(CFBASE) -g -DALIGN8
XCFLAGS  = $(CFLAGS) -fvisibility=hidden

LD       = gcc
LDFLAGS  = $(ADDFLAGS)
LIBS     = -lm $(ADDLIBS)

# ADDOBJS  = $(UTILDIR)/storage.o

HDRS     = $(UTILDIR)/fntypes.h   $(UTILDIR)/arrays.h   \
           $(UTILDIR)/memsys.h    $(UTILDIR)/symtab.h   \
           $(UTILDIR)/random.h    $(UTILDIR)/sigint.h   \
           $(MATHDIR)/ruleval.h   $(TRACTDIR)/tract.h   \
           $(TRACTDIR)/patspec.h  $(TRACTDIR)/clomax.h  \
           $(TRACTDIR)/report.h   $(TRACTDIR)/patred.h  \
           $(APRIDIR)/istree.h    $(APRIDIR)/apriori.h  \
           $(ECLATDIR)/eclat.h    $(FPGDIR)/fpgrowth.h  \
           $(FPGDIR)/fpgpsp.h     $(SAMDIR)/sam.h       \
           $(RELIMDIR)/relim.h    $(CARPDIR)/repotree.h \
           $(CARPDIR)/carpenter.h $(ISTADIR)/ista.h     \
           $(ACCDIR)/accretion.h
OBJS     = arrays.o memsys.o idmap.o random.o sigint.o \
           chi2.o gamma.o ruleval.o \
           tatree.o fim16.o patspec.o clomax.o report.o patred.o \
           istree.o apriori.o eclat.o fpgrowth.o sam.o relim.o \
           repotree.o carpenter.o pfxtree.o pattree.o ista.o \
           accretion.o fpgpsp.o 
OBJS2    = $(OBJS) py2fim.o $(ADDOBJS)
OBJS3    = $(OBJS) py3fim.o $(ADDOBJS)

#-----------------------------------------------------------------------
# Build Shared Library
#-----------------------------------------------------------------------
all:          fim.so

fim.so:       $(OBJS2) makefile
	$(LD) $(LDFLAGS) -shared $(OBJS2) -o $@

py3:          fim3.so

fim3.so:      $(OBJS3) makefile
	$(LD) $(LDFLAGS) -shared $(OBJS3) -o $@
	mv fim3.so fim.so

#-----------------------------------------------------------------------
# Array Operations
#-----------------------------------------------------------------------
arrays.o:     $(UTILDIR)/fntypes.h
arrays.o:     $(UTILDIR)/arrays.h   $(UTILDIR)/arrays.c makefile
	$(CC) $(XCFLAGS) $(INCS) $(UTILDIR)/arrays.c -o $@

#-----------------------------------------------------------------------
# Memory Management System for Objects of Equal Size
#-----------------------------------------------------------------------
memsys.o:     $(UTILDIR)/memsys.h   $(UTILDIR)/memsys.c makefile
	$(CC) $(XCFLAGS) $(INCS) $(UTILDIR)/memsys.c -o $@

#-----------------------------------------------------------------------
# Symbol Table Management
#-----------------------------------------------------------------------
idmap.o:      $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h
idmap.o:      $(UTILDIR)/symtab.h   $(UTILDIR)/symtab.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DIDMAPFN $(UTILDIR)/symtab.c -o $@

#-----------------------------------------------------------------------
# Random Number Generator Management
#-----------------------------------------------------------------------
random.o:     $(UTILDIR)/random.h   $(UTILDIR)/random.c makefile
	$(CC) $(XCFLAGS) $(INCS) $(UTILDIR)/random.c -o $@

#-----------------------------------------------------------------------
# Interrupt Signal Handling
#-----------------------------------------------------------------------
sigint.o:     $(UTILDIR)/sigint.h   $(UTILDIR)/sigint.c makefile
	$(CC) $(XCFLAGS) $(INCS) $(UTILDIR)/sigint.c -o $@

#-----------------------------------------------------------------------
# Mathematical Functions
#-----------------------------------------------------------------------
gamma.o:      $(MATHDIR)/gamma.h    $(MATHDIR)/gamma.c makefile
	$(CC) $(XCFLAGS) $(INCS)    $(MATHDIR)/gamma.c -o $@

chi2.o:       $(MATHDIR)/gamma.h
chi2.o:       $(MATHDIR)/chi2.h     $(MATHDIR)/chi2.c makefile
	$(CC) $(XCFLAGS) $(INCS) $(MATHDIR)/chi2.c -o $@

ruleval.o:    $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h  \
              $(UTILDIR)/symtab.h   $(MATHDIR)/gamma.h   \
              $(MATHDIR)/chi2.h     $(TRACTDIR)/tract.h
ruleval.o:    $(MATHDIR)/ruleval.h  $(MATHDIR)/ruleval.c makefile
	$(CC) $(XCFLAGS) $(INCS) $(MATHDIR)/ruleval.c -o $@

#-----------------------------------------------------------------------
# 16 Items Machine
#-----------------------------------------------------------------------
fim16.o:      $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h  \
              $(UTILDIR)/symtab.h   $(UTILDIR)/memsys.h  \
              $(MATHDIR)/gamma.h    $(MATHDIR)/chi2.h    \
              $(TRACTDIR)/tract.h   $(TRACTDIR)/fim16.h  \
              $(TRACTDIR)/report.h  $(TRACTDIR)/clomax.h
fim16.o:      $(TRACTDIR)/fim16.h   $(TRACTDIR)/fim16.c makefile
	$(CC) $(XCFLAGS) $(INCS) $(TRACTDIR)/fim16.c -o $@

fim16.d:      $(TRACTDIR)/fim16.c
	$(CC) -MM $(XCFLAGS) $(INCS) $(TRACTDIR)/fim16.c > fim16.d

#-----------------------------------------------------------------------
# Item and Transaction Management
#-----------------------------------------------------------------------
tatree.o:     $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h \
              $(UTILDIR)/symtab.h   $(UTILDIR)/random.h
tatree.o:     $(TRACTDIR)/tract.h   $(TRACTDIR)/tract.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DTATREEFN -DTA_SURR \
              $(TRACTDIR)/tract.c -o $@

tatree.d:     $(TRACTDIR)/tract.c
	$(CC) -MM $(XCFLAGS) $(INCS) -DTATREEFN -DTA_SURR \
              $(TRACTDIR)/tract.c > tatree.d

#-----------------------------------------------------------------------
# Item Set Reporter Management
#-----------------------------------------------------------------------
patspec.o:    $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h \
              $(UTILDIR)/symtab.h   $(UTILDIR)/random.h \
              $(MATHDIR)/gamma.h    $(TRACTDIR)/tract.h
patspec.o:    $(TRACTDIR)/patspec.h $(TRACTDIR)/patspec.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DPSP_ESTIM \
              $(TRACTDIR)/patspec.c -o $@

clomax.o:     $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h \
              $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h \
              $(TRACTDIR)/tract.h
clomax.o:     $(TRACTDIR)/clomax.h  $(TRACTDIR)/clomax.c makefile
	$(CC) $(XCFLAGS) $(INCS) $(TRACTDIR)/clomax.c -o $@

report.o:     $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h \
              $(UTILDIR)/symtab.h   $(TRACTDIR)/tract.h \
              $(TRACTDIR)/patspec.h
report.o:     $(TRACTDIR)/report.h  $(TRACTDIR)/report.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DISR_PATSPEC -DISR_CLOMAX \
              -DISR_NONAMES $(TRACTDIR)/report.c -o $@

#-----------------------------------------------------------------------
# Pattern Set Reduction Functions
#-----------------------------------------------------------------------
patred.o:     $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h \
              $(UTILDIR)/symtab.h   $(TRACTDIR)/tract.h \
              $(TRACTDIR)/report.h 
patred.o:     $(TRACTDIR)/patred.h  $(TRACTDIR)/patred.c makefile
	$(CC) $(XCFLAGS) $(INCS) $(TRACTDIR)/patred.c -o $@

#-----------------------------------------------------------------------
# Apriori
#-----------------------------------------------------------------------
istree.o:     $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h  \
              $(UTILDIR)/symtab.h   $(MATHDIR)/gamma.h   \
              $(MATHDIR)/chi2.h     $(MATHDIR)/ruleval.h \
              $(TRACTDIR)/tract.h   $(TRACTDIR)/report.h 
istree.o:     $(APRIDIR)/istree.h   $(APRIDIR)/istree.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DTATREEFN $(APRIDIR)/istree.c -o $@

apriori.o:    $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h   \
              $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h   \
              $(UTILDIR)/sigint.h   $(MATHDIR)/ruleval.h  \
              $(TRACTDIR)/tract.h   $(TRACTDIR)/report.h  \
              $(TRACTDIR)/patspec.h $(APRIDIR)/istree.h
apriori.o:    $(APRIDIR)/apriori.h  $(APRIDIR)/apriori.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DISR_CLOMAX -DAPR_ABORT \
              $(APRIDIR)/apriori.c -o $@

apriori.d:    $(APRIDIR)/apriori.c
	$(CC) -MM $(XCFLAGS) $(INCS) -DISR_CLOMAX -DAPR_ABORT \
              $(APRIDIR)/apriori.c > apriori.d

#-----------------------------------------------------------------------
# Eclat
#-----------------------------------------------------------------------
eclat.o:      $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h  \
              $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h  \
              $(UTILDIR)/sigint.h   $(MATHDIR)/ruleval.h \
              $(TRACTDIR)/tract.h   $(TRACTDIR)/fim16.h  \
              $(TRACTDIR)/report.h  $(TRACTDIR)/clomax.h \
              $(TRACTDIR)/patspec.h $(APRIDIR)/istree.h
eclat.o:      $(ECLATDIR)/eclat.h   $(ECLATDIR)/eclat.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DECL_ABORT $(ECLATDIR)/eclat.c -o $@

eclat.d:       $(ECLATDIR)/eclat.c
	$(CC) -MM $(XCFLAGS) $(INCS) -DECL_ABORT \
              $(ECLATDIR)/eclat.c > eclat.d

#-----------------------------------------------------------------------
# FP-growth
#-----------------------------------------------------------------------
fpgrowth.o:   $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h  \
              $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h  \
              $(UTILDIR)/sigint.h   $(MATHDIR)/ruleval.h \
              $(TRACTDIR)/tract.h   $(TRACTDIR)/fim16.h  \
              $(TRACTDIR)/report.h  $(TRACTDIR)/clomax.h \
              $(TRACTDIR)/patspec.h $(APRIDIR)/istree.h
fpgrowth.o:   $(FPGDIR)/fpgrowth.h  $(FPGDIR)/fpgrowth.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DFPG_ABORT $(FPGDIR)/fpgrowth.c -o $@

fpgrowth.d:   $(FPGDIR)/fpgrowth.c
	$(CC) -MM $(XCFLAGS) $(INCS) -DFPG_ABORT \
              $(FPGDIR)/fpgrowth.c > fpgrowth.d

fpgpsp.o:     $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h  \
              $(UTILDIR)/random.h   $(UTILDIR)/memsys.h  \
              $(UTILDIR)/symtab.h   $(UTILDIR)/sigint.h  \
              $(MATHDIR)/ruleval.h  $(TRACTDIR)/tract.h  \
              $(TRACTDIR)/report.h  $(TRACTDIR)/clomax.h \
              $(TRACTDIR)/patspec.h $(APRIDIR)/istree.h  \
              $(FPGDIR)/fpgrowth.h
fpgpsp.o:     $(FPGDIR)/fpgpsp.h $(FPGDIR)/fpgpsp.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DFPG_ABORT $(FPGDIR)/fpgpsp.c -o $@

fpgpsp.d:     $(FPGDIR)/fpgpsp.c
	$(CC) -MM $(XCFLAGS) $(INCS) -DFPG_ABORT \
              $(FPGDIR)/fpgpsp.c > fpgpsp.d

#-----------------------------------------------------------------------
# SaM
#-----------------------------------------------------------------------
sam.o:        $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h  \
              $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h  \
              $(UTILDIR)/sigint.h   $(TRACTDIR)/tract.h  \
              $(TRACTDIR)/fim16.h   $(TRACTDIR)/report.h \
              $(TRACTDIR)/clomax.h  $(TRACTDIR)/patspec.h
sam.o:        $(SAMDIR)/sam.h $(SAMDIR)/sam.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DSAM_ABORT $(SAMDIR)/sam.c -o $@

sam.d:        $(SAMDIR)/sam.c
	$(CC) -MM $(XCFLAGS) $(INCS) -DSAM_ABORT $(SAMDIR)/sam.c > sam.d

#-----------------------------------------------------------------------
# RElim
#-----------------------------------------------------------------------
relim.o:      $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h  \
              $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h  \
              $(UTILDIR)/sigint.h   $(TRACTDIR)/tract.h  \
              $(TRACTDIR)/fim16.h   $(TRACTDIR)/report.h \
              $(TRACTDIR)/clomax.h  $(TRACTDIR)/patspec.h
relim.o:      $(RELIMDIR)/relim.h $(RELIMDIR)/relim.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DRELIM_ABORT $(RELIMDIR)/relim.c -o $@

relim.d:      $(RELIMDIR)/relim.c
	$(CC) -MM $(XCFLAGS) $(INCS) -DRELIM_ABORT \
              $(RELIMDIR)/relim.c > relim.d

#-----------------------------------------------------------------------
# Carpenter
#-----------------------------------------------------------------------
repotree.o:   $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h  \
              $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h  \
              $(TRACTDIR)/tract.h   $(TRACTDIR)/report.h \
              $(TRACTDIR)/clomax.h
repotree.o:   $(CARPDIR)/repotree.h $(CARPDIR)/repotree.c makefile
	$(CC) $(XCFLAGS) $(INCS) $(CARPDIR)/repotree.c -o $@

repotree.d:   $(CARPDIR)/repotree.c
	$(CC) -MM $(XCFLAGS) $(INCS) $(CARPDIR)/repotree.c > repotree.d

carpenter.o:  $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h  \
              $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h  \
              $(UTILDIR)/sigint.h   $(TRACTDIR)/tract.h  \
              $(TRACTDIR)/report.h  $(TRACTDIR)/clomax.h \
              $(TRACTDIR)/patspec.h $(CARPDIR)/repotree.h
carpenter.o:  $(CARPDIR)/carpenter.h $(CARPDIR)/carpenter.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DCARP_ABORT \
              $(CARPDIR)/carpenter.c -o $@

carpenter.d:  $(CARPDIR)/carpenter.c
	$(CC) -MM $(XCFLAGS) $(INCS) -DCARP_ABORT \
              $(CARPDIR)/carpenter.c > carpenter.d

#-----------------------------------------------------------------------
# IsTa
#-----------------------------------------------------------------------
pfxtree.o:    $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h  \
              $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h  \
              $(TRACTDIR)/tract.h   $(TRACTDIR)/report.h \
              $(TRACTDIR)/clomax.h
pfxtree.o:    $(ISTADIR)/pfxtree.h $(ISTADIR)/pfxtree.c makefile
	$(CC) $(XCFLAGS) $(INCS) $(ISTADIR)/pfxtree.c -o $@

pfxtree.d:    $(ISTADIR)/pfxtree.c
	$(CC) -MM $(XCFLAGS) $(INCS) $(ISTADIR)/pfxtree.c > pfxtree.d

pattree.o:    $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h  \
              $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h  \
              $(TRACTDIR)/tract.h   $(TRACTDIR)/report.h \
              $(TRACTDIR)/clomax.h
pattree.o:    $(ISTADIR)/pattree.h $(ISTADIR)/pattree.c makefile
	$(CC) $(XCFLAGS) $(INCS) $(ISTADIR)/pattree.c -o $@

pattree.d:    $(ISTADIR)/pattree.c
	$(CC) -MM $(XCFLAGS) $(INCS) $(ISTADIR)/pattree.c > pattree.d

ista.o:       $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h   \
              $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h   \
              $(UTILDIR)/sigint.h   $(TRACTDIR)/tract.h   \
              $(TRACTDIR)/report.h  $(TRACTDIR)/clomax.h  \
              $(TRACTDIR)/patspec.h $(ISTADIR)/pfxtree.h  \
              $(ISTADIR)/pattree.h
ista.o:       $(ISTADIR)/ista.h $(ISTADIR)/ista.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DISTA_ABORT $(ISTADIR)/ista.c -o $@

ista.d:       $(ISTADIR)/ista.c
	$(CC) -MM $(XCFLAGS) $(INCS) -DISTA_ABORT \
              $(ISTADIR)/ista.c > ista.d

#-----------------------------------------------------------------------
# Accretion
#-----------------------------------------------------------------------
accretion.o:  $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h   \
              $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h   \
              $(UTILDIR)/sigint.h   $(MATHDIR)/ruleval.h  \
              $(TRACTDIR)/tract.h   $(TRACTDIR)/report.h  \
              $(TRACTDIR)/clomax.h  $(TRACTDIR)/patspec.h 
accretion.o:  $(ACCDIR)/accretion.h $(ACCDIR)/accretion.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DISR_CLOMAX -DACC_ABORT \
              $(ACCDIR)/accretion.c -o $@

accretion.d:  $(ACCDIR)/accretion.c
	$(CC) -MM $(XCFLAGS) $(INCS) -DISR_CLOMAX -DACC_ABORT \
              $(ACCDIR)/accretion.c > accretion.d

#-----------------------------------------------------------------------
# Python Stuff
#-----------------------------------------------------------------------
py2fim.o:     $(HDRS)
py2fim.o:     pyfim.c makefile
	$(CC) $(CFLAGS) $(INCS) $(PY2INC) -fno-strict-aliasing \
              pyfim.c -o $@

py2fim.d:     pyfim.c makefile
	$(CC) -MM $(CFLAGS) $(INCS) $(PY2INC) -fno-strict-aliasing \
              pyfim.c > $@

py3fim.o:     $(HDRS)
py3fim.o:     pyfim.c makefile
	$(CC) $(CFLAGS) $(INCS) $(PY3INC) -fno-strict-aliasing \
              pyfim.c -o $@

py3fim.d:     pyfim.c makefile
	$(CC) -MM $(CFLAGS) $(INCS) $(PY3INC) -fno-strict-aliasing \
              pyfim.c > $@

#-----------------------------------------------------------------------
# Source Distribution Packages
#-----------------------------------------------------------------------
dist:
	$(MAKE) clean
	cd ../..; rm -f pyfim.zip pyfim.tar.gz; \
	cp setup_fim.py setup.py; \
        zip -rq pyfim.zip pyfim/{src,ex,doc} \
          apriori/{src,ex,doc} eclat/{src,ex,doc} \
          fpgrowth/{src,ex,doc} sam/{src,ex,doc} relim/{src,ex,doc} \
          carpenter/{src,ex,doc} ista/{src,ex,doc} \
          accretion/{src,ex,doc} \
          tract/src/{tract.[ch],fim16.[ch],patspec.[ch]} \
          tract/src/{clomax.[ch],report.[ch],patred.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          math/src/{gamma.[ch],chi2.[ch],ruleval.[ch]} \
          math/src/{makefile,math.mak} math/doc \
          util/src/{fntypes.h,arrays.[ch],symtab.[ch]} \
          util/src/{memsys.[ch],random.[ch],sigint.[ch]} \
          util/src/{makefile,util.mak} util/doc \
          setup_fim.py setup.py vs*.bat; \
        tar cfz pyfim.tar.gz pyfim/{src,ex,doc} \
          apriori/{src,ex,doc} eclat/{src,ex,doc} \
          fpgrowth/{src,ex,doc} sam/{src,ex,doc} relim/{src,ex,doc} \
          carpenter/{src,ex,doc} ista/{src,ex,doc} \
          accretion/{src,ex,doc} \
          tract/src/{tract.[ch],fim16.[ch],patspec.[ch]} \
          tract/src/{clomax.[ch],report.[ch],patred.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          math/src/{gamma.[ch],chi2.[ch],ruleval.[ch]} \
          math/src/{makefile,math.mak} math/doc \
          util/src/{fntypes.h,arrays.[ch],symtab.[ch]} \
          util/src/{memsys.[ch],random.[ch],sigint.[ch]} \
          util/src/{makefile,util.mak} util/doc \
          setup_fim.py setup.py vs*.bat; \
	rm -f setup.py

#-----------------------------------------------------------------------
# Install
#-----------------------------------------------------------------------
install:
	mv fim.so $(HOME)/lib

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
clean:
	rm -f *.d *.o *.so *~ *.flc core
