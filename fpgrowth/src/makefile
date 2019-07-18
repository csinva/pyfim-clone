#-----------------------------------------------------------------------
# File    : makefile
# Contents: build fpgrowth program (on Unix systems)
# Author  : Christian Borgelt
# History : 2004.11.21 file created from eclat makefile
#           2004.12.10 general memory management system added
#           2008.08.11 adapted to name change from vecops to arrays
#           2009.10.19 external module pfxtree added
#           2010.02.08 external module pfxtree replaced by module clomax
#           2010.08.22 module escape added (for module tabread)
#           2010.10.08 changed standard from -ansi to -std=c99
#           2013.03.20 extended the requested warnings in CFBASE
#           2014.08.21 extended by module istree from apriori source
#           2016.04.20 creation of dependency files added
#-----------------------------------------------------------------------
# For large file support (> 2GB) compile with
#   make ADDFLAGS=-D_FILE_OFFSET_BITS=64
#-----------------------------------------------------------------------
SHELL    = /bin/bash
THISDIR  = ../../fpgrowth/src
UTILDIR  = ../../util/src
MATHDIR  = ../../math/src
TRACTDIR = ../../tract/src
APRIDIR  = ../../apriori/src

CC       = gcc -std=c99
CFBASE   = -Wall -Wextra -Wno-unused-parameter -Wconversion \
           -pedantic -c $(ADDFLAGS)
CFLAGS   = $(CFBASE) -DNDEBUG -O3 -funroll-loops
# CFLAGS   = $(CFBASE) -DNDEBUG -O3 -funroll-loops -DBENCH
# CFLAGS   = $(CFBASE) -g
# CFLAGS   = $(CFBASE) -g -DSTORAGE
INCS     = -I$(UTILDIR) -I$(MATHDIR) -I$(TRACTDIR) -I$(APRIDIR)

LD       = gcc
LDFLAGS  = $(ADDFLAGS)
LIBS     = -lm $(ADDLIBS)

# ADDOBJS  = $(UTILDIR)/storage.o

HDRS     = $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h   \
           $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h   \
           $(UTILDIR)/error.h    $(MATHDIR)/ruleval.h  \
           $(TRACTDIR)/tract.h   $(TRACTDIR)/patspec.h \
           $(TRACTDIR)/clomax.h  $(TRACTDIR)/report.h  \
           $(APRIDIR)/istree.h   fpgrowth.h
OBJS     = $(UTILDIR)/memsys.o   $(UTILDIR)/arrays.o   \
           $(UTILDIR)/idmap.o    $(UTILDIR)/escape.o   \
           $(UTILDIR)/tabread.o  $(UTILDIR)/tabwrite.o \
           $(UTILDIR)/scform.o   $(MATHDIR)/gamma.o    \
           $(MATHDIR)/chi2.o     $(MATHDIR)/ruleval.o  \
           $(TRACTDIR)/clomax.o  $(TRACTDIR)/repcm.o   \
           $(TRACTDIR)/fim16.o   $(APRIDIR)/istree.o $(ADDOBJS)

FPGOBJS  = $(OBJS)               $(TRACTDIR)/taread.o  \
           $(TRACTDIR)/patspec.o fpgmain.o
PSPOBJS  = $(OBJS)               $(UTILDIR)/random.o   \
           $(TRACTDIR)/tars.o    $(TRACTDIR)/pspest.o  \
           fpgrowth.o pspmain.o

PRGS     = fpgrowth fpgpsp

#-----------------------------------------------------------------------
# Build Program
#-----------------------------------------------------------------------
all:          $(PRGS)

fpgrowth:     $(FPGOBJS) makefile
	$(LD) $(LDFLAGS) $(FPGOBJS) $(LIBS) -o $@

fpgpsp:       $(PSPOBJS) makefile
	$(LD) $(LDFLAGS) $(PSPOBJS) $(LIBS) -lpthread -o $@

#-----------------------------------------------------------------------
# Main Programs
#-----------------------------------------------------------------------
fpgmain.o:    $(HDRS)               $(UTILDIR)/tabread.h \
              $(UTILDIR)/tabwrite.h $(TRACTDIR)/fim16.h
fpgmain.o:    fpgrowth.c makefile
	$(CC) $(CFLAGS) $(INCS) -DFPG_MAIN fpgrowth.c -o $@

fpgmain.d:    fpgrowth.c
	$(CC) -MM $(CFLAGS) $(INCS) -DFPG_MAIN fpgrowth.c > fpgmain.d

pspmain.o:    $(HDRS)               $(UTILDIR)/tabread.h \
              $(UTILDIR)/tabwrite.h $(UTILDIR)/random.h fpgpsp.h
pspmain.o:    fpgpsp.c makefile
	$(CC) $(CFLAGS) $(INCS) -DFPGPSP_MAIN fpgpsp.c -o $@

pspmain.d:    fpgpsp.c
	$(CC) -MM $(CFLAGS) $(INCS) -DFPGPSP_MAIN fpgpsp.c > pspmain.d

#-----------------------------------------------------------------------
# FP-growth as a module
#-----------------------------------------------------------------------
fpgrowth.o:   $(HDRS)               $(TRACTDIR)/fim16.h
fpgrowth.o:   fpgrowth.h fpgrowth.c makefile
	$(CC) $(CFLAGS) $(INCS) fpgrowth.c -o $@

fpgrowth.d:   fpgrowth.c
	$(CC) -MM $(CFLAGS) $(INCS) fpgrowth.c > fpgrowth.d

#-----------------------------------------------------------------------
# External Modules
#-----------------------------------------------------------------------
$(UTILDIR)/memsys.o:
	cd $(UTILDIR);  $(MAKE) memsys.o   ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/arrays.o:
	cd $(UTILDIR);  $(MAKE) arrays.o   ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/idmap.o:
	cd $(UTILDIR);  $(MAKE) idmap.o    ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/escape.o:
	cd $(UTILDIR);  $(MAKE) escape.o   ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/tabread.o:
	cd $(UTILDIR);  $(MAKE) tabread.o  ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/tabwrite.o:
	cd $(UTILDIR);  $(MAKE) tabwrite.o ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/scform.o:
	cd $(UTILDIR);  $(MAKE) scform.o   ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/random.o:
	cd $(UTILDIR);  $(MAKE) random.o   ADDFLAGS="$(ADDFLAGS)"
$(MATHDIR)/gamma.o:
	cd $(MATHDIR);  $(MAKE) gamma.o    ADDFLAGS="$(ADDFLAGS)"
$(MATHDIR)/chi2.o:
	cd $(MATHDIR);  $(MAKE) chi2.o     ADDFLAGS="$(ADDFLAGS)"
$(MATHDIR)/ruleval.o:
	cd $(MATHDIR);  $(MAKE) ruleval.o  ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/taread.o:
	cd $(TRACTDIR); $(MAKE) taread.o   ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/tars.o:
	cd $(TRACTDIR); $(MAKE) tars.o     ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/patspec.o:
	cd $(TRACTDIR); $(MAKE) patspec.o  ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/pspest.o:
	cd $(TRACTDIR); $(MAKE) pspest.o   ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/clomax.o:
	cd $(TRACTDIR); $(MAKE) clomax.o   ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/repcm.o:
	cd $(TRACTDIR); $(MAKE) repcm.o    ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/fim16.o:
	cd $(TRACTDIR); $(MAKE) fim16.o    ADDFLAGS="$(ADDFLAGS)"
$(APRIDIR)/istree.o:
	cd $(APRIDIR);  $(MAKE) istree.o   ADDFLAGS="$(ADDFLAGS)"

#-----------------------------------------------------------------------
# Source Distribution Packages
#-----------------------------------------------------------------------
dist:
	$(MAKE) clean
	cd ../..; rm -f fpgrowth.zip fpgrowth.tar.gz; \
        zip -rq fpgrowth.zip fpgrowth/{src,ex,doc} \
          apriori/src/{istree.[ch],makefile,apriori.mak} \
          tract/src/{tract.[ch],fim16.[ch]} \
          tract/src/{patspec.[ch],clomax.[ch],report.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          math/src/{gamma.[ch],chi2.[ch],ruleval.[ch]} \
          math/src/{makefile,math.mak} math/doc \
          util/src/{fntypes.h,error.h} \
          util/src/{arrays.[ch],memsys.[ch],symtab.[ch]} \
          util/src/{escape.[ch],tabread.[ch],tabwrite.[ch]} \
          util/src/{scanner.[ch],random.[ch]} \
          util/src/{makefile,util.mak} util/doc; \
        tar cfz fpgrowth.tar.gz fpgrowth/{src,ex,doc} \
          apriori/src/{istree.[ch],makefile,apriori.mak} \
          tract/src/{tract.[ch],fim16.[ch]} \
          tract/src/{patspec.[ch],clomax.[ch],report.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          math/src/{gamma.[ch],chi2.[ch],ruleval.[ch]} \
          math/src/{makefile,math.mak} math/doc \
          util/src/{fntypes.h,error.h} \
          util/src/{arrays.[ch],memsys.[ch],symtab.[ch]} \
          util/src/{escape.[ch],tabread.[ch],tabwrite.[ch]} \
          util/src/{scanner.[ch],random.[ch]} \
          util/src/{makefile,util.mak} util/doc

#-----------------------------------------------------------------------
# Installation
#-----------------------------------------------------------------------
install:
	cp $(PRGS) $(HOME)/bin

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
localclean:
	rm -f *.d *.o *~ *.flc core $(PRGS)

clean:
	$(MAKE) localclean
	cd $(APRIDIR);  $(MAKE) localclean
	cd $(TRACTDIR); $(MAKE) localclean
	cd $(MATHDIR);  $(MAKE) clean
	cd $(UTILDIR);  $(MAKE) clean
