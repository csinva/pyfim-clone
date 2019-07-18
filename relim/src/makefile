#-----------------------------------------------------------------------
# File    : makefile
# Contents: build relim program (on Unix systems)
# Author  : Christian Borgelt
# History : 2004.11.05 file created from eclat makefile
#           2004.12.12 program relx added
#           2008.08.11 adapted to name change from vecops to arrays
#           2009.10.16 external modules memsys and pfxtree added
#           2010.02.12 module pfxtree replaced by module clomax
#           2010.08.22 module escape added (for module tabread)
#           2010.10.08 changed standard from -ansi to -std=c99
#           2011.08.31 external module fim16 added (16 items machine)
#           2013.03.20 extended the requested warnings in CFBASE
#           2016.04.20 creation of dependency files added
#-----------------------------------------------------------------------
# For large file support (> 2GB) compile with
#   make ADDFLAGS=-D_FILE_OFFSET_BITS=64
#-----------------------------------------------------------------------
SHELL    = /bin/bash
THISDIR  = ../../relim/src
UTILDIR  = ../../util/src
TRACTDIR = ../../tract/src

CC       = gcc -std=c99
CFBASE   = -Wall -Wextra -Wno-unused-parameter -Wconversion \
           -pedantic -c $(ADDFLAGS)
CFLAGS   = $(CFBASE) -DNDEBUG -O3 -funroll-loops
# CFLAGS   = $(CFBASE) -g
# CFLAGS   = $(CFBASE) -g -DSTORAGE
INCS     = -I$(UTILDIR) -I$(TRACTDIR)

LD       = gcc
LDFLAGS  = $(ADDFLAGS)
LIBS     = -lm $(ADDLIBS)

# ADDOBJS  = $(UTILDIR)/storage.o

HDRS     = $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h   \
           $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h   \
           $(UTILDIR)/tabread.h  $(UTILDIR)/tabwrite.h \
           $(UTILDIR)/error.h    $(TRACTDIR)/tract.h   \
           $(TRACTDIR)/patspec.h $(TRACTDIR)/clomax.h  \
           $(TRACTDIR)/report.h  $(TRACTDIR)/fim16.h
OBJS     = $(UTILDIR)/arrays.o   $(UTILDIR)/memsys.o   \
           $(UTILDIR)/idmap.o    $(UTILDIR)/escape.o   \
           $(UTILDIR)/tabread.o  $(UTILDIR)/tabwrite.o \
           $(UTILDIR)/scform.o   $(TRACTDIR)/taread.o  \
           $(TRACTDIR)/patspec.o $(TRACTDIR)/clomax.o  \
           $(TRACTDIR)/repcm.o   $(TRACTDIR)/fim16.o   \
           relim.o $(ADDOBJS)
PRGS     = relim

#-----------------------------------------------------------------------
# Build Program
#-----------------------------------------------------------------------
all:          $(PRGS)

relim:        $(OBJS) makefile
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

#-----------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------
relim.o:      $(HDRS)
relim.o:      relim.h relim.c makefile
	$(CC) $(CFLAGS) $(INCS) -DRELIM_MAIN relim.c -o $@

relim.d:      relim.c
	$(CC) -MM $(CFLAGS) $(INCS) -DRELIM_MAIN relim.c > relim.d

#-----------------------------------------------------------------------
# External Modules
#-----------------------------------------------------------------------
$(UTILDIR)/arrays.o:
	cd $(UTILDIR);  $(MAKE) arrays.o   ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/memsys.o:
	cd $(UTILDIR);  $(MAKE) memsys.o   ADDFLAGS="$(ADDFLAGS)"
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
$(TRACTDIR)/taread.o:
	cd $(TRACTDIR); $(MAKE) taread.o   ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/patspec.o:
	cd $(TRACTDIR); $(MAKE) patspec.o  ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/clomax.o:
	cd $(TRACTDIR); $(MAKE) clomax.o   ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/repcm.o:
	cd $(TRACTDIR); $(MAKE) repcm.o    ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/fim16.o:
	cd $(TRACTDIR); $(MAKE) fim16.o    ADDFLAGS="$(ADDFLAGS)"

#-----------------------------------------------------------------------
# Source Distribution Packages
#-----------------------------------------------------------------------
dist:
	$(MAKE) clean
	cd ../..; rm -f relim.zip relim.tar.gz; \
        zip -rq relim.zip relim/{src,ex,doc} \
          tract/src/{tract.[ch],fim16.[ch]} \
          tract/src/{patspec.[ch],clomax.[ch],report.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          util/src/{fntypes.h,error.h} \
          util/src/{arrays.[ch],memsys.[ch],symtab.[ch]} \
          util/src/{escape.[ch],tabread.[ch],tabwrite.[ch]} \
          util/src/scanner.[ch] \
          util/src/{makefile,util.mak} util/doc; \
        tar cfz relim.tar.gz relim/{src,ex,doc} \
          tract/src/{tract.[ch],fim16.[ch]} \
          tract/src/{patspec.[ch],clomax.[ch],report.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          util/src/{fntypes.h,error.h} \
          util/src/{arrays.[ch],memsys.[ch],symtab.[ch]} \
          util/src/{escape.[ch],tabread.[ch],tabwrite.[ch]} \
          util/src/scanner.[ch] \
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
	cd $(UTILDIR);  $(MAKE) clean
	cd $(TRACTDIR); $(MAKE) localclean
