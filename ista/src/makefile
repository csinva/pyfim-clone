#-----------------------------------------------------------------------
# File    : makefile
# Contents: build ista program (intersect transactions) (Unix systems)
# Author  : Christian Borgelt
# History : 2010.10.26 file created
#           2010.07.12 module clomax added (for module report)
#           2010.08.22 module escape added (for module tabread)
#           2010.10.08 changed standard from -ansi to -std=c99
#           2012.07.16 patricia tree module added
#           2013.03.20 extended the requested warnings in CFBASE
#           2016.04.20 creation of dependency files added
#-----------------------------------------------------------------------
# For large file support (> 2GB) compile with
#   make ADDFLAGS=-D_FILE_OFFSET_BITS=64
#-----------------------------------------------------------------------
SHELL    = /bin/bash
THISDIR  = ../../ista/src
UTILDIR  = ../../util/src
TRACTDIR = ../../tract/src

CC       = gcc -std=c99
CFBASE   = -Wall -Wextra -Wno-unused-parameter -Wconversion \
           -pedantic -c $(ADDFLAGS)
CFLAGS   = $(CFBASE) -DNDEBUG -O3 -funroll-loops
# CFLAGS   = $(CFBASE) -DNDEBUG -O3 -funroll-loops -DBENCH
# CFLAGS   = $(CFBASE) -g
# CFLAGS   = $(CFBASE) -g -DSTORAGE
INCS     = -I$(UTILDIR) -I$(TRACTDIR)

LD       = gcc
LDFLAGS  = $(ADDFLAGS)
LIBS     = -lm $(ADDLIBS)

# ADDOBJS  = $(UTILDIR)/storage.o

HDRS_T   = $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h   \
           $(UTILDIR)/memsys.h   $(UTILDIR)/symtab.h   \
           $(TRACTDIR)/tract.h   $(TRACTDIR)/clomax.h  \
           $(TRACTDIR)/report.h
HDRS     = $(HDRS_T)             $(UTILDIR)/error.h    \
           $(UTILDIR)/tabread.h  $(UTILDIR)/tabwrite.h \
           $(TRACTDIR)/patspec.h pfxtree.h pattree.h ista.h
OBJS     = $(UTILDIR)/arrays.o   $(UTILDIR)/memsys.o   \
           $(UTILDIR)/idmap.o    $(UTILDIR)/escape.o   \
           $(UTILDIR)/tabread.o  $(UTILDIR)/tabwrite.o \
           $(UTILDIR)/scform.o   $(TRACTDIR)/taread.o  \
           $(TRACTDIR)/patspec.o $(TRACTDIR)/clomax.o  \
           $(TRACTDIR)/repcm.o   pfxtree.o pattree.o   \
           ista.o $(ADDOBJS)
PRGS     = ista

#-----------------------------------------------------------------------
# Build Program
#-----------------------------------------------------------------------
all:          $(PRGS)

ista:         $(OBJS) makefile
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

#-----------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------
ista.o:       $(HDRS)
ista.o:       ista.c makefile
	$(CC) $(CFLAGS) $(INCS) -DISTA_MAIN ista.c -o $@

ista.d:       ista.c
	$(CC) -MM $(CFLAGS) $(INCS) -DISTA_MAIN ista.c > ista.d

#-----------------------------------------------------------------------
# Prefix Tree Management
#-----------------------------------------------------------------------
pfxtree.o:    $(HDRS_T)
pfxtree.o:    pfxtree.h pfxtree.c makefile
	$(CC) $(CFLAGS) $(INCS) pfxtree.c -o $@

pfxtree.d:    pfxtree.c
	$(CC) -MM $(CFLAGS) $(INCS) pfxtree.c > pfxtree.d

#-----------------------------------------------------------------------
# Patricia Tree Management
#-----------------------------------------------------------------------
pattree.o:    $(HDRS_T)
pattree.o:    pattree.h pattree.c makefile
	$(CC) $(CFLAGS) $(INCS) pattree.c -o $@

pattree.d:    pattree.c
	$(CC) -MM $(CFLAGS) $(INCS) pattree.c > pattree.d

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

#-----------------------------------------------------------------------
# Source Distribution Packages
#-----------------------------------------------------------------------
dist:
	$(MAKE) clean
	cd ../..; rm -f ista.zip ista.tar.gz; \
        zip -rq ista.zip ista/{src,ex,doc} \
          tract/src/{tract.[ch],patspec.[ch],clomax.[ch],report.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          util/src/{fntypes.h,error.h} \
          util/src/{arrays.[ch],memsys.[ch],symtab.[ch]} \
          util/src/{escape.[ch],tabread.[ch],tabwrite.[ch]} \
          util/src/scanner.[ch] \
          util/src/{makefile,util.mak} util/doc; \
        tar cfz ista.tar.gz ista/{src,ex,doc} \
          tract/src/{tract.[ch],patspec.[ch],clomax.[ch],report.[ch]} \
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
	cd $(TRACTDIR); $(MAKE) localclean
	cd $(UTILDIR);  $(MAKE) clean
