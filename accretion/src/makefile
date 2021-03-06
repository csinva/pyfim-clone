#-----------------------------------------------------------------------
# File    : makefile
# Contents: build accretion program (on Unix systems)
# Author  : Christian Borgelt
# History : 2011.06.22 file created from eclat makefile
#           2011.07.22 module ruleval added (rule evaluation)
#           2013.03.20 extended the requested warnings in CFBASE
#           2013.10.19 modules tabread and patspec added
#           2016.04.20 creation of dependency files added
#-----------------------------------------------------------------------
# For large file support (> 2GB) compile with
#   make ADDFLAGS=-D_FILE_OFFSET_BITS=64
#-----------------------------------------------------------------------
SHELL    = /bin/bash
THISDIR  = ../../accretion/src
UTILDIR  = ../../util/src
MATHDIR  = ../../math/src
TRACTDIR = ../../tract/src

CC       = gcc -std=c99
CFBASE   = -Wall -Wextra -Wno-unused-parameter -Wconversion \
           -pedantic -c $(ADDFLAGS)
CFLAGS   = $(CFBASE) -DNDEBUG -O3 -funroll-loops
# CFLAGS   = $(CFBASE) -g
# CFLAGS   = $(CFBASE) -g -DSTORAGE
INCS     = -I$(UTILDIR) -I$(MATHDIR) -I$(TRACTDIR)

LD       = gcc
LDFLAGS  = $(ADDFLAGS)
LIBS     = -lm $(ADDLIBS)

# ADDOBJS  = $(UTILDIR)/storage.o

HDRS     = $(UTILDIR)/fntypes.h  $(UTILDIR)/arrays.h   \
           $(UTILDIR)/symtab.h   $(UTILDIR)/escape.h   \
           $(UTILDIR)/tabread.h  $(UTILDIR)/tabwrite.h \
           $(UTILDIR)/scanner.h  $(UTILDIR)/error.h    \
           $(MATHDIR)/ruleval.h  $(TRACTDIR)/tract.h   \
           $(TRACTDIR)/patspec.h $(TRACTDIR)/report.h  \
           accretion.h
OBJS     = $(UTILDIR)/arrays.o   $(UTILDIR)/idmap.o    \
           $(UTILDIR)/escape.o   $(UTILDIR)/tabread.o  \
           $(UTILDIR)/tabwrite.o $(UTILDIR)/scform.o   \
           $(MATHDIR)/gamma.o    $(MATHDIR)/chi2.o     \
           $(MATHDIR)/ruleval.o  $(TRACTDIR)/taread.o  \
           $(TRACTDIR)/patspec.o $(TRACTDIR)/report.o  \
           accretion.o $(ADDOBJS)
PRGS     = accretion

#-----------------------------------------------------------------------
# Build Program
#-----------------------------------------------------------------------
all:          $(PRGS)

accretion:    $(OBJS) makefile
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

#-----------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------
accretion.o:  $(HDRS)
accretion.o:  accretion.c makefile
	$(CC) $(CFLAGS) $(INCS) -DACC_MAIN accretion.c -o $@

accretion.d:  accretion.c
	$(CC) -MM $(CFLAGS) $(INCS) -DACC_MAIN accretion.c > accretion.d

#-----------------------------------------------------------------------
# External Modules
#-----------------------------------------------------------------------
$(UTILDIR)/arrays.o:
	cd $(UTILDIR);  $(MAKE) arrays.o  ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/idmap.o:
	cd $(UTILDIR);  $(MAKE) idmap.o   ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/escape.o:
	cd $(UTILDIR);  $(MAKE) escape.o  ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/tabread.o:
	cd $(UTILDIR);  $(MAKE) tabread.o ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/scform.o:
	cd $(UTILDIR);  $(MAKE) scform.o  ADDFLAGS="$(ADDFLAGS)"
$(MATHDIR)/gamma.o:
	cd $(MATHDIR);  $(MAKE) gamma.o   ADDFLAGS="$(ADDFLAGS)"
$(MATHDIR)/chi2.o:
	cd $(MATHDIR);  $(MAKE) chi2.o    ADDFLAGS="$(ADDFLAGS)"
$(MATHDIR)/ruleval.o:
	cd $(MATHDIR);  $(MAKE) ruleval.o ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/taread.o:
	cd $(TRACTDIR); $(MAKE) taread.o  ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/patspec.o:
	cd $(TRACTDIR); $(MAKE) patspec.o ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/report.o:
	cd $(TRACTDIR); $(MAKE) report.o  ADDFLAGS="$(ADDFLAGS)"

#-----------------------------------------------------------------------
# Source Distribution Packages
#-----------------------------------------------------------------------
dist:
	$(MAKE) clean
	cd ../..; rm -f accretion.zip accretion.tar.gz; \
        zip -rq accretion.zip accretion/{src,ex,doc} \
          tract/src/{tract.[ch],patspec.[ch],report.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          math/src/{gamma.[ch],chi2.[ch],ruleval.[ch]} \
          math/src/{makefile,math.mak} math/doc \
          util/src/{fntypes.h,error.h} \
          util/src/{arrays.[ch],symtab.[ch],escape.[ch]} \
          util/src/{tabread.[ch],tabwrite.[ch],scanner.[ch]} \
          util/src/{makefile,util.mak} util/doc; \
        tar cfz accretion.tar.gz accretion/{src,ex,doc} \
          tract/src/{tract.[ch],patspec.[ch],report.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          math/src/{gamma.[ch],chi2.[ch],ruleval.[ch]} \
          math/src/{makefile,math.mak} math/doc \
          util/src/{fntypes.h,error.h} \
          util/src/{arrays.[ch],symtab.[ch],escape.[ch]} \
          util/src/{tabread.[ch],tabwrite.[ch],scanner.[ch]} \
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
	cd $(MATHDIR);  $(MAKE) clean
	cd $(UTILDIR);  $(MAKE) clean
