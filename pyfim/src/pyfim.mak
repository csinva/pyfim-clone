#-----------------------------------------------------------------------
# File    : pyfim.mak
# Contents: build pyfim dynamic link library (on Windows systems)
# Author  : Christian Borgelt
# History : 2012.08.03 file created
#           2013.10.19 module patspec added
#           2013.11.08 properly adapted to Windows/Microsoft Visual C
#           2014.08.25 module ista added
#           2015.03.04 interrupt signal handling added
#           2015.03.05 preprocessor definitions XXX_ABORT added 
#           2015.08.19 module patred added
#           2016.04.20 completed dependencies on header files
#-----------------------------------------------------------------------
!IFNDEF PY2DIR
PY2DIR   = "C:\Program Files\Python27"
!ENDIF
!IFNDEF PY3DIR
PY3DIR   = "C:\Program Files\Python35"
!ENDIF
THISDIR  = ..\..\pyfim\src
UTILDIR  = ..\..\util\src
MATHDIR  = ..\..\math\src
TRACTDIR = ..\..\tract\src
APRIDIR  = ..\..\apriori\src
ECLATDIR = ..\..\eclat\src
FPGDIR   = ..\..\fpgrowth\src
SAMDIR   = ..\..\sam\src
RELIMDIR = ..\..\relim\src
CARPDIR  = ..\..\carpenter\src
ISTADIR  = ..\..\ista\src
ACCDIR   = ..\..\accretion\src

CC       = cl.exe
DEFS     = /D WIN32 /D NDEBUG /D _CONSOLE /D _CRT_SECURE_NO_WARNINGS \
           /D QUIET
CFLAGS   = /nologo /MD /W3 /O2 /GS- $(DEFS) /c
PY2INC   = /I $(PY2DIR)\include
PY3INC   = /I $(PY3DIR)\include
INCS     = /I $(UTILDIR) /I $(MATHDIR)  /I $(TRACTDIR) \
           /I $(APRIDIR) /I $(ECLATDIR) /I $(FPGDIR)   \
           /I $(SAMDIR)  /I $(RELIMDIR) /I $(CARPDIR)  \
           /I $(ISTADIR) /I $(ACCDIR)

LD       = link.exe
LDFLAGS  = /DLL /nologo /incremental:no
LIBS2    = /LIBPATH:$(PY2DIR)\libs
LIBS3    = /LIBPATH:$(PY3DIR)\libs

HDRS     = $(UTILDIR)\fntypes.h   $(UTILDIR)\arrays.h   \
           $(UTILDIR)\memsys.h    $(UTILDIR)\symtab.h   \
           $(UTILDIR)\random.h    $(UTILDIR)\sigint.h   \
           $(MATHDIR)\ruleval.h   $(TRACTDIR)\tract.h   \
           $(TRACTDIR)\patspec.h  $(TRACTDIR)\clomax.h  \
           $(TRACTDIR)\report.h   $(TRACTDIR)\patred.h  \
           $(APRIDIR)\istree.h    $(APRIDIR)\apriori.h  \
           $(ECLATDIR)\eclat.h    $(FPGDIR)\fpgrowth.h  \
           $(FPGDIR)\fpgpsp.h     $(SAMDIR)\sam.h       \
           $(RELIMDIR)\relim.h    $(CARPDIR)\repotree.h \
           $(CARPDIR)\carpenter.h $(ISTADIR)\ista.h     \
           $(ACCDIR)\accretion.h
OBJS     = arrays.obj memsys.obj idmap.obj random.obj sigint.obj \
           chi2.obj gamma.obj ruleval.obj tatree.obj \
           fim16.obj patspec.obj clomax.obj report.obj patred.obj \
           istree.obj apriori.obj eclat.obj fpgrowth.obj \
           sam.obj relim.obj repotree.obj carpenter.obj \
           pfxtree.obj pattree.obj ista.obj accretion.obj \
           fpgpsp.obj
OBJS2    = $(OBJS) py2fim.obj $(ADDOBJS)
OBJS3    = $(OBJS) py3fim.obj $(ADDOBJS)

#-----------------------------------------------------------------------
# Build Dynamic Link Library
#-----------------------------------------------------------------------
all:          fim.pyd

fim.pyd:      $(OBJS2) pyfim.mak
	$(LD) $(LDFLAGS) $(LIBS2) $(OBJS2) /out:$@ /IMPLIB:fim.lib

py3:          fim3.pyd

fim3.pyd:     $(OBJS3) pyfim.mak
	$(LD) $(LDFLAGS) $(LIBS3) $(OBJS3) /out:$@ /IMPLIB:fim.lib
	-@del fim.pyd
	-@ren fim3.pyd fim.pyd

#-----------------------------------------------------------------------
# Array Operations
#-----------------------------------------------------------------------
arrays.obj:   $(UTILDIR)\fntypes.h
arrays.obj:   $(UTILDIR)\arrays.h   $(UTILDIR)\arrays.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(UTILDIR)\arrays.c /Fo$@

#-----------------------------------------------------------------------
# Memory Management System for Objects of Equal Size
#-----------------------------------------------------------------------
memsys.obj:   $(UTILDIR)\memsys.h   $(UTILDIR)\memsys.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(UTILDIR)\memsys.c /Fo$@

#-----------------------------------------------------------------------
# Symbol Table Management
#-----------------------------------------------------------------------
idmap.obj:    $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h
idmap.obj:    $(UTILDIR)\symtab.h $(UTILDIR)\symtab.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D IDMAPFN $(UTILDIR)\symtab.c /Fo$@

#-----------------------------------------------------------------------
# Random Number Generator Management
#-----------------------------------------------------------------------
random.obj:   $(UTILDIR)\random.h   $(UTILDIR)\random.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(UTILDIR)\random.c /Fo$@

#-----------------------------------------------------------------------
# Interrupt Signal Handling
#-----------------------------------------------------------------------
sigint.obj:   $(UTILDIR)\sigint.h   $(UTILDIR)\sigint.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(UTILDIR)\sigint.c /Fo$@

#-----------------------------------------------------------------------
# Mathematical Functions
#-----------------------------------------------------------------------
gamma.obj:    $(MATHDIR)\gamma.h    $(MATHDIR)\gamma.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(MATHDIR)\gamma.c /Fo$@

chi2.obj:     $(MATHDIR)\gamma.h
chi2.obj:     $(MATHDIR)\chi2.h     $(MATHDIR)\chi2.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(MATHDIR)\chi2.c /Fo$@

ruleval.obj:  $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h  \
              $(UTILDIR)\symtab.h   $(MATHDIR)\gamma.h   \
              $(MATHDIR)\chi2.h     $(TRACTDIR)\tract.h
ruleval.obj:  $(MATHDIR)\ruleval.h  $(MATHDIR)\ruleval.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(MATHDIR)\ruleval.c /Fo$@

#-----------------------------------------------------------------------
# 16 Items Machine
#-----------------------------------------------------------------------
fim16.obj:    $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h  \
              $(UTILDIR)\symtab.h   $(UTILDIR)\memsys.h  \
              $(MATHDIR)\gamma.h    $(MATHDIR)\chi2.h    \
              $(TRACTDIR)\tract.h   $(TRACTDIR)\fim16.h  \
              $(TRACTDIR)\report.h  $(TRACTDIR)\clomax.h
fim16.obj:    $(TRACTDIR)\fim16.h   $(TRACTDIR)\fim16.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(TRACTDIR)\fim16.c /Fo$@

#-----------------------------------------------------------------------
# Item and Transaction Management
#-----------------------------------------------------------------------
tatree.obj:   $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h \
              $(UTILDIR)\random.h   $(UTILDIR)\symtab.h
tatree.obj:   $(TRACTDIR)\tract.h $(TRACTDIR)\tract.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D TATREEFN /D TA_SURR \
              $(TRACTDIR)\tract.c /Fo$@

#-----------------------------------------------------------------------
# Item Set Reporter Management
#-----------------------------------------------------------------------
patspec.obj:  $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h \
              $(UTILDIR)\random.h   $(UTILDIR)\symtab.h \
              $(MATHDIR)\gamma.h    $(TRACTDIR)\tract.h
patspec.obj:  $(TRACTDIR)\patspec.h $(TRACTDIR)\patspec.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D PSP_ESTIM \
              $(TRACTDIR)\patspec.c /Fo$@

clomax.obj:   $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h \
              $(UTILDIR)\memsys.h   $(UTILDIR)\symtab.h \
              $(TRACTDIR)\tract.h
clomax.obj:   $(TRACTDIR)\clomax.h  $(TRACTDIR)\clomax.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(TRACTDIR)\clomax.c /Fo$@

report.obj:   $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h  \
              $(UTILDIR)\symtab.h   $(UTILDIR)\scanner.h \
              $(TRACTDIR)\tract.h   $(TRACTDIR)\patspec.h
report.obj:   $(TRACTDIR)\report.h  $(TRACTDIR)\report.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D ISR_PATSPEC /D ISR_CLOMAX \
              /D ISR_NONAMES $(TRACTDIR)\report.c /Fo$@

#-----------------------------------------------------------------------
# Pattern Set Reduction Functions
#-----------------------------------------------------------------------
patred.obj:   $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h \
              $(UTILDIR)\symtab.h   $(TRACTDIR)\tract.h \
              $(TRACTDIR)\report.h 
patred.obj:   $(TRACTDIR)\patred.h  $(TRACTDIR)\patred.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(TRACTDIR)\patred.c /Fo$@

#-----------------------------------------------------------------------
# Apriori
#-----------------------------------------------------------------------
istree.obj:   $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h  \
              $(UTILDIR)\symtab.h   $(MATHDIR)\chi2.h    \
              $(MATHDIR)\gamma.h    $(MATHDIR)\ruleval.h \
              $(TRACTDIR)\tract.h   $(TRACTDIR)\report.h 
istree.obj:   $(APRIDIR)\istree.h   $(APRIDIR)\istree.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D TATREEFN $(APRIDIR)\istree.c /Fo$@

apriori.obj:  $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h   \
              $(UTILDIR)\memsys.h   $(UTILDIR)\symtab.h   \
              $(UTILDIR)\sigint.h   $(MATHDIR)\ruleval.h  \
              $(TRACTDIR)\tract.h   $(TRACTDIR)\report.h  \
              $(TRACTDIR)\patspec.h $(APRIDIR)\istree.h
apriori.obj:  $(APRIDIR)\apriori.h  $(APRIDIR)\apriori.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D APR_ABORT $(APRIDIR)\apriori.c /Fo$@

#-----------------------------------------------------------------------
# Eclat
#-----------------------------------------------------------------------
eclat.obj:    $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h  \
              $(UTILDIR)\memsys.h   $(UTILDIR)\symtab.h  \
              $(UTILDIR)\sigint.h   $(MATHDIR)\ruleval.h \
              $(TRACTDIR)\tract.h   $(TRACTDIR)\fim16.h  \
              $(TRACTDIR)\report.h  $(TRACTDIR)\clomax.h \
              $(TRACTDIR)\patspec.h $(APRIDIR)\istree.h
eclat.obj:    $(ECLATDIR)\eclat.h   $(ECLATDIR)\eclat.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D ECL_ABORT $(ECLATDIR)\eclat.c /Fo$@

#-----------------------------------------------------------------------
# FP-growth
#-----------------------------------------------------------------------
fpgrowth.obj: $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h  \
              $(UTILDIR)\memsys.h   $(UTILDIR)\symtab.h  \
              $(UTILDIR)\sigint.h   $(MATHDIR)\ruleval.h \
              $(TRACTDIR)\tract.h   $(TRACTDIR)\fim16.h  \
              $(TRACTDIR)\report.h  $(TRACTDIR)\clomax.h \
              $(TRACTDIR)\patspec.h $(APRIDIR)\istree.h
fpgrowth.obj: $(FPGDIR)\fpgrowth.h  $(FPGDIR)\fpgrowth.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D FPG_ABORT $(FPGDIR)\fpgrowth.c /Fo$@

fpgpsp.obj:   $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h  \
              $(UTILDIR)\random.h   $(UTILDIR)\memsys.h  \
              $(UTILDIR)\symtab.h   $(UTILDIR)\sigint.h  \
              $(MATHDIR)\ruleval.h  $(TRACTDIR)\tract.h  \
              $(TRACTDIR)\report.h  $(TRACTDIR)\clomax.h \
              $(TRACTDIR)\patspec.h $(APRIDIR)\istree.h  \
              $(FPGDIR)\fpgrowth.h
fpgpsp.obj:   $(FPGDIR)\fpgpsp.h $(FPGDIR)\fpgpsp.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D FPG_ABORT $(FPGDIR)\fpgpsp.c /Fo$@

#-----------------------------------------------------------------------
# SaM
#-----------------------------------------------------------------------
sam.obj:      $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h  \
              $(UTILDIR)\memsys.h   $(UTILDIR)\symtab.h  \
              $(UTILDIR)\sigint.h   $(TRACTDIR)\tract.h  \
              $(TRACTDIR)\fim16.h   $(TRACTDIR)\report.h \
              $(TRACTDIR)\clomax.h  $(TRACTDIR)\patspec.h
sam.obj:      $(SAMDIR)\sam.h       $(SAMDIR)\sam.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D SAM_ABORT $(SAMDIR)\sam.c /Fo$@

#-----------------------------------------------------------------------
# RElim
#-----------------------------------------------------------------------
relim.obj:    $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h  \
              $(UTILDIR)\memsys.h   $(UTILDIR)\symtab.h  \
              $(UTILDIR)\sigint.h   $(TRACTDIR)\tract.h  \
              $(TRACTDIR)\fim16.h   $(TRACTDIR)\report.h \
              $(TRACTDIR)\clomax.h  $(TRACTDIR)\patspec.h
relim.obj:    $(RELIMDIR)\relim.h $(RELIMDIR)\relim.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D RELIM_ABORT $(RELIMDIR)\relim.c /Fo$@

#-----------------------------------------------------------------------
# Carpenter
#-----------------------------------------------------------------------
repotree.obj: $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h  \
              $(UTILDIR)\memsys.h   $(UTILDIR)\symtab.h  \
              $(TRACTDIR)\tract.h   $(TRACTDIR)\report.h \
              $(TRACTDIR)\clomax.h
repotree.obj: $(CARPDIR)\repotree.h $(CARPDIR)\repotree.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(CARPDIR)\repotree.c /Fo$@

carpenter.obj:  $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h  \
                $(UTILDIR)\memsys.h   $(UTILDIR)\symtab.h  \
                $(UTILDIR)\sigint.h   $(TRACTDIR)\tract.h  \
                $(TRACTDIR)\report.h  $(TRACTDIR)\clomax.h \
                $(TRACTDIR)\patspec.h $(CARPDIR)\repotree.h
carpenter.obj:  $(CARPDIR)\carpenter.h $(CARPDIR)\carpenter.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D CARP_ABORT \
              $(CARPDIR)\carpenter.c /Fo$@

#-----------------------------------------------------------------------
# IsTa
#-----------------------------------------------------------------------
pfxtree.obj:  $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h  \
              $(UTILDIR)\memsys.h   $(UTILDIR)\symtab.h  \
              $(TRACTDIR)\tract.h   $(TRACTDIR)\report.h \
              $(TRACTDIR)\clomax.h
pfxtree.obj:  $(ISTADIR)\pfxtree.h $(ISTADIR)\pfxtree.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(ISTADIR)\pfxtree.c /Fo$@

pattree.obj:  $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h  \
              $(UTILDIR)\memsys.h   $(UTILDIR)\symtab.h  \
              $(TRACTDIR)\tract.h   $(TRACTDIR)\report.h \
              $(TRACTDIR)\clomax.h
pattree.obj:  $(ISTADIR)\pattree.h $(ISTADIR)\pattree.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(ISTADIR)\pattree.c /Fo$@

ista.obj:     $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h   \
              $(UTILDIR)\memsys.h   $(UTILDIR)\symtab.h   \
              $(UTILDIR)\sigint.h   $(TRACTDIR)\tract.h   \
              $(TRACTDIR)\report.h  $(TRACTDIR)\clomax.h  \
              $(TRACTDIR)\patspec.h $(ISTADIR)\pfxtree.h  \
              $(ISTADIR)\pattree.h
ista.obj:     $(ISTADIR)\ista.h     $(ISTADIR)\ista.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D ISTA_ABORT $(ISTADIR)\ista.c /Fo$@

#-----------------------------------------------------------------------
# Accretion
#-----------------------------------------------------------------------
accretion.obj:  $(UTILDIR)\fntypes.h  $(UTILDIR)\arrays.h   \
                $(UTILDIR)\memsys.h   $(UTILDIR)\symtab.h   \
                $(UTILDIR)\sigint.h   $(MATHDIR)\ruleval.h  \
                $(TRACTDIR)\tract.h   $(TRACTDIR)\report.h  \
                $(TRACTDIR)\clomax.h  $(TRACTDIR)\patspec.h 
accretion.obj:  $(ACCDIR)\accretion.h $(ACCDIR)\accretion.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D ACC_ABORT \
              $(ACCDIR)\accretion.c /Fo$@

#-----------------------------------------------------------------------
# Python Stuff
#-----------------------------------------------------------------------
py2fim.obj:   $(HDRS)
py2fim.obj:   pyfim.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(PY2INC) pyfim.c /Fo$@

py3fim.obj:   $(HDRS)
py3fim.obj:   pyfim.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(PY3INC) pyfim.c /Fo$@

#-----------------------------------------------------------------------
# Install
#-----------------------------------------------------------------------
install:
	-@copy *.pyd ..\..\..\bin

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
clean:
	-@erase /Q *~ *.pyd *.obj *.idb *.pch *.exp *.lib $(PRGS)
