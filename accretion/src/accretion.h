/*----------------------------------------------------------------------
  File    : accretion.h
  Contents: accretion algorithm for identifying neural assemblies
  Author  : Christian Borgelt
  History : 2011.07.14 file created
            2011.07.22 adapted to new module ruleval (rule evaluation)
            2013.01.31 definition of flag ACC_INVBXS added
            2013.03.29 adapted to type changes in module tract
            2014.08.28 functions acc_data() and acc_report() added
            2016.11.15 accretion miner object and interface introduced
            2017.05.30 optional output compression with zlib added
----------------------------------------------------------------------*/
#ifndef __ACCRETION__
#define __ACCRETION__
#include "report.h"
#include "ruleval.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- evaluation measures --- */
/* most definitions in ruleval.h */
#define ACC_INVBXS    INT_MIN   /* invalidate stat. below exp. supp. */

/* --- operation modes --- */
#define ACC_PREFMT    0x1000    /* pre-format integer numbers */
#ifdef USE_ZLIB                 /* if optional output compression */
#define ACC_ZLIB      0x4000    /* flag for output compression */
#endif
#define ACC_DEFAULT   0x0000    /* default operation mode */
#ifdef NDEBUG
#define ACC_NOCLEAN   0x8000    /* do not clean up memory */
#else                           /* in function accretion() */
#define ACC_NOCLEAN   0         /* in debug version */
#endif                          /* always clean up memory */
#define ACC_VERBOSE   INT_MIN   /* verbose message output */

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct _accret          /* accretion miner */
ACCRET;                         /* (opaque structure) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern ACCRET* accret_create (int target, double smin, double smax,
                              ITEM zmin, ITEM zmax,
                              int stat, double siglvl, int mode);
extern void    accret_delete (ACCRET *accret, int deldar);
extern int     accret_data   (ACCRET *accret, TABAG *tabag, int sort);
extern int     accret_report (ACCRET *accret, ISREPORT *report);
extern int     accret_mine   (ACCRET *accret, ITEM maxext);
#endif
