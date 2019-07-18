/*----------------------------------------------------------------------
  File    : relim.h
  Contents: relim algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2013.11.20 file created
            2014.08.22 interface of function relim() changed
            2014.08.28 functions rel_data() and rel_report() added
            2017.05.30 optional output compression with zlib added
----------------------------------------------------------------------*/
#ifndef __RELIM__
#define __RELIM__
#include "report.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- target pattern types --- */
#define REL_FREQ      ISR_FREQUENT  /* frequent item sets */
#define REL_FREQUENT  ISR_FREQUENT  /* frequent item sets */
#define REL_CLOSED    ISR_CLOSED    /* closed  frequent item sets */
#define REL_MAXIMAL   ISR_MAXIMAL   /* maximal frequent item sets */

/* --- triangular norms (t-norms) --- */
#define REL_MIN       0         /* minimum */
#define REL_NILP      1         /* nil-potent minimum */
#define REL_PROD      2         /* product */
#define REL_LUKA      3         /* Lukasiewicz */
#define REL_HAMA      4         /* Hamacher product */

/* --- evaluation measures --- */
#define REL_NONE      0         /* no measure */
#define REL_LDRATIO   1         /* binary log. of support quotient */

/* --- algorithm variants --- */
#define REL_BASIC     0         /* basic algorithm (list based) */
#define REL_TREE      1         /* tree-based algorithm */
#define REL_AUTO      0         /* automatic algorithm choice */

/* --- operation modes --- */
#define REL_FIM16     0x001f    /* use 16 items machine (bit rep.) */
#define REL_PERFECT   0x0020    /* prune with perfect extensions */
#define REL_PREFMT    0x1000    /* pre-format integer numbers */
#ifdef USE_ZLIB                 /* if optional output compression */
#define REL_ZLIB      0x4000    /* flag for output compression */
#endif
#define REL_DEFAULT   (REL_PERFECT|REL_FIM16|REL_PREFMT)
#ifdef NDEBUG
#define REL_NOCLEAN   0x8000    /* do not clean up memory */
#else                           /* in function relim() */
#define REL_NOCLEAN   0         /* in debug version */
#endif                          /* always clean up memory */
#define REL_VERBOSE   INT_MIN   /* verbose message output */

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct _relim           /* recursive elimination miner */
RELIM;                          /* (opaque structure) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern RELIM* relim_create (int target, double smin, double sins,
                            ITEM zmin, ITEM zmax,
                            int tnorm, double twgt,
                            int eval, double thresh,
                            int algo, int mode);
extern void   relim_delete (RELIM *relim, int deldar);
extern int    relim_data   (RELIM *relim, TABAG *tabag, int sort);
extern int    relim_report (RELIM *relim, ISREPORT *report);
extern int    relim_mine   (RELIM *relim, ITEM sort);
#endif
