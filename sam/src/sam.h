/*----------------------------------------------------------------------
  File    : sam.h
  Contents: sam algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2013.11.20 file created
            2014.08.22 interface of function sam() changed
            2014.08.28 functions sam_data() and sam_report() added
            2016.11.23 eclat miner object and interface introduced
            2017.05.30 optional output compression with zlib added
----------------------------------------------------------------------*/
#ifndef __SAM__
#define __SAM__
#include "report.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- target pattern types --- */
#define SAM_FREQ      ISR_FREQUENT  /* frequent item sets */
#define SAM_FREQUENT  ISR_FREQUENT  /* frequent item sets */
#define SAM_CLOSED    ISR_CLOSED    /* closed  frequent item sets */
#define SAM_MAXIMAL   ISR_MAXIMAL   /* maximal frequent item sets */

/* --- triangular norms (t-norms) --- */
#define SAM_MIN       0         /* minimum */
#define SAM_NILP      1         /* nil-potent minimum */
#define SAM_PROD      2         /* product */
#define SAM_LUKA      3         /* Lukasiewicz */
#define SAM_HAMA      4         /* Hamacher product */

/* --- evaluation measures --- */
#define SAM_NONE      0         /* no measure */
#define SAM_LDRATIO   1         /* binary log. of support quotient */

/* --- sam variants --- */
#define SAM_BASIC     0         /* basic split and merge */
#define SAM_BSEARCH   1         /* allow binary search merge */
#define SAM_DOUBLE    2         /* use double source buffering */
#define SAM_TREE      3         /* transaction prefix trees */
#define SAM_AUTO      4         /* automatic algorithm choice */

/* --- operation modes --- */
#define SAM_FIM16     0x001f    /* use 16 items machine (bit rep.) */
#define SAM_PERFECT   0x0020    /* prune with perfect extensions */
#define SAM_PREFMT    0x1000    /* pre-format integer numbers */
#ifdef USE_ZLIB                 /* if optional output compression */
#define SAM_ZLIB      0x4000    /* flag for output compression */
#endif
#define SAM_DEFAULT   (SAM_PERFECT|SAM_FIM16)
#ifdef NDEBUG
#define SAM_NOCLEAN   0x8000    /* do not clean up memory */
#else                           /* in function sam() */
#define SAM_NOCLEAN   0         /* in debug version */
#endif                          /* always clean up memory */
#define SAM_VERBOSE   INT_MIN   /* verbose message output */

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct _sam             /* split and merge miner */
SAM;                            /* (opaque structure) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern SAM* sam_create (int target, double smin, double sins,
                        ITEM zmin, ITEM zmax, int tnorm, double twgt,
                        int eval, double thresh, int algo, int mode);
extern void sam_delete (SAM *sam, int deldar);
extern int  sam_data   (SAM *sam, TABAG *tabag, int sort);
extern int  sam_report (SAM *sam, ISREPORT *report);
extern int  sam_mine   (SAM *sam, TID merge);
#endif
