/*----------------------------------------------------------------------
  File    : ista.h
  Contents: finding frequent item sets by intersecting transactions
  Author  : Christian Borgelt
  History : 2014.08.24 file created
            2014.08.28 functions ista_data() and ista_report() added
            2017.03.24 ista miner object and interface introduced
            2017.05.30 optional output compression with zlib added
----------------------------------------------------------------------*/
#ifndef __ISTA__
#define __ISTA__
#include "report.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- target pattern types --- */
#define ISTA_CLOSED   ISR_CLOSED    /* closed  frequent item sets */
#define ISTA_MAXIMAL  ISR_MAXIMAL   /* maximal frequent item sets */

/* --- evaluation measures --- */
#define ISTA_NONE     0         /* no measure */
#define ISTA_LDRATIO  1         /* binary log. of support quotient */

/* --- algorithm variants --- */
#define ISTA_PREFIX   0         /* use a prefix   tree */
#define ISTA_PATRICIA 1         /* use a patricia tree */
#define ISTA_AUTO     2         /* automatic choice */

/* --- operation modes --- */
#define ISTA_PRUNE    0x0010    /* prune the prefix/patricia tree */
#define ISTA_FILTER   0x0020    /* filter maximal sets with repo. */
#define ISTA_MAXONLY  0x0040    /* add only maximal sets to repo. */
#define ISTA_PREFMT   0x1000    /* pre-format integer numbers */
#ifdef USE_ZLIB                 /* if optional output compression */
#define ISTA_ZLIB     0x4000    /* flag for output compression */
#endif
#define ISTA_DEFAULT  ISTA_PRUNE
#ifdef NDEBUG
#define ISTA_NOCLEAN  0x8000    /* do not clean up memory */
#else                           /* in function ista() */
#define ISTA_NOCLEAN  0         /* in debug version */
#endif                          /* always clean up memory */
#define ISTA_VERBOSE  INT_MIN   /* verbose message output */

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct _ista            /* ista miner */
ISTA;                           /* (opaque structure) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern ISTA* ista_create (int target, double smin, double smax,
                          ITEM zmin, ITEM zmax, int eval, double thresh,
                          int algo, int mode);
extern void  ista_delete (ISTA *ista, int deldar);
extern int   ista_data   (ISTA *ista, TABAG *tabag, int sort);
extern int   ista_report (ISTA *ista, ISREPORT *report);
extern int   ista_mine   (ISTA *ista);
#endif
