/*----------------------------------------------------------------------
  File    : carpenter.h
  Contents: carpenter algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2013.10.31 file created from carpenter.c
            2014.08.23 interface of function carpenter() changed
            2014.08.28 functions carp_data() and carp_report() added
            2017.05.30 optional output compression with zlib added
----------------------------------------------------------------------*/
#ifndef __CARPENTER__
#define __CARPENTER__
#include "report.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- target pattern types --- */
#define CARP_CLOSED   ISR_CLOSED    /* closed  frequent item sets */
#define CARP_MAXIMAL  ISR_MAXIMAL   /* maximal frequent item sets */

/* --- evaluation measures --- */
#define CARP_NONE     0         /* no measure */
#define CARP_LDRATIO  1         /* binary log. of support quotient */

/* --- carpenter variants --- */
#define CARP_AUTO     0         /* auto. choice based on table size */
#define CARP_TABLE    1         /* item occurrence counter table */
#define CARP_TIDLIST  2         /* transaction identifier lists */

/* --- operation modes --- */
#define CARP_PERFECT  0x0010    /* prune with perfect extensions */
#define CARP_FILTER   0x0020    /* filter maximal sets with repo. */
#define CARP_MAXONLY  0x0040    /* add only maximal sets to repo. */
#define CARP_COLLATE  0x0080    /* flag for collating transactions */
#define CARP_PREFMT   0x1000    /* pre-format integer numbers */
#ifdef USE_ZLIB                 /* if optional output compression */
#define CARP_ZLIB     0x4000    /* flag for output compression */
#endif
#define CARP_DEFAULT  (CARP_COLLATE|CARP_PERFECT)
#ifdef NDEBUG
#define CARP_NOCLEAN  0x8000    /* do not clean up memory */
#else                           /* in function carpenter() */
#define CARP_NOCLEAN  0         /* in debug version */
#endif                          /* always clean up memory */
#define CARP_VERBOSE  INT_MIN   /* verbose message output */

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct _carp            /* carpenter miner */
CARP;                           /* (opaque structure) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern CARP* carp_create (int target, double smin, double smax,
                          ITEM zmin, ITEM zmax, int eval, double thresh,
                          int algo, int mode);
extern void  carp_delete (CARP *carp, int deldar);
extern int   carp_data   (CARP *carp, TABAG *tabag, int sort);
extern int   carp_report (CARP *carp, ISREPORT *report);
extern int   carp_mine   (CARP *carp);
#endif
