/*----------------------------------------------------------------------
  File    : fpgrowth.h
  Contents: fpgrowth algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2011.08.22 file created
            2011.09.21 available variants and modes reorganized
            2014.08.07 association rule generation/evaluation added
            2014.08.19 adapted to modified item set reporter interface
            2014.08.21 parameter 'body' added to function fpgrowth()
            2014.08.28 functions fpg_data() and fpg_report() added
            2016.11.20 fpgrowth miner object and interface introduced
            2017.05.30 optional output compression with zlib added
----------------------------------------------------------------------*/
#ifndef __FPGROWTH__
#define __FPGROWTH__
#include "report.h"
#include "ruleval.h"
#include "istree.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- target pattern types --- */
#define FPG_FREQ      ISR_FREQUENT  /* frequent item sets */
#define FPG_FREQUENT  ISR_FREQUENT  /* frequent item sets */
#define FPG_CLOSED    ISR_CLOSED    /* closed  frequent item sets */
#define FPG_MAXIMAL   ISR_MAXIMAL   /* maximal frequent item sets */
#define FPG_GENERAS   ISR_GENERAS   /* generators */
#define FPG_RULES     ISR_RULES     /* association rules */

/* --- data preparation modes --- */
#define FPG_NORECODE  0x0001    /* do not sort and recode items */
#define FPG_NOFILTER  0x0002    /* do not filter transactions by size */
#define FPG_NOSORT    0x0004    /* do not sort items and transactions */
#define FPG_NOREDUCE  0x0008    /* do not reduce transactions */
#define FPG_NOPACK    0x0010    /* do not pack most frequent items */
#define FPG_SURR      (FPG_NORECODE|FPG_NOFILTER|FPG_NOREDUCE)

/* --- evaluation measures --- */
/* most definitions in ruleval.h */
#define FPG_LDRATIO   RE_FNCNT  /* binary log. of support quotient */
#define FPG_INVBXS    IST_INVBXS/* inval. eval. below exp. supp. */

/* --- aggregation modes --- */
#define FPG_NONE      IST_NONE  /* no aggregation (use first value) */
#define FPG_FIRST     IST_FIRST /* no aggregation (use first value) */
#define FPG_MIN       IST_MIN   /* minimum of measure values */
#define FPG_MAX       IST_MAX   /* maximum of measure values */
#define FPG_AVG       IST_AVG   /* average of measure values */

/* --- algorithm variants --- */
#define FPG_SIMPLE    0         /* simple  nodes (parent/link) */
#define FPG_COMPLEX   1         /* complex nodes (children/sibling) */
#define FPG_SINGLE    2         /* top-down processing on single tree */
#define FPG_TOPDOWN   3         /* top-down processing of the tree */
#define FPG_AUTO      4         /* automatic choice */

/* --- operation modes --- */
#define FPG_FIM16     0x001f    /* use 16 items machine (bit rep.) */
#define FPG_PERFECT   0x0020    /* perfect extension pruning */
#define FPG_REORDER   0x0040    /* reorder items in cond. databases */
#define FPG_ORIGSUPP  0x0080    /* use original support definition */
#define FPG_TAIL      0x0100    /* head union tail pruning */
#define FPG_PREFMT    0x1000    /* pre-format integer numbers */
#ifdef USE_ZLIB                 /* if optional output compression */
#define FPG_ZLIB      0x4000    /* flag for output compression */
#endif
#define FPG_DEFAULT   (FPG_PERFECT|FPG_REORDER|FPG_TAIL|FPG_FIM16)
#ifdef NDEBUG
#define FPG_NOCLEAN   0x8000    /* do not clean up memory */
#else                           /* in function fpgrowth() */
#define FPG_NOCLEAN   0         /* in debug version */
#endif                          /* always clean up memory */
#define FPG_VERBOSE   INT_MIN   /* verbose message output */

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct _fpgrowth        /* fpgrowth miner */
FPGROWTH;                       /* (opaque structure) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern FPGROWTH* fpg_create (int target, double smin, double smax,
                             double conf, ITEM zmin, ITEM zmax,
                             int eval, int agg, double thresh,
                             int algo, int mode);
extern void      fpg_delete (FPGROWTH *fpg, int deldar);
extern int       fpg_data   (FPGROWTH *fpg, TABAG *tabag,
                             int mode, int sort);
extern int       fpg_report (FPGROWTH *fpg, ISREPORT *report);
extern int       fpg_mine   (FPGROWTH *fpg, ITEM prune, int order);
#endif
