/*----------------------------------------------------------------------
  File    : apriori.h
  Contents: apriori algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2011.07.18 file created
            2011.10.18 several mode flags added
            2013.03.30 adapted to type changes in module tract
            2014.08.21 parameter 'body' added to function apriori()
            2014.08.28 functions apr_data() and apr_report() added
            2016.11.04 apriori miner object and interface introduced
            2017.05.30 optional output compression with zlib added
----------------------------------------------------------------------*/
#ifndef __APRIORI__
#define __APRIORI__
#include "report.h"
#include "ruleval.h"
#include "istree.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- target pattern types --- */
#define APR_FREQ      ISR_FREQUENT  /* frequent item sets */
#define APR_FREQUENT  ISR_FREQUENT  /* frequent item sets */
#define APR_CLOSED    ISR_CLOSED    /* closed  frequent item sets */
#define APR_MAXIMAL   ISR_MAXIMAL   /* maximal frequent item sets */
#define APR_GENERAS   ISR_GENERAS   /* generators */
#define APR_RULES     ISR_RULES     /* association rules */

/* --- data preparation modes --- */
#define APR_NORECODE  0x0001    /* do not sort and recode items */
#define APR_NOFILTER  0x0002    /* do not filter transactions by size */
#define APR_NOSORT    0x0004    /* do not sort items and transactions */
#define APR_NOREDUCE  0x0008    /* do not reduce transactions */

/* --- evaluation measures --- */
/* most measure definitions in ruleval.h */
#define APR_LDRATIO   RE_FNCNT  /* binary log. of support quotient */
#define APR_INVBXS    IST_INVBXS/* inval. eval. below exp. support */

/* --- aggregation modes --- */
#define APR_NONE      IST_NONE  /* no aggregation (use first value) */
#define APR_FIRST     IST_FIRST /* no aggregation (use first value) */
#define APR_MIN       IST_MIN   /* minimum of measure values */
#define APR_MAX       IST_MAX   /* maximum of measure values */
#define APR_AVG       IST_AVG   /* average of measure values */

/* --- algorithm variants --- */
#define APR_BASIC     0         /* basic algorithm (dummy) */
#define APR_AUTO      0         /* automatic algorithm choice (dummy) */

/* --- operation modes --- */
#define APR_ORIGSUPP  0x0080    /* use original support definition */
#define APR_PERFECT   IST_PERFECT  /* perfect extension pruning */
#define APR_TATREE    0x0200    /* use transaction tree */
#define APR_POST      0x0400    /* use a-posteriori pruning */
#define APR_PREFMT    0x1000    /* pre-format integer numbers */
#ifdef USE_ZLIB                 /* if optional output compression */
#define APR_ZLIB      0x4000    /* flag for output compression */
#endif
#define APR_DEFAULT   (APR_PERFECT|APR_TATREE)
#ifdef NDEBUG
#define APR_NOCLEAN   0x8000    /* do not clean up memory */
#else                           /* in function apriori() */
#define APR_NOCLEAN   0         /* in debug version */
#endif                          /* always clean up memory */
#define APR_VERBOSE   INT_MIN   /* verbose message output */

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct _apriori         /* apriori miner */
APRIORI;                        /* (opaque structure) */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern APRIORI* apriori_create (int target, double smin, double smax,
                                double conf, ITEM zmin, ITEM zmax,
                                int eval, int agg, double thresh,
                                int algo, int mode);
extern void     apriori_delete (APRIORI *apriori, int deldar);
extern int      apriori_data   (APRIORI *apriori, TABAG *tabag,
                                int mode, int sort);
extern int      apriori_report (APRIORI *apriori, ISREPORT *report);
extern int      apriori_mine   (APRIORI *apriori, ITEM prune,
                                double filter, int order);
#endif
