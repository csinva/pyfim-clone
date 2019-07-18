/*----------------------------------------------------------------------
  File    : ista.c
  Contents: finding frequent item sets by intersecting transactions
  Author  : Christian Borgelt
  History : 2009.10.26 file created
            2009.10.30 first version for closed item sets completed
            2009.11.05 default item sorting direction reversed
            2009.11.06 processing order of transactions optimized
            2009.11.07 item order and transaction order made options
            2009.11.18 pruning for higher minimum support added
            2010.06.25 support-based pruning added to intersection
            2010.07.14 output file made optional (for benchmarking)
            2010.08.19 item selection file added as optional input
            2010.08.22 adapted to modified modules tabread and tract
            2010.10.15 adapted to modified interface of module report
            2010.11.24 adapted to modified error reporting (tract)
            2010.12.11 adapted to a generic error reporting function
            2010.12.20 adapted to function tbg_ifrqs() (filter problem)
            2011.03.20 optional integer transaction weights added
            2011.07.08 adapted to modified function tbg_recode()
            2011.08.28 output of item set counters per size added
            2011.11.23 output/filtering of maximal item sets added
            2012.04.29 alternative filtering of maximal item sets added
            2012.06.13 bug resulting from exchange of n and k fixed
            2012.07.09 pruning of the prefix tree made optional (-p)
            2012.07.16 optional processing with patricia tree added
            2013.04.01 adapted to type changes in module tract
            2013.10.18 optional pattern spectrum collection added
            2013.11.12 item selection file changed to option -R#
            2014.05.12 option -F# added (support border for filtering)
            2014.08.24 adapted to modified item set reporter interface
            2014.08.28 functions ista_data() and ista_report() added
            2014.10.24 changed from LGPL license to MIT license
            2016.02.19 added pre-formatting for some integer numbers
            2017.03.24 ista miner object and interface introduced
            2017.05.30 optional output compression with zlib added
            2017.06.13 bug in reporting mode fixed (ISR_NOFILTER)
------------------------------------------------------------------------
  Reference for the IsTa algorithm:
    C. Borgelt, X. Yang, R. Nogales-Cadenas,
    P. Carmona-Saez, and A. Pascual-Montano.
    Finding Closed Frequent Item Sets by Intersecting Transactions.
    Proc. 14th Int. Conf. on Extending Database Technology
    (EDBT 2011, Uppsala, Sweden), 367-376.
    ACM Press, New York, NY, USA 2011
  Reference for the basic idea:
    T. Mielik"ainen.
    Intersecting Data to Closed Sets with Constraints.
    Proc. Workshop Frequent Item Set Mining Implementations
    (FIMI 2003, Melbourne, FL).
    CEUR Workshop Proceedings 90, Aachen, Germany 2003
    http://www.ceur-ws.org/Vol-90/
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#ifndef ISR_PATSPEC
#define ISR_PATSPEC
#endif
#ifndef ISR_CLOMAX
#define ISR_CLOMAX
#endif
#ifdef  ISTA_MAIN
#ifndef PSP_REPORT
#define PSP_REPORT
#endif
#ifndef TA_READ
#define TA_READ
#endif
#endif
#ifdef ISTA_ABORT
#include "sigint.h"
#endif
#include "ista.h"
#include "pfxtree.h"
#include "pattree.h"
#ifdef ISTA_MAIN
#include "error.h"
#endif
#ifdef STORAGE
#include "storage.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define PRGNAME     "ista"
#define DESCRIPTION "find closed/maximal frequent item sets " \
                    "by intersecting transactions"
#define VERSION     "version 4.21 (2017.06.13)        " \
                    "(c) 2009-2017   Christian Borgelt"

/* --- error codes --- */
/* error codes   0 to  -4 defined in tract.h */
#define E_STDIN      (-5)       /* double assignment of stdin */
#define E_OPTION     (-6)       /* unknown option */
#define E_OPTARG     (-7)       /* missing option argument */
#define E_ARGCNT     (-8)       /* too few/many arguments */
#define E_TARGET     (-9)       /* invalid target type */
#define E_SIZE      (-10)       /* invalid item set size */
#define E_SUPPORT   (-11)       /* invalid item set support */
#define E_MEASURE   (-13)       /* invalid evaluation measure */
#define E_VARIANT   (-14)       /* invalid algorithm variant */
/* error codes -15 to -25 defined in tract.h */

#ifndef QUIET                   /* if not quiet version, */
#define MSG         fprintf     /* print messages */
#define XMSG        if (ista->mode & ISTA_VERBOSE) fprintf
#define CLOCK(t)    ((t) = clock())
#else                           /* if quiet version, */
#define MSG(...)    ((void)0)   /* suppress messages */
#define XMSG(...)   ((void)0)
#define CLOCK(t)    ((void)0)
#endif

#define SEC_SINCE(t)  ((double)(clock()-(t)) /(double)CLOCKS_PER_SEC)

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
struct _ista {                  /* --- ista miner --- */
  int      target;              /* target type (e.g. closed/maximal) */
  double   smin;                /* minimum support of an item set */
  double   smax;                /* maximum support of an item set */
  SUPP     supp;                /* minimum support of an item set */
  ITEM     zmin;                /* minimum size of a rule/item set */
  ITEM     zmax;                /* maximum size of a rule/item set */
  int      eval;                /* additional evaluation measure */
  double   thresh;              /* threshold for evaluation measure */
  int      algo;                /* variant of ista algorithm */
  int      mode;                /* search mode (e.g. pruning) */
  TABAG    *tabag;              /* transaction bag/multiset */
  ISREPORT *report;             /* item set reporter */
  PFXTREE  *pxt;                /* prefix tree for intersections */
  PATTREE  *pat;                /* patricia tree for intersections */
  SUPP     *frqs;               /* (remaining) item frequencies */
};                              /* (ista miner) */

/*----------------------------------------------------------------------
  Constants
----------------------------------------------------------------------*/
#if !defined QUIET && defined ISTA_MAIN
/* --- error messages --- */
static const char *errmsgs[] = {
  /* E_NONE      0 */  "no error",
  /* E_NOMEM    -1 */  "not enough memory",
  /* E_FOPEN    -2 */  "cannot open file %s",
  /* E_FREAD    -3 */  "read error on file %s",
  /* E_FWRITE   -4 */  "write error on file %s",
  /* E_STDIN    -5 */  "double assignment of standard input",
  /* E_OPTION   -6 */  "unknown option -%c",
  /* E_OPTARG   -7 */  "missing option argument",
  /* E_ARGCNT   -8 */  "wrong number of arguments",
  /* E_TARGET   -9 */  "invalid target type '%c'",
  /* E_SIZE    -10 */  "invalid item set size %"ITEM_FMT,
  /* E_SUPPORT -11 */  "invalid minimum support %g",
  /*           -12 */  NULL,
  /* E_MEASURE -13 */  "invalid evaluation measure '%c'",
  /* E_VARIANT -14 */  "invalid IsTa variant '%c'",
  /* E_NOITEMS -15 */  "no (frequent) items found",
  /*           -16 */  "unknown error"
};
#endif

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
#ifdef ISTA_MAIN
#ifndef QUIET
static CCHAR    *prgname;       /* program name for error messages */
#endif
static TABREAD  *tread  = NULL; /* table/transaction reader */
static ITEMBASE *ibase  = NULL; /* item base */
static TABAG    *tabag  = NULL; /* transaction bag/multiset */
static ISREPORT *report = NULL; /* item set reporter */
static TABWRITE *twrite = NULL; /* table writer for pattern spectrum */
static double   *border = NULL; /* support border for filtering */
static ISTA     *ista   = NULL; /* ista miner object */
#endif

/*----------------------------------------------------------------------
  Intersecting Transactions
----------------------------------------------------------------------*/

ISTA* ista_create (int target, double smin, double smax,
                   ITEM zmin, ITEM zmax, int eval, double thresh,
                   int algo, int mode)
{                               /* --- create an ista miner */
  ISTA *ista;                   /* created ista miner */

  /* --- make parameters consistent --- */
  if (target & ISTA_MAXIMAL) target = ISR_MAXIMAL;
  else                       target = ISR_CLOSED;

  /* --- create an ista miner --- */
  ista = (ISTA*)malloc(sizeof(ISTA));
  if (!ista) return NULL;       /* create an ista miner */
  ista->target = target;        /* and store all parameters */
  ista->smin   = smin;
  ista->smax   = smax;
  ista->supp   = 1;
  ista->zmin   = zmin;
  ista->zmax   = zmax;
  ista->eval   = eval;
  ista->thresh = thresh/100.0;
  ista->algo   = algo;
  ista->mode   = mode;
  ista->tabag  = NULL;
  ista->report = NULL;
  ista->pxt    = NULL;
  ista->pat    = NULL;
  return ista;                  /* return the created ista miner */
}  /* ista_create() */

/*--------------------------------------------------------------------*/

static int cleanup (ISTA *ista)
{                               /* --- clean up on error */
  if (ista->mode & ISTA_NOCLEAN)
    return E_NOMEM;             /* if not to clean up memory, abort */
  if (ista->pxt) {              /* free prefix tree */
    pxt_delete(ista->pxt, 1); ista->pxt = NULL; }
  if (ista->pat) {              /* free patricia tree */
    pat_delete(ista->pat);    ista->pat = NULL; }
  if (ista->frqs) {             /* free item frequencies */
    free(ista->frqs);         ista->frqs = NULL; }
  return E_NOMEM;               /* return an error indicator */
}  /* cleanup() */

/*--------------------------------------------------------------------*/

void ista_delete (ISTA *ista, int deldar)
{                               /* --- delete an ista miner */
  cleanup(ista);                /* clean up temporary data */
  if (deldar) {                 /* if to delete data and reporter */
    if (ista->report) isr_delete(ista->report, 0);
    if (ista->tabag)  tbg_delete(ista->tabag,  1);
  }                             /* delete if existing */
  free(ista);                   /* delete the base structure */
}  /* ista_delete() */

/*--------------------------------------------------------------------*/

int ista_data (ISTA *ista, TABAG *tabag, int sort)
{                               /* --- prepare data for IsTa */
  ITEM    m;                    /* number of items */
  double  smin;                 /* absolute minimum support */
  SUPP    w;                    /* total transaction weight */
  #ifndef QUIET                 /* if to print messages */
  TID     n;                    /* number of transactions */
  clock_t t;                    /* timer for measurements */
  #endif                        /* (only needed for messages) */

  assert(ista && tabag);        /* check the function arguments */
  ista->tabag = tabag;          /* note the transaction bag */

  /* --- compute data-specific parameters --- */
  w = tbg_wgt(tabag);           /* compute absolute minimum support */
  smin = ceilsupp((ista->smin < 0) ? -ista->smin
                : (ista->smin/100.0) *(double)w *(1-DBL_EPSILON));
  ista->supp = (SUPP)smin;
  if (ista->algo == ISTA_AUTO)  /* if automatic variant choice, */
    ista->algo = ISTA_PREFIX;   /* fix to prefix tree */

  /* --- sort and recode items --- */
  CLOCK(t);                     /* start timer, print log message */
  XMSG(stderr, "filtering, sorting and recoding items ... ");
  m = tbg_recode(tabag, ista->supp, -1, -1, -sort);
  if (m < 0) return E_NOMEM;    /* recode items and transactions */
  if (m < 1) return E_NOITEMS;  /* and check the number of items */
  XMSG(stderr, "[%"ITEM_FMT" item(s)]", m);
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));

  /* --- filter and sort transactions --- */
  CLOCK(t);                     /* start timer, print log message */
  XMSG(stderr, "filtering and sorting transactions ... ");
  tbg_filter(tabag, ista->zmin, NULL, 0);
  tbg_itsort(tabag, -1, 0);   /* remove items of short transactions */
  tbg_sortsz(tabag, -1, 0);   /* and sort items and transactions */
  /* The sorting direction is inverted here, because the transaction */
  /* identifiers are traversed backwards in both algorithm variants. */
  tbg_reduce(tabag, 0);         /* reduce transactions to unique ones */
  #ifndef QUIET                 /* if to print messages */
  n = tbg_cnt(tabag);           /* get the number of transactions */
  w = tbg_wgt(tabag);           /* and the transaction weight */
  XMSG(stderr, "[%"TID_FMT, n); /* print number of transactions */
  if (w != (SUPP)n) { XMSG(stderr, "/%"SUPP_FMT, w); }
  XMSG(stderr, " transaction(s)] done [%.2fs].\n", SEC_SINCE(t));
  #endif
  return 0;                     /* return 'ok' */
}  /* ista_data() */

/*--------------------------------------------------------------------*/

int ista_report (ISTA *ista, ISREPORT *report)
{                               /* --- prepare reporter for IsTa */
  TID    n;                     /* number of transactions */
  SUPP   w;                     /* total transaction weight */
  double smax;                  /* absolute maximum support */
  int    mrep;                  /* mode for item set reporter */

  assert(ista && report);       /* check the function arguments */
  ista->report = report;        /* note the item set reporter */

  /* --- get reporting mode --- */
  mrep = 0;                     /* initialize reporting mode */
  if ((ista->target & ISR_MAXIMAL) && !(ista->mode & ISTA_FILTER))
       mrep |= ISR_MAXIMAL;     /* build reporting mode */
  else mrep |= ISR_NOFILTER;    /* (no filtering if possible) */
  #ifdef USE_ZLIB               /* if optional output compression */
  if (ista->mode & ISTA_ZLIB)   /* if the compression flag is set, */
    mrep |= ISR_ZLIB;           /* transfer it to the report mode */
  #endif

  /* --- configure item set reporter --- */
  w = tbg_wgt(ista->tabag);     /* set support and size range */
  smax = (ista->smax < 0) ? -ista->smax
       : (ista->smax/100.0) *(double)w *(1-DBL_EPSILON);
  isr_setsupp(report, (RSUPP)ista->supp, (RSUPP)floorsupp(smax));
  isr_setsize(report, ista->zmin, ista->zmax);
  if (ista->eval == ISTA_LDRATIO)  /* set add. evaluation function */
    isr_seteval(report, isr_logrto, NULL, +1, ista->thresh);
  n = (ista->mode & ISTA_PREFMT)/* get range of numbers to preformat */
    ? (TID)ib_maxfrq(tbg_base(ista->tabag)) : -1;
  if ((isr_prefmt(report, (TID)ista->supp, n)      != 0)
  ||  (isr_settarg(report, ista->target, mrep, -1) != 0))
    return E_NOMEM;             /* set pre-format and target type */
  return 0;                     /* return 'ok' */
}  /* ista_report() */

/*--------------------------------------------------------------------*/

int ista_mine (ISTA *ista)
{                               /* --- intersecting transactions */
  int      r;                   /* result of function call, buffer */
  ITEM     m, k, z;             /* number of items, buffers */
  TID      n;                   /* number of transactions */
  SUPP     w;                   /* total transaction weight */
  TRACT    *tract;              /* to traverse the transactions */
  const ITEM *items;            /* to access the transaction items */
  const SUPP *ifs;              /* to access the item frequencies */
  #ifndef QUIET                 /* if to print messages */
  size_t   zc, zm;              /* number of tree nodes */
  clock_t  t;                   /* timer for measurements */
  #endif                        /* (only needed for messages) */

  assert(ista);                 /* check the function arguments */

  /* --- find closed/maximal frequent item sets --- */
  CLOCK(t);                     /* start timer, print log message */
  XMSG(stderr, "intersecting transactions ... ");
  ifs = tbg_ifrqs(ista->tabag, 0); /* get the item frequencies */
  if (!ifs)  return E_NOMEM;       /* in the transaction bag */
  m = tbg_itemcnt(ista->tabag);    /* get the number of items */
  ista->frqs = (SUPP*)malloc((size_t)m *sizeof(SUPP));
  if (!ista->frqs) return E_NOMEM; /* copy the item frequencies */
  memcpy(ista->frqs, ifs,    (size_t)m *sizeof(SUPP));
  if (ista->algo == ISTA_PATRICIA) {  /* if to use a patricia tree */
    ista->pat = pat_create(m, -1);    /* create a patricia tree */
    if (!ista->pat) return cleanup(ista); }
  else {                        /* if to use a prefix tree */
    ista->pxt = pxt_create(m, -1, NULL);
    if (!ista->pxt) return cleanup(ista);
  }                             /* create a prefix tree */
  n = tbg_cnt(ista->tabag);     /* get the number of transactions */
  for (k = 0; --n >= 0; ) {     /* traverse the transactions */
    #ifdef ISTA_ABORT           /* if to check for interrupt */
    if (sig_aborted()) break;   /* if execution was aborted, */
    #endif                      /* abort the loop */
    tract = tbg_tract(ista->tabag, n);
    items = ta_items(tract); z = ta_size(tract); w = ta_wgt(tract);
    r = (ista->pat)             /* intersect transaction and tree */
      ? pat_isect(ista->pat, items, z, w, ista->supp, ista->frqs)
      : pxt_isect(ista->pxt, items, z, w, ista->supp, ista->frqs);
    if (r < 0) return cleanup(ista);
    while (*items >= 0)         /* count newly prunable items */
      if ((ista->frqs[*items++] -= w) < ista->supp) k++;
    if ((ista->mode & ISTA_PRUNE) /* if to prune prefix/patricia tree */
    &&  (ista->supp >= 4) && (k > 0) && ((n & 0x0f) == 0x0f)) {
      r = (ista->pat) ? pat_prunex(ista->pat, ista->supp, ista->frqs)
                      : pxt_prunex(ista->pxt, ista->supp, ista->frqs);
      if (r < 0) return cleanup(ista); /* if item sets can be pruned, */
      k = 0;                    /* prune the prefix/patricia tree and */
    }                           /* clear the prunable item counter */
    #ifndef QUIET               /* if to print messages */
    if (((n & 0xff) == 0)       /* print number of rem. transactions */
    ||  ((n < 0xff) && ((n & 0x0f) == 0)) || (n <= 0x0f))
      XMSG(stderr, "%12"TID_FMT"\b\b\b\b\b\b\b\b\b\b\b\b", n);
    #endif
  }
  free(ista->frqs);             /* delete the item frequency array */
  ista->frqs = NULL;
  #ifndef QUIET                 /* if to print messages */
  if (ista->pat) { zc = pat_nodecnt(ista->pat);
                   zm = pat_nodemax(ista->pat); }
  else           { zc = pxt_nodecnt(ista->pxt);
                   zm = pxt_nodemax(ista->pxt); }
  XMSG(stderr, "[%"SIZE_FMT"/%"SIZE_FMT, zc, zm);
  XMSG(stderr, " node(s)] done [%.2fs].\n", SEC_SINCE(t));
  #endif
  #ifdef ISTA_ABORT             /* if to check for interrupt */
  if (sig_aborted()) { cleanup(ista); return -1; }
  #endif                        /* abort the function if requested */

  /* --- write closed/maximal frequent item sets --- */
  r = (ista->target & ISTA_MAXIMAL) ? 1 : 0;
  if (ista->mode & ISTA_FILTER)
    r = -r;                     /* get the maximal item set mode */
  if ((ista->mode & ISTA_PRUNE)
  && (r < 0)) {                 /* if to filter with repository */
    CLOCK(t);                   /* start timer, print log message */
    XMSG(stderr, "pruning item set repository ... ");
    if (ista->pat) pat_prune(ista->pat, ista->supp);
    else           pxt_prune(ista->pxt, ista->supp);
    #ifndef QUIET               /* if to print messages */
    if (ista->pat) { zc = pat_nodecnt(ista->pat);
                     zm = pat_nodemax(ista->pat); }
    else           { zc = pxt_nodecnt(ista->pxt);
                     zm = pxt_nodemax(ista->pxt); }
    XMSG(stderr, "[%"SIZE_FMT"/%"SIZE_FMT, zc, zm);
    XMSG(stderr, " node(s)] done [%.2fs].\n", SEC_SINCE(t));
    #endif                      /* prune the item set repository */
  }
  #ifdef ISTA_ABORT             /* if to check for interrupt */
  if (sig_aborted()) { cleanup(ista); return -1; }
  #endif                        /* abort the function if requested */
  CLOCK(t);                     /* start timer, print log message */
  XMSG(stderr, "writing %s ... ", isr_name(ista->report));
  k = (ista->pat)               /* report the found item sets */
    ? pat_report(ista->pat, r, ista->supp, ista->report)
    : pxt_report(ista->pxt, r, ista->supp, ista->report);
  if (k < 0) return cleanup(ista); /* check for a reporting error */
  XMSG(stderr, "[%"SIZE_FMT" set(s)]", isr_repcnt(ista->report));
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  cleanup(ista);                /* clean up allocated memory */
  return 0;                     /* return 'ok' */
}  /* ista_mine() */

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/
#ifdef ISTA_MAIN

static void help (void)
{                               /* --- print add. option information */
  #ifndef QUIET
  fprintf(stderr, "\n");        /* terminate startup message */
  printf("additional evaluation measures (option -e#)\n");
  printf("  x   no measure (default)\n");
  printf("  b   binary logarithm of support quotient\n");
  printf("\n");
  printf("information output format characters (option -v#)\n");
  printf("  %%%%  a percent sign\n");
  printf("  %%a  absolute item set support\n");
  printf("  %%s  relative item set support as a fraction\n");
  printf("  %%S  relative item set support as a percentage\n");
  printf("  %%e  additional evaluation measure\n");
  printf("  %%E  additional evaluation measure as a percentage\n");
  printf("  %%Q  total transaction weight (database size)\n");
  printf("All format characters can be preceded by the number\n");
  printf("of significant digits to be printed (at most 32 digits),\n");
  printf("even though this value is ignored for integer numbers.\n");
  #endif                        /* print help information */
  exit(0);                      /* abort the program */
}  /* help() */

/*--------------------------------------------------------------------*/

static ITEM getbdr (char *s, char **end, double **border)
{                               /* --- get the support border */
  ITEM   i, k;                  /* loop variables */
  double *b;                    /* support border */

  assert(s && end && border);   /* check the function arguments */
  for (i = k = 0; s[i]; i++)    /* traverse the string and */
    if (s[i] == ':') k++;       /* count the number separators */
  *border = b = (double*)malloc((size_t)++k *sizeof(double));
  if (!b) return -1;            /* allocate a support border */
  for (i = 0; i < k; i++) {     /* traverse the parameters */
    b[i] = strtod(s, end);      /* get the next parameter and */
    if (*end == s) break;       /* check for an empty parameter */
    s = *end; if (*s++ != ':') break;
  }                             /* check for a colon */
  if (++i < k)                  /* shrink support array if possible */
    *border = (double*)realloc(b, (size_t)i *sizeof(double));
  return i;                     /* return number of support values */
}  /* getbdr() */

/*--------------------------------------------------------------------*/

static int setbdr (ISREPORT *report, SUPP w, ITEM zmin,
                   double **border, ITEM n)
{                               /* --- set the support border */
  double s;                     /* to traverse the support values */

  assert(report                 /* check the function arguments */
  &&    (w > 0) && (zmin >= 0) && border && (*border || (n <= 0)));
  while (--n >= 0) {            /* traverse the support values */
    s = (*border)[n];           /* transform to absolute count */
    s = ceilsupp((s >= 0) ? s/100.0 *(double)w *(1-DBL_EPSILON) : -s);
    if (isr_setbdr(report, n+zmin, (RSUPP)s) < 0) return -1;
  }                             /* set support in item set reporter */
  if (*border) { free(*border); *border = NULL; }
  return 0;                     /* return 'ok' */
}  /* setbdr() */

/*--------------------------------------------------------------------*/

#ifndef NDEBUG                  /* if debug version */
  #undef  CLEANUP               /* clean up memory and close files */
  #define CLEANUP \
  if (ista)   ista_delete(ista,  0); \
  if (twrite) twr_delete(twrite, 1); \
  if (report) isr_delete(report, 0); \
  if (tabag)  tbg_delete(tabag,  0); \
  if (tread)  trd_delete(tread,  1); \
  if (ibase)  ib_delete (ibase);     \
  if (border) free(border);
#endif

GENERROR(error, exit)           /* generic error reporting function */

/*--------------------------------------------------------------------*/

int main (int argc, char *argv[])
{                               /* --- main function */
  int     i, k = 0;             /* loop variables, counters */
  char    *s;                   /* to traverse the options */
  CCHAR   **optarg = NULL;      /* option argument */
  CCHAR   *fn_inp  = NULL;      /* name of input  file */
  CCHAR   *fn_out  = NULL;      /* name of output file */
  CCHAR   *fn_sel  = NULL;      /* name of item selection file */
  CCHAR   *fn_psp  = NULL;      /* name of pattern spectrum file */
  CCHAR   *recseps = NULL;      /* record  separators */
  CCHAR   *fldseps = NULL;      /* field   separators */
  CCHAR   *blanks  = NULL;      /* blank   characters */
  CCHAR   *comment = NULL;      /* comment characters */
  CCHAR   *hdr     = "";        /* record header  for output */
  CCHAR   *sep     = " ";       /* item separator for output */
  CCHAR   *dflt    = " (%S)";   /* default format for check */
  CCHAR   *info    = dflt;      /* format for information output */
  int     target   = 'c';       /* target type (closed/maximal) */
  ITEM    zmin     = 1;         /* minimum size of an item set */
  ITEM    zmax     = ITEM_MAX;  /* maximum size of an item set */
  double  smin     = 10;        /* minimum support (in percent) */
  double  smax     = 100;       /* maximum support of an item set */
  int     eval     = 'x';       /* additional evaluation measure */
  double  thresh   = 10;        /* threshold for evaluation measure */
  int     sort     = -2;        /* flag for item sorting and recoding */
  int     algo     = ISTA_PREFIX;  /* variant of IsTa algorithm */
  int     mode     = ISTA_DEFAULT|ISTA_PREFMT;   /* search mode */
  int     mtar     = 0;         /* mode for transaction reading */
  int     scan     = 0;         /* flag for scanable item output */
  int     bdrcnt   = 0;         /* number of support values in border */
  int     stats    = 0;         /* flag for item set statistics */
  PATSPEC *psp;                 /* collected pattern spectrum */
  ITEM    m;                    /* number of items */
  TID     n;                    /* number of transactions */
  SUPP    w;                    /* total transaction weight */
  #ifndef QUIET                 /* if not quiet version */
  clock_t t;                    /* timer for measurements */

  prgname = argv[0];            /* get program name for error msgs. */

  /* --- print usage message --- */
  if (argc > 1) {               /* if arguments are given */
    fprintf(stderr, "%s - %s\n", argv[0], DESCRIPTION);
    fprintf(stderr, VERSION); } /* print a startup message */
  else {                        /* if no arguments given */
    printf("usage: %s [options] infile [outfile]\n", argv[0]);
    printf("%s\n", DESCRIPTION);
    printf("%s\n", VERSION);
    printf("-t#      target type                              "
                    "(default: %c)\n", target);
    printf("         (c: closed item sets, m: maximal item sets)\n");
    printf("-m#      minimum number of items per item set     "
                    "(default: %"ITEM_FMT")\n", zmin);
    printf("-n#      maximum number of items per item set     "
                    "(default: no limit)\n");
    printf("-s#      minimum support of an item set           "
                    "(default: %g%%)\n", smin);
    printf("-S#      maximum support of an item set/rule      "
                    "(default: %g%%)\n", smax);
    printf("         (positive: percentage, "
                     "negative: absolute number)\n");
    printf("-e#      additional evaluation measure            "
                    "(default: none)\n");
    printf("-d#      threshold for add. evaluation measure    "
                    "(default: %g%%)\n", thresh);
    printf("-q#      sort items w.r.t. their frequency        "
                    "(default: %d)\n", sort);
    printf("         (1: ascending, -1: descending, 0: do not sort,\n"
           "          2: ascending, -2: descending w.r.t. "
                    "transaction size sum)\n");
    printf("-i       use a patricia tree (or patricia trie)   ");
    printf(         "(default: prefix)\n");
    printf("         (may be faster for very few transactions "
                    "and very many items)\n");
    printf("-p       do not prune the prefix/patricia tree    "
                    "(default: prune)\n");
    printf("-j       filter maximal item sets with repository "
                    "(default: extra)\n");
    printf("         (needs less memory, but is usually slower)\n");
    printf("-F#:#..  support border for filtering item sets   "
                    "(default: none)\n");
    printf("         (list of minimum support values, "
                    "one per item set size,\n");
    printf("         starting at the minimum size, "
                    "as given with option -m#)\n");
    printf("-R#      read an item selection from a file\n");
    printf("-P#      write a pattern spectrum to a file\n");
    printf("-Z       print item set statistics "
                    "(number of item sets per size)\n");
    printf("-N       do not pre-format some integer numbers   "
                    "(default: do)\n");
    printf("-g       write output in scanable form "
                    "(quote certain characters)\n");
    #ifdef USE_ZLIB             /* if optional output compression */
    printf("-z       compress output with zlib (deflate)      "
                    "(default: plain text)\n");
    #endif                      /* print compression option */
    printf("-h#      record header  for output                "
                    "(default: \"%s\")\n", hdr);
    printf("-k#      item separator for output                "
                    "(default: \"%s\")\n", sep);
    printf("-v#      output format for item set information   "
                    "(default: \"%s\")\n", info);
    printf("-w       integer transaction weight in last field "
                    "(default: only items)\n");
    printf("-r#      record/transaction separators            "
                    "(default: \"\\n\")\n");
    printf("-f#      field /item        separators            "
                    "(default: \" \\t,\")\n");
    printf("-b#      blank   characters                       "
                    "(default: \" \\t\\r\")\n");
    printf("-C#      comment characters                       "
                    "(default: \"#\")\n");
    printf("-!       print additional option information\n");
    printf("infile   file to read transactions from           "
                    "[required]\n");
    printf("outfile  file to write frequent item sets to      "
                    "[optional]\n");
    return 0;                   /* print a usage message */
  }                             /* and abort the program */
  #endif  /* #ifndef QUIET */
  /* free option characters: aclouxy [A-Z]\[CFNPRZ] */

  /* --- evaluate arguments --- */
  for (i = 1; i < argc; i++) {  /* traverse arguments */
    s = argv[i];                /* get option argument */
    if (optarg) { *optarg = s; optarg = NULL; continue; }
    if ((*s == '-') && *++s) {  /* -- if argument is an option */
      while (*s) {              /* traverse options */
        switch (*s++) {         /* evaluate switches */
          case '!': help();                          break;
          case 't': target = (*s) ? *s++ : 'c';      break;
          case 'm': zmin   = (ITEM)strtol(s, &s, 0); break;
          case 'n': zmax   = (ITEM)strtol(s, &s, 0); break;
          case 's': smin   =       strtod(s, &s);    break;
          case 'S': smax   =       strtod(s, &s);    break;
          case 'e': eval   = (*s) ? *s++ : 0;        break;
          case 'd': thresh =       strtod(s, &s);    break;
          case 'q': sort   = (int) strtol(s, &s, 0); break;
          case 'i': algo   =  ISTA_PATRICIA;         break;
          case 'p': mode  &= ~ISTA_PRUNE;            break;
          case 'j': mode  |=  ISTA_FILTER;           break;
          case 'F': bdrcnt = getbdr(s, &s, &border); break;
          case 'R': optarg = &fn_sel;                break;
          case 'P': optarg = &fn_psp;                break;
          case 'Z': stats  = 1;                      break;
          case 'N': mode  &= ~ISTA_PREFMT;           break;
          case 'g': scan   = 1;                      break;
          #ifdef USE_ZLIB       /* if optional output compression */
          case 'z': mode  |= ISTA_ZLIB;              break;
          #endif                /* set the compression flag */
          case 'h': optarg = &hdr;                   break;
          case 'k': optarg = &sep;                   break;
          case 'v': optarg = &info;                  break;
          case 'w': mtar  |= TA_WEIGHT;              break;
          case 'r': optarg = &recseps;               break;
          case 'f': optarg = &fldseps;               break;
          case 'b': optarg = &blanks;                break;
          case 'C': optarg = &comment;               break;
          default : error(E_OPTION, *--s);           break;
        }                       /* set option variables */
        if (optarg && *s) { *optarg = s; optarg = NULL; break; }
      } }                       /* get option argument */
    else {                      /* -- if argument is no option */
      switch (k++) {            /* evaluate non-options */
        case  0: fn_inp = s;      break;
        case  1: fn_out = s;      break;
        default: error(E_ARGCNT); break;
      }                         /* note filenames */
    }
  }
  if (optarg)       error(E_OPTARG);     /* check option arguments */
  if (k      < 1)   error(E_ARGCNT);     /* and number of arguments */
  if (zmin   < 0)   error(E_SIZE, zmin); /* check the size limits */
  if (zmax   < 0)   error(E_SIZE, zmax); /* and the minimum support */
  if (smin   > 100) error(E_SUPPORT, smin);
  if (bdrcnt < 0)   error(E_NOMEM);
  switch (target) {             /* check and translate target type */
    case 'c': target = ISTA_CLOSED;          break;
    case 'm': target = ISTA_MAXIMAL;         break;
    default : error(E_TARGET, (char)target); break;
  }                             /* (get target type code) */
  switch (eval) {               /* check and translate measure */
    case 'x': eval = ISTA_NONE;              break;
    case 'b': eval = ISTA_LDRATIO;           break;
    default : error(E_MEASURE, (char)eval);  break;
  }                             /* (get evaluation measure code) */
  if (info == dflt)             /* adapt the default info. format */
    info = (smin < 0) ? " (%a)" : " (%S)";
  mode |= ISTA_VERBOSE|ISTA_NOCLEAN;
  MSG(stderr, "\n");            /* terminate the startup message */

  /* --- read item selection --- */
  ibase = ib_create(0, 0);      /* create an item base */
  if (!ibase) error(E_NOMEM);   /* to manage the items */
  tread = trd_create();         /* create a transaction reader */
  if (!tread) error(E_NOMEM);   /* and configure the characters */
  trd_allchs(tread, recseps, fldseps, blanks, "", comment);
  if (fn_sel) {                 /* if item appearances are given */
    CLOCK(t);                   /* start timer, open input file */
    if (trd_open(tread, NULL, fn_sel) != 0)
      error(E_FOPEN, trd_name(tread));
    MSG(stderr, "reading %s ... ", trd_name(tread));
    m = ib_readsel(ibase,tread);/* read the given item selection */
    if (m < 0) error((int)-m, ib_errmsg(ibase, NULL, 0));
    trd_close(tread);           /* close the input file */
    MSG(stderr, "[%"ITEM_FMT" item(s)]", m);
    MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  }                             /* print a log message */

  /* --- read transaction database --- */
  tabag = tbg_create(ibase);    /* create a transaction bag */
  if (!tabag) error(E_NOMEM);   /* to store the transactions */
  CLOCK(t);                     /* start timer, open input file */
  if (trd_open(tread, NULL, fn_inp) != 0)
    error(E_FOPEN, trd_name(tread));
  MSG(stderr, "reading %s ... ", trd_name(tread));
  k = tbg_read(tabag, tread, mtar);
  if (k < 0) error(-k, tbg_errmsg(tabag, NULL, 0));
  trd_delete(tread, 1);         /* read the transaction database, */
  tread = NULL;                 /* then delete the table reader */
  m = ib_cnt(ibase);            /* get the number of items, */
  n = tbg_cnt(tabag);           /* the number of transactions, */
  w = tbg_wgt(tabag);           /* the total transaction weight */
  MSG(stderr, "[%"ITEM_FMT" item(s), %"TID_FMT, m, n);
  if (w != (SUPP)n) { MSG(stderr, "/%"SUPP_FMT, w); }
  MSG(stderr, " transaction(s)] done [%.2fs].", SEC_SINCE(t));
  if ((m <= 0) || (n <= 0))     /* check for at least one item */
    error(E_NOITEMS);           /* and at least one transaction */
  MSG(stderr, "\n");            /* terminate the log message */

  /* --- find maximal/closed item sets --- */
  ista = ista_create(target, smin, smax, zmin, zmax,
                     eval, thresh, algo, mode);
  if (!ista) error(E_NOMEM);    /* create an ista miner */
  k = ista_data(ista, tabag, sort);
  if (k) error(k);              /* prepare data for IsTa */
  report = isr_create(ibase);   /* create an item set reporter */
  if (!report) error(E_NOMEM);  /* and configure it */
  k = ista_report(ista, report);
  if (k) error(k);              /* prepare reporter for IsTa */
  if (setbdr(report, w, zmin, &border, bdrcnt) != 0)
    error(E_NOMEM);             /* set the support border */
  if (fn_psp && (isr_addpsp(report, NULL) < 0))
    error(E_NOMEM);             /* set a pattern spectrum if req. */
  if (isr_setfmt(report, scan, hdr, sep, NULL, info) != 0)
    error(E_NOMEM);             /* set the output format strings */
  k = isr_open(report, NULL, fn_out);
  if (k) error(k, isr_name(report));
  if (isr_setup(report) < 0)    /* open the item set file and */
    error(E_NOMEM);             /* set up the item set reporter */
  k = ista_mine(ista);          /* find frequent item sets */
  if (k) error(k);              /* with the IsTa algorithm */
  if (stats)                    /* print item set statistics */
    isr_prstats(report, stdout, 0);
  if (isr_close(report) != 0)   /* close item set output file */
    error(E_FWRITE, isr_name(report));

  /* --- write pattern spectrum --- */
  if (fn_psp) {                 /* if to write a pattern spectrum */
    CLOCK(t);                   /* start timer, create table write */
    psp    = isr_getpsp(report);/* get the pattern spectrum */
    twrite = twr_create();      /* create a table writer and */
    if (!twrite) error(E_NOMEM);/* open the output file */
    if (twr_open(twrite, NULL, fn_psp) != 0)
      error(E_FOPEN,  twr_name(twrite));
    MSG(stderr, "writing %s ... ", twr_name(twrite));
    if (psp_report(psp, twrite, 1.0) != 0)
      error(E_FWRITE, twr_name(twrite));
    twr_delete(twrite, 1);      /* write the pattern spectrum */
    twrite = NULL;              /* and delete the table writer */
    MSG(stderr, "[%"SIZE_FMT" signature(s)]", psp_sigcnt(psp));
    MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  }                             /* write a log message */

  /* --- clean up --- */
  CLEANUP;                      /* clean up memory and close files */
  SHOWMEM;                      /* show (final) memory usage */
  return 0;                     /* return 'ok' */
}  /* main() */

#endif  /* #ifdef ISTA_MAIN ... */
