/*----------------------------------------------------------------------
  File    : apriori.c
  Contents: apriori algorithm for finding association rules
  Author  : Christian Borgelt
  History : 1996.02.14 file created
            1996.07.26 output precision reduced
            1996.11.22 options -b, -f, and -r added
            1996.11.24 option -e added (add. evaluation measures)
            1997.08.18 normalized chi^2 measure and option -m added
            1997.10.13 quiet version (no output to stdout or stderr)
            1998.01.27 adapted to changed ist_create() function
            1998.08.08 optional input file (item appearances) added
            1998.09.07 hyperedge mode added (option -h)
            1998.12.08 output of absolute support added (option -a)
            1998.12.09 conversion of names to a scanable form added
            1999.02.09 input from stdin, output to stdout added
            1999.08.09 bug in check of support parameter fixed (<= 0)
            1999.11.05 rule evaluation measure IST_LIFTDIFF added
            1999.11.08 output of add. rule eval. measure value added
            2000.03.16 optional use of original rule support definition
            2001.04.01 option -h replaced by option -t (target type)
            2001.05.26 extended support output added
            2001.08.15 module scan used for item name formatting
            2001.11.18 item and transaction functions made a module
            2001.11.19 options -C, -l changed, option -y removed
            2001.12.28 adapted to module tract, some improvements
            2002.01.11 evaluation measures codes changed to letters
            2002.02.10 option -q extended by a direction parameter
            2002.02.11 memory usage minimization option added
            2002.06.09 arbitrary supp./conf. formats made possible
            2003.01.09 option -k added (user-defined item separator)
            2003.01.14 check for empty transaction set added
            2003.07.17 item filtering w.r.t. usage added (option -u)
            2003.07.17 sorting w.r.t. transaction size sum added
            2003.07.18 maximal item set filter added
            2003.08.11 closed  item set filter added
            2003.08.15 item filtering for transaction tree added
            2003.08.16 parameter for transaction filtering added
            2003.08.18 dynamic filtering decision based on times added
            2004.03.25 option -S added (maximal support of a set/rule)
            2004.05.09 additional selection measures for sets added
            2004.10.28 two unnecessary assignments removed
            2004.11.23 absolute/relative support output changed
            2005.01.25 bug in output of absolute/relative support fixed
            2005.01.31 another bug in this output fixed
            2005.06.20 use of flag for "no item sorting" corrected
            2007.02.13 adapted to modified module tabread
            2008.08.12 option -l removed (do not load transactions)
            2008.08.14 adapted to extension of function tbg_filter()
            2008.08.18 more flexible output format control (numbers)
            2008.09.01 option -I added (implication sign for rules)
            2008.09.07 optional a-posteriori pruning added
            2008.09.09 more flexible information output control
            2008.09.10 item set extraction and evaluation redesigned
            2008.10.30 adapted to changes in item set reporting
            2008.11.13 adapted to changes in transaction management
            2008.11.19 adapted to modified transaction tree interface
            2008.12.06 perfect extension pruning added (optional)
            2009.05.28 adapted to modified function tbg_filter()
            2010.03.02 bug concerning size-sorted output fixed
            2010.03.04 several additional evaluation measures added
            2010.06.18 filtering for increase of evaluation added (-i)
            2010.07.14 output file made optional (for benchmarking)
            2010.08.22 adapted to modified modules tract and tabread
            2010.08.30 Fisher's exact test added as evaluation measure
            2010.10.15 adapted to modified interface of module report
            2010.10.22 chi^2 measure with Yates correction added
            2010.11.24 adapted to modified error reporting (tract)
            2010.12.11 adapted to a generic error reporting function
            2011.03.20 optional integer transaction weights added
            2011.07.08 adapted to modified function tbg_recode()
            2011.07.21 bug in function apriori fixed (closed/maximal)
            2011.07.25 threshold inverted for measures yielding p-values
            2011.08.04 weak forward pruning with evaluation added
            2011.08.16 new target type 'generators' added (option -tg)
            2011.08.28 output of item set counters per size added
            2011.10.18 all processing shifted to apriori function
            2012.01.23 rounding the mininum support for mining fixed
            2012.04.30 bug in apriori() fixed ("rule(s)" output)
            2012.06.13 bug in apriori() fixed (APR_INVBXS in eval)
            2012.11.29 bug in apriori() fixed (perf. ext. and ISR_RULE)
            2013.03.30 adapted to type changes in module tract
            2013.06.03 real-valued support made possible (SUPP==double)
            2013.10.18 optional pattern spectrum collection added
            2013.11.12 item selection file changed to option -R#
            2014.05.12 option -F# added (support border for filtering)
            2014.05.15 evaluation measures "cprob" and "import" added
            2014.05.31 bug with number of items in selection file fixed
            2014.07.03 missing terminal output of item set sizes added
            2014.08.18 preprocessing of transactions moved to apriori()
            2014.08.21 adapted to modified item set reporter interface
            2014.08.28 functions apriori_data(), apriori_report() added
            2014.10.24 changed from LGPL license to MIT license
            2016.02.19 added pre-formatting for some integer numbers
            2016.11.04 apriori miner object and interface introduced
            2017.05.30 optional output compression with zlib added
            2017.08.01 bug in calls to apriori_data() fixed (arg. order)
------------------------------------------------------------------------
  Reference for the Apriori algorithm:
    R. Agrawal and R. Srikant.
    Fast Algorithms for Mining Association Rules.
    Proc. 20th Int. Conf. on Very Large Databases
    (VLDB 1994, Santiago de Chile), 487-499.
    Morgan Kaufmann, San Mateo, CA, USA 1994
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#ifndef ISR_PATSPEC
#define ISR_PATSPEC
#endif
#ifdef  APR_MAIN
#ifndef PSP_REPORT
#define PSP_REPORT
#endif
#ifndef TA_READ
#define TA_READ
#endif
#endif
#ifndef TATREEFN
#define TATREEFN
#endif
#ifdef APR_ABORT
#include "sigint.h"
#endif
#ifndef TATREEFN
#define TATREEFN
#endif
#include "apriori.h"
#ifdef APR_MAIN
#include "error.h"
#endif
#ifdef STORAGE
#include "storage.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define PRGNAME     "apriori"
#define DESCRIPTION "find frequent item sets with the apriori algorithm"
#define VERSION     "version 6.28 (2019.03.31)        " \
                    "(c) 1996-2019   Christian Borgelt"

/* --- error codes --- */
/* error codes   0 to  -4 defined in tract.h */
#define E_STDIN      (-5)       /* double assignment of stdin */
#define E_OPTION     (-6)       /* unknown option */
#define E_OPTARG     (-7)       /* missing option argument */
#define E_ARGCNT     (-8)       /* too few/many arguments */
#define E_TARGET     (-9)       /* invalid target type */
#define E_SIZE      (-10)       /* invalid set/rule size */
#define E_SUPPORT   (-11)       /* invalid support */
#define E_CONF      (-12)       /* invalid confidence */
#define E_MEASURE   (-13)       /* invalid evaluation measure */
#define E_AGGMODE   (-14)       /* invalid aggregation mode */
#define E_STAT      (-16)       /* invalid test statistic */
#define E_SIGLVL    (-17)       /* invalid significance level */
/* error codes -15 to -25 defined in tract.h */

#ifndef QUIET                   /* if not quiet version, */
#define MSG         fprintf     /* print messages */
#define XMSG        if (apriori->mode & APR_VERBOSE) fprintf
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
struct _apriori {               /* --- apriori miner --- */
  int      target;              /* target type (e.g. closed/maximal) */
  double   smin;                /* minimum support of an item set */
  double   smax;                /* maximum support of an item set */
  SUPP     supp;                /* minimum support of an item set */
  SUPP     body;                /* minimum support of a rule body */
  double   conf;                /* minimum confidence of a rule */
  ITEM     zmin;                /* minimum size of a rule/item set */
  ITEM     zmax;                /* maximum size of a rule/item set */
  int      eval;                /* additional evaluation measure */
  int      agg;                 /* aggregation mode for eval. measure */
  double   thresh;              /* threshold for evaluation measure */
  int      algo;                /* variant of apriori algorithm */
  int      mode;                /* search mode (e.g. pruning) */
  TABAG    *tabag;              /* transaction bag/multiset */
  ISREPORT *report;             /* item set reporter */
  TATREE   *tatree;             /* transaction tree */
  ISTREE   *istree;             /* item set tree (for counting) */
  ITEM     *map;                /* identifier map for filtering */
};                              /* (apriori miner) */

/*----------------------------------------------------------------------
  Constants
----------------------------------------------------------------------*/
#if !defined QUIET && defined APR_MAIN
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
  /* E_SIZE    -10 */  "invalid item set or rule size %"ITEM_FMT,
  /* E_SUPPORT -11 */  "invalid minimum support %g",
  /* E_CONF    -12 */  "invalid minimum confidence %g",
  /* E_MEASURE -13 */  "invalid evaluation measure '%c'",
  /* E_AGGMODE -14 */  "invalid aggregation mode '%c'",
  /* E_NOITEMS -15 */  "no (frequent) items found",
  /* E_STAT    -16 */  "invalid test statistic '%c'",
  /* E_SIGLVL  -17 */  "invalid significance level/p-value %g",
  /*           -18 */  "unknown error"
};
#endif

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
#ifdef APR_MAIN
#ifndef QUIET
static CCHAR    *prgname;       /* program name for error messages */
#endif
static TABREAD  *tread  = NULL; /* table/transaction reader */
static ITEMBASE *ibase  = NULL; /* item base */
static TABAG    *tabag  = NULL; /* transaction bag/multiset */
static ISREPORT *report = NULL; /* item set reporter */
static TABWRITE *twrite = NULL; /* table writer for pattern spectrum */
static double   *border = NULL; /* support border for filtering */
static APRIORI *apriori = NULL; /* apriori miner object */
#endif

/*----------------------------------------------------------------------
  Apriori Algorithm
----------------------------------------------------------------------*/

APRIORI* apriori_create (int target, double smin, double smax,
                         double conf, ITEM zmin, ITEM zmax,
                         int eval, int agg, double thresh,
                         int algo, int mode)
{                               /* --- create an apriori miner */
  APRIORI *apriori;             /* created apriori miner */

  /* --- make parameters consistent --- */
  if      (target & APR_RULES)   target = ISR_RULES;
  else if (target & APR_GENERAS) target = ISR_GENERAS;
  else if (target & APR_MAXIMAL) target = ISR_MAXIMAL;
  else if (target & APR_CLOSED)  target = ISR_CLOSED;
  else                           target = ISR_FREQUENT;
  if (!(target & ISR_RULES))    /* if not to find association rules, */
    conf = 100.0;               /* set a neutral confidence */

  /* --- create an apriori miner --- */
  apriori = (APRIORI*)malloc(sizeof(APRIORI));
  if (!apriori) return NULL;    /* create an apriori miner */
  apriori->target = target;     /* store all parameters */
  apriori->smin   = smin;
  apriori->smax   = smax;
  apriori->supp   = 1;
  apriori->body   = 1;
  apriori->conf   = conf/100.0;
  apriori->zmin   = zmin;
  apriori->zmax   = zmax;
  apriori->eval   = eval;
  apriori->agg    = agg;
  apriori->thresh = thresh/100.0;
  apriori->algo   = algo;
  apriori->mode   = mode;
  apriori->tabag  = NULL;
  apriori->report = NULL;
  apriori->tatree = NULL;
  apriori->istree = NULL;
  apriori->map    = NULL;
  return apriori;               /* return the created apriori miner */
}  /* apriori_create() */

/*--------------------------------------------------------------------*/

static int cleanup (APRIORI *apriori)
{                               /* --- clean up on error */
  if (apriori->mode & APR_NOCLEAN)
    return E_NOMEM;             /* if not to clean up memory, abort */
  if (apriori->map) {           /* free identifier map for filtering */
    free(apriori->map);             apriori->map    = NULL; }
  if (apriori->istree) {        /* free item set tree (for counting) */
    ist_delete(apriori->istree);    apriori->istree = NULL; }
  if (apriori->tatree) {        /* free the transaction tree */
    tat_delete(apriori->tatree, 0); apriori->tatree = NULL; }
  return E_NOMEM;               /* return an error indicator */
}  /* cleanup() */

/*--------------------------------------------------------------------*/

void apriori_delete (APRIORI *apriori, int deldar)
{                               /* --- delete an apriori miner */
  cleanup(apriori);             /* clean up temporary data */
  if (deldar) {                 /* if to delete data and reporter */
    if (apriori->report) isr_delete(apriori->report, 0);
    if (apriori->tabag)  tbg_delete(apriori->tabag,  1);
  }                             /* delete if existing */
  free(apriori);                /* delete the base structure */
}  /* apriori_delete() */

/*--------------------------------------------------------------------*/

int apriori_data (APRIORI *apriori, TABAG *tabag, int mode, int sort)
{                               /* --- prepare data for Apriori */
  ITEM    m;                    /* number of items */
  double  smin;                 /* absolute minimum support */
  SUPP    w;                    /* total transaction weight */
  int     e;                    /* evaluation without flags */
  #ifndef QUIET                 /* if to print messages */
  TID     n;                    /* number of transactions */
  clock_t t;                    /* timer for measurements */
  #endif                        /* (only needed for messages) */

  assert(apriori && tabag);     /* check the function arguments */
  apriori->tabag = tabag;       /* note the transaction bag */

  /* --- compute data-specific parameters --- */
  w = tbg_wgt(tabag);           /* compute absolute minimum support */
  smin = ceilsupp((apriori->smin < 0) ? -apriori->smin
                : (apriori->smin/100.0) *(double)w *(1-DBL_EPSILON));
  apriori->body = (SUPP)smin;   /* compute body and body&head support */
  if ((apriori->target & ISR_RULES) && !(apriori->mode & APR_ORIGSUPP))
    smin *= apriori->conf *(1-DBL_EPSILON);
  apriori->supp = (SUPP)ceilsupp(smin);

  /* --- sort and recode items --- */
  if (!(mode & APR_NORECODE)) { /* if to sort and recode the items */
    CLOCK(t);                   /* start timer, print log message */
    XMSG(stderr, "filtering, sorting and recoding items ... ");
    m = tbg_recode(tabag, apriori->supp, -1, -1, sort);
    if (m < 0) return E_NOMEM;  /* recode items and transactions */
    if (m < 1) return E_NOITEMS;/* and check the number of items */
    XMSG(stderr, "[%"ITEM_FMT" item(s)]", m);
    XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  }                             /* print a log message */

  /* --- sort and reduce transactions --- */
  CLOCK(t);                     /* start timer, print log message */
  XMSG(stderr, "sorting and reducing transactions ... ");
  e = apriori->eval & ~APR_INVBXS;
  if (!(mode & APR_NOFILTER)    /* filter transactions if possible */
  &&  !(apriori->target & ISR_RULES)
  &&  ((e <= RE_NONE) || (e >= RE_FNCNT)))
    tbg_filter(tabag, apriori->zmin, NULL, 0);
  if (!(mode & APR_NOSORT)) {   /* if to sort items and transactions, */
    tbg_itsort(tabag, +1, 0);   /* sort items in transactions and */
    tbg_sort(tabag, +1, 0);     /* sort the trans. lexicographically */
    if (!(mode & APR_NOREDUCE)) /* if to combine equal transactions, */
      tbg_reduce(tabag, 0);     /* reduce transactions to unique ones */
  }                             /* (need sorting for reduction) */
  #ifndef QUIET                 /* if to print messages */
  n = tbg_cnt(tabag);           /* get the number of transactions */
  w = tbg_wgt(tabag);           /* and the transaction weight */
  XMSG(stderr, "[%"TID_FMT, n); /* print number of transactions */
  if (w != (SUPP)n) { XMSG(stderr, "/%"SUPP_FMT, w); }
  XMSG(stderr, " transaction(s)] done [%.2fs].\n", SEC_SINCE(t));
  #endif
  return 0;                     /* return 'ok' */
}  /* apriori_data() */

/*--------------------------------------------------------------------*/

int apriori_report (APRIORI *apriori, ISREPORT *report)
{                               /* --- prepare reporter for Apriori */
  TID    n;                     /* number of transactions */
  SUPP   w;                     /* total transaction weight */
  double smax;                  /* absolute maximum support */
  int    mrep;                  /* mode for item set reporter */

  assert(apriori && report);    /* check the function arguments */
  apriori->report = report;     /* note the item set reporter */

  /* --- get reporting mode --- */
  mrep = ISR_NOFILTER;          /* do not filter for closed/maximal */
  #ifdef USE_ZLIB               /* if optional output compression */
  if (apriori->mode & APR_ZLIB) /* if the compression flag is set, */
    mrep |= ISR_ZLIB;           /* transfer it to the report mode */
  #endif

  /* --- configure item set reporter --- */
  w = tbg_wgt(apriori->tabag);  /* set support and size range */
  smax = (apriori->smax < 0) ? -apriori->smax
       : (apriori->smax/100.0) *(double)w *(1-DBL_EPSILON);
  isr_setsupp(report, (RSUPP)apriori->supp, (RSUPP)floorsupp(smax));
  isr_setsize(report, apriori->zmin, apriori->zmax);
  if ((apriori->eval & ~APR_INVBXS) == APR_LDRATIO)
    isr_seteval(report, isr_logrto, NULL, +1, apriori->thresh);
  n = (apriori->mode & APR_PREFMT)/* get range of nums. to preformat */
    ? (TID)ib_maxfrq(tbg_base(apriori->tabag)) : -1;
  if ((isr_prefmt(report, (TID)apriori->supp, n)      != 0)
  ||  (isr_settarg(report, apriori->target, mrep, -1) != 0))
    return E_NOMEM;             /* set pre-format and target type */
  return 0;                     /* return 'ok' */
}  /* apriori_report() */

/*--------------------------------------------------------------------*/

int apriori_mine (APRIORI *apriori, ITEM prune, double filter,int order)
{                               /* --- apriori algorithm */
  ITEM    m, i, k;              /* number of items, loop variables */
  ITEM    size;                 /* number of items in set/rule */
  ITEM    xmax;                 /* maximum size for extensions */
  int     e, mode;              /* evaluation without flags, mode */
  clock_t t, tt, tc, x;         /* timers for measurements */

  assert(apriori);              /* check the function arguments */
  e = apriori->eval & ~APR_INVBXS; /* check and adapt evaluation */
  if (e <= RE_NONE) prune = ITEM_MIN;

  /* --- create transaction tree --- */
  tt = 0;                       /* init. the tree construction time */
  if (apriori->mode & APR_TATREE) { /* if to use a transaction tree */
    t = clock();                /* start the timer for construction */
    XMSG(stderr, "building transaction tree ... ");
    apriori->tatree = tat_create(apriori->tabag);
    if (!apriori->tatree)       /* create a transaction tree */
      return E_NOMEM;           /* as a compressed representation */
    XMSG(stderr, "[%"SIZE_FMT" node(s)]", tat_size(apriori->tatree));
    XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
    tt = clock() -t;            /* note the time for the construction */
  }                             /* of the transaction tree */

  /* --- create item set tree --- */
  if ((apriori->target & (ISR_CLOSED|ISR_MAXIMAL|ISR_RULES))
  ||  ((e > RE_NONE) && (e < RE_FNCNT)) || order)
    apriori->mode &= ~IST_PERFECT; /* remove perfect ext. pruning */
  t = clock(); tc = 0;          /* start the timer for the search */
  mode = apriori->mode & ~(IST_PARTIAL|IST_REVERSE);
  apriori->istree = ist_create(tbg_base(apriori->tabag), mode,
                         apriori->supp, apriori->body, apriori->conf);
  if (!apriori->istree) return cleanup(apriori);
  xmax = ((apriori->target & (ISR_CLOSED|ISR_MAXIMAL))
      && !(apriori->target & ISR_RULES)
      &&  (apriori->zmax   < ITEM_MAX))
       ? apriori->zmax+1 : apriori->zmax;
  m = tbg_max(apriori->tabag);  /* compute maximum extension size */
  if (xmax > m) xmax = m;       /* and limit it to transaction size */
  if (e == APR_LDRATIO)         /* set additional evaluation measure */
       isr_seteval(apriori->report, isr_logrto, NULL,
                   +1, apriori->thresh);
  else ist_seteval(apriori->istree, apriori->eval, apriori->agg,
                   apriori->thresh, prune);

  /* --- check item subsets --- */
  XMSG(stderr, "checking subsets of size 1");
  m = tbg_itemcnt(apriori->tabag); /* create an item map for pruning */
  apriori->map = (ITEM*)malloc((size_t)m *sizeof(ITEM));
  if (!apriori->map) return cleanup(apriori);
  for (i = m; 1; ) {            /* traverse the item set sizes */
    #ifdef APR_ABORT            /* if to check for interrupt: */
    if (sig_aborted()) break;   /* if execution was aborted, */
    #endif                      /* abort the loop */
    size = ist_height(apriori->istree);
    if (size >= xmax)           /* get the current item set size and */
      break;                    /* abort if maximal size is reached */
    if ((filter != 0)           /* if to filter w.r.t. item usage */
    && ((i = ist_check(apriori->istree, (int*)apriori->map)) <= size))
      break;                    /* check which items are still used */
    if (apriori->mode & APR_POST)    /* if a-posteriori pruning, */
      ist_prune(apriori->istree);    /* prune infrequent item sets */
    k = ist_addlvl(apriori->istree); /* add a level to the tree */
    if (k < 0) return cleanup(apriori);
    if (k > 0) break;           /* if no level was added, abort */
    if (((filter < 0)           /* if to filter w.r.t. item usage and */
    &&   ((double)i < -filter *(double)m))    /* enough items removed */
    ||  ((filter > 0) && (i < m)/* or counting time is long enough */
    && ((double)i*(double)tt < (double)filter*(double)m*(double)tc))) {
      m = i;                    /* note the new number of items */
      x = clock();              /* start the timer for filtering */
      if (apriori->tatree) {    /* if a transaction tree was created */
        if (tat_filter(apriori->tatree, size+1, (int*)apriori->map, 0))
          return cleanup(apriori); }   /* filter the transaction tree */
      else {                    /* if there is only a transaction bag */
        tbg_filter(apriori->tabag, size+1, (int*)apriori->map, 0);
        tbg_sort  (apriori->tabag, 0, 0); /* remove unnecessary items */
        tbg_reduce(apriori->tabag, 0);    /* and trans. and reduce */
      }                                   /* trans. to unique ones */
      tt = clock() -x;          /* note the filter/rebuild time */
    }
    size += 1;                  /* increment the item set size */
    XMSG(stderr, " %"ITEM_FMT, size);          /* and print it */
    x = clock();                /* start the timer for counting */
    if (apriori->tatree) ist_countx(apriori->istree, apriori->tatree);
    else                 ist_countb(apriori->istree, apriori->tabag);
    ist_commit(apriori->istree);/* count the transaction tree/bag */
    tc = clock() -x;            /* compute the new counting time */
  }
  free(apriori->map);           /* delete the filter map */
  apriori->map = NULL;          /* and the transaction tree */
  if (apriori->tatree && !(apriori->mode & APR_NOCLEAN)) {
    tat_delete(apriori->tatree, 0); apriori->tatree = NULL; }
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  #ifdef APR_ABORT              /* if to check for interrupt */
  if (sig_aborted()) { cleanup(apriori); return -1; }
  #endif                        /* abort the function if requested */

  /* --- filter found item sets --- */
  if ((prune >  ITEM_MIN)       /* if to filter with evaluation */
  &&  (prune <= 0)) {           /* (backward and weak forward) */
    CLOCK(t);                   /* start the timer for filtering */
    XMSG(stderr, "filtering with evaluation ... ");
    ist_filter(apriori->istree, prune);
    XMSG(stderr, "done [%.2fs].\n", SEC_SINCE(t));
  }                             /* mark non-qualifying sets */
  #ifdef APR_ABORT              /* if to check for interrupt */
  if (sig_aborted()) { cleanup(apriori); return -1; }
  #endif                        /* abort the function if requested */
  if (apriori->target & (ISR_CLOSED|ISR_MAXIMAL|ISR_GENERAS)) {
    CLOCK(t);                   /* start the timer for filtering */
    XMSG(stderr, "filtering for %s item sets ... ",
         (apriori->target & ISR_GENERAS) ? "generator" :
         (apriori->target & ISR_MAXIMAL) ? "maximal" : "closed");
    ist_clomax(apriori->istree, /* filter closed/maximal/generators */
               apriori->target | ((prune > ITEM_MIN) ? IST_SAFE : 0));
    XMSG(stderr, "done [%.2fs].\n", SEC_SINCE(t));
  }
  #ifdef APR_ABORT              /* if to check for interrupt */
  if (sig_aborted()) { cleanup(apriori); return -1; }
  #endif                        /* abort the function if requested */

  /* --- report item sets/association rules --- */
  CLOCK(t);                     /* start the output timer */
  XMSG(stderr, "writing %s ... ", isr_name(apriori->report));
  ist_init(apriori->istree, order); /* initialize the extraction */
  if (ist_report(apriori->istree, apriori->report, apriori->target) < 0)
    return cleanup(apriori);    /* report item sets/association rules */
  XMSG(stderr, "[%"SIZE_FMT" %s(s)]", isr_repcnt(apriori->report),
               (apriori->target == ISR_RULES) ? "rule" : "set");
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  #ifdef BENCH                  /* if benchmark version, */
  ist_stats(apriori->istree);   /* show the search statistics */
  #endif                        /* (especially memory usage) */
  cleanup(apriori);             /* clean up the allocated memory */
  return 0;                     /* return 'ok' */
}  /* apriori_mine() */

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/
#ifdef APR_MAIN
#ifdef APRIACC

static void help (void)
{                               /* --- print add. option information */
  #ifndef QUIET
  fprintf(stderr, "\n");        /* terminate startup message */
  printf("test statistics for p-value computation (option -t#)\n");
  printf("  x      no statistic / zero\n");
  printf("  c/p/n  chi^2 measure (default)\n");
  printf("  y/t    chi^2 measure with Yates' correction\n");
  printf("  i/g    mutual information / G statistic\n");
  printf("  f      Fisher's exact test (table probability)\n");
  printf("  h      Fisher's exact test (chi^2 measure)\n");
  printf("  m      Fisher's exact test (mutual information)\n");
  printf("  s      Fisher's exact test (support)\n");
  printf("\n");
  printf("information output format characters (option -v#)\n");
  printf("  %%%%    a percent sign\n");
  printf("  %%i    number of items (item set size)\n");
  printf("  %%a    absolute item set support\n");
  printf("  %%s    relative item set support as a fraction\n");
  printf("  %%S    relative item set support as a percentage\n");
  printf("  %%p    p-value of item set test as a fraction\n");
  printf("  %%P    p-value of item set test as a percentage\n");
  printf("  %%Q    total transaction weight (database size)\n");
  printf("All format characters can be preceded by the number\n");
  printf("of significant digits to be printed (at most 32 digits),\n");
  printf("even though this value is ignored for integer numbers.\n");
  #endif                        /* print help information */
  exit(0);                      /* abort the program */
}  /* help() */

/*--------------------------------------------------------------------*/
#else

static void help (void)
{                               /* --- print add. option information */
  #ifndef QUIET
  fprintf(stderr, "\n");        /* terminate startup message */
  printf("additional evaluation measures (option -e#)\n");
  printf("frequent item sets:\n");
  printf("  x   no measure (default)\n");
  printf("  b   binary logarithm of support quotient            (+)\n");
  printf("association rules:\n");
  printf("  x   no measure (default)\n");
  printf("  o   rule support (original def.: body & head)       (+)\n");
  printf("  c   rule confidence                                 (+)\n");
  printf("  d   absolute confidence difference to prior         (+)\n");
  printf("  l   lift value (confidence divided by prior)        (+)\n");
  printf("  a   absolute difference of lift value to 1          (+)\n");
  printf("  q   difference of lift quotient to 1                (+)\n");
  printf("  v   conviction (inverse lift for negated head)      (+)\n");
  printf("  e   absolute difference of conviction to 1          (+)\n");
  printf("  r   difference of conviction quotient to 1          (+)\n");
  printf("  k   conditional probability ratio                   (+)\n");
  printf("  j   importance (binary log. of cond. prob. ratio)   (+)\n");
  printf("  z   certainty factor (relative confidence change)   (+)\n");
  printf("  n   normalized chi^2 measure                        (+)\n");
  printf("  p   p-value from (unnormalized) chi^2 measure       (-)\n");
  printf("  y   normalized chi^2 measure with Yates' correction (+)\n");
  printf("  t   p-value from Yates-corrected chi^2 measure      (-)\n");
  printf("  i   information difference to prior                 (+)\n");
  printf("  g   p-value from G statistic/information difference (-)\n");
  printf("  f   Fisher's exact test (table probability)         (-)\n");
  printf("  h   Fisher's exact test (chi^2 measure)             (-)\n");
  printf("  m   Fisher's exact test (information gain)          (-)\n");
  printf("  s   Fisher's exact test (support)                   (-)\n");
  printf("All measures for association rules are also applicable\n");
  printf("to item sets and are then aggregated over all possible\n");
  printf("association rules with a single item in the consequent.\n");
  printf("The aggregation mode can be set with the option -a#.\n");
  printf("Measures marked with (+) "
         "must meet or exceed the threshold,\n");
  printf("measures marked with (-) "
         "must not exceed the threshold\n");
  printf("in order for the rule or item set to be reported.\n");
  printf("\n");
  printf("evaluation measure aggregation modes (option -a#)\n");
  printf("  x   no aggregation (use first value)\n");
  printf("  m   minimum of individual measure values\n");
  printf("  n   maximum of individual measure values\n");
  printf("  a   average of individual measure values\n");
  printf("\n");
  printf("information output format characters (option -v#)\n");
  printf("  %%%%  a percent sign\n");
  printf("  %%i  number of items (item set size)\n");
  printf("  %%a  absolute item set  support\n");
  printf("  %%s  relative item set  support as a fraction\n");
  printf("  %%S  relative item set  support as a percentage\n");
  printf("  %%b  absolute body set  support\n");
  printf("  %%x  relative body set  support as a fraction\n");
  printf("  %%X  relative body set  support as a percentage\n");
  printf("  %%h  absolute head item support\n");
  printf("  %%y  relative head item support as a fraction\n");
  printf("  %%Y  relative head item support as a percentage\n");
  printf("  %%c  rule confidence as a fraction\n");
  printf("  %%C  rule confidence as a percentage\n");
  printf("  %%l  lift value of a rule (confidence/prior)\n");
  printf("  %%L  lift value of a rule as a percentage\n");
  printf("  %%e  additional evaluation measure\n");
  printf("  %%E  additional evaluation measure as a percentage\n");
  printf("  %%q  total transaction weight (database size)\n");
  printf("  %%Q  total transaction weight (database size)\n");
  printf("All format characters can be preceded by the number\n");
  printf("of significant digits to be printed (at most 32 digits),\n");
  printf("even though this value is ignored for integer numbers.\n");
  #endif                        /* print help information */
  exit(0);                      /* abort the program */
}  /* help() */

#endif
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
  if (apriori) apriori_delete(apriori, 0); \
  if (twrite)  twr_delete(twrite, 1);      \
  if (report)  isr_delete(report, 0);      \
  if (tabag)   tbg_delete(tabag,  0);      \
  if (tread)   trd_delete(tread,  1);      \
  if (ibase)   ib_delete (ibase);          \
  if (border)  free(border);
#endif

GENERROR(error, exit)           /* generic error reporting function */

/*--------------------------------------------------------------------*/
#ifdef APRIACC                  /* Accretion-style Apriori */

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
  CCHAR   *dflt    = " (%a,%4P)";    /* default format for check */
  CCHAR   *info    = dflt;      /* format for information output */
  ITEM    zmin     = 2;         /* minimum size of an item set */
  ITEM    zmax     = ITEM_MAX;  /* maximum size of an item set */
  double  smin     = 10;        /* minimum support of an item set */
  int     stat     = 'p';       /* test statistic to use */
  int     sflgs    = 0;         /* test statistic flags */
  double  siglvl   = 1;         /* significance level (in percent) */
  int     sort     = 2;         /* item sorting/processing order */
  int     algo     = APR_BASIC; /* variant of Apriori algorithm */
  int     mode     = APR_DEFAULT;  /* search mode (e.g. pruning) */
  ITEM    prune    = 0;         /* (min. size for) evaluation pruning */
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
    printf("-s#      minimum support of an item set           "
                    "(default: %g)\n", smin);
    printf("         (positive: percentage, "
                     "negative: absolute number)\n");
    printf("-m#      minimum number of items per item set     "
                    "(default: %"ITEM_FMT")\n", zmin);
    printf("-n#      maximum number of items per item set     "
                    "(default: no limit)\n");
    printf("-e#      test statistic for item set evaluation   "
                    "(default: '%c')\n", stat);
    printf("-d#      significance level (maximum p-value)     "
                    "(default: %g%%)\n", siglvl);
    printf("-i       invalidate stat. below expected support  "
                    "(default: evaluate all)\n");
    printf("-p#      minimum set size for subset filtering    "
                    "(default: 0)\n");
    printf("         (0: backward filtering       "
                        "(no subset checks),\n");
    printf("         <0: weak   forward filtering "
                        "(one subset  must qualify),\n");
    printf("         >0: strong forward filtering "
                        "(all subsets must qualify))\n");
    printf("-q#      sort items w.r.t. their frequency        "
                    "(default: %d)\n", sort);
    printf("         (1: ascending, -1: descending, 0: do not sort,\n"
           "          2: ascending, -2: descending w.r.t. "
                    "transaction size sum)\n");
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
    printf("outfile  file to write found item sets to         "
                    "[optional]\n");
    return 0;                   /* print a usage message */
  }                             /* and abort the program */
  #endif  /* #ifndef QUIET */
  /* free option characters: acijloptuxy [A-Z]\[CFPRZ] */

  /* --- evaluate arguments --- */
  for (i = 1; i < argc; i++) {  /* traverse arguments */
    s = argv[i];                /* get option argument */
    if (optarg) { *optarg = s; optarg = NULL; continue; }
    if ((*s == '-') && *++s) {  /* -- if argument is an option */
      while (*s) {              /* traverse options */
        switch (*s++) {         /* evaluate switches */
          case '!': help();                          break;
          case 's': smin   =       strtod(s, &s);    break;
          case 'm': zmin   = (ITEM)strtol(s, &s, 0); break;
          case 'n': zmax   = (ITEM)strtol(s, &s, 0); break;
          case 'e': stat   = (*s) ? *s++ : 'x';      break;
          case 'd': siglvl =       strtod(s, &s);    break;
          case 'i': sflgs  = APR_INVBXS;             break;
          case 'p': prune  = (ITEM)strtol(s, &s, 0); break;
          case 'q': sort   = (int) strtol(s, &s, 0); break;
          case 'F': bdrcnt = getbdr(s, &s, &border); break;
          case 'R': optarg = &fn_sel;                break;
          case 'P': optarg = &fn_psp;                break;
          case 'Z': stats  = 1;                      break;
          case 'g': scan   = 1;                      break;
          #ifdef USE_ZLIB       /* if optional output compression */
          case 'z': mode  |= APR_ZLIB;               break;
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
  if (optarg)       error(E_OPTARG);     /* check (option) arguments */
  if (k      < 1)   error(E_ARGCNT);     /* and number of arguments */
  if (zmin   < 0)   error(E_SIZE, zmin); /* check the size limits */
  if (zmax   < 0)   error(E_SIZE, zmax); /* and the minimum support */
  if (smin   > 100) error(E_SUPPORT, smin);
  if (siglvl > 100) error(E_SIGLVL, siglvl);
  if (bdrcnt < 0)   error(E_NOMEM);
  if ((!fn_inp || !*fn_inp) && (fn_sel && !*fn_sel))
    error(E_STDIN);             /* stdin must not be used twice */
  switch (stat) {               /* check and translate evaluation */
    case 'x': stat = RE_NONE;            break;
    case 'c': stat = RE_CHI2PVAL;        break;
    case 'p': stat = RE_CHI2PVAL;        break;
    case 'n': stat = RE_CHI2PVAL;        break;
    case 'y': stat = RE_YATESPVAL;       break;
    case 't': stat = RE_YATESPVAL;       break;
    case 'i': stat = RE_INFOPVAL;        break;
    case 'g': stat = RE_INFOPVAL;        break;
    case 'f': stat = RE_FETPROB;         break;
    case 'h': stat = RE_FETCHI2;         break;
    case 'm': stat = RE_FETINFO;         break;
    case 's': stat = RE_FETSUPP;         break;
    default : error(E_STAT, (char)stat); break;
  }                             /* (get target type code) */
  stat |= sflgs;                /* add test statistic flags */
  if ((info == dflt) && (smin >= 0))
    info = " (%3S,%4P)";        /* adapt the default info. format */
  mode |= APR_VERBOSE|APR_NOCLEAN;
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
    MSG(stderr, "[%"ITEM_FMT" item(s)]", ib_cnt(ibase));
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
  MSG(stderr, "\n");            /* compute absolute support value */

  /* --- find frequent item sets --- */
  apriori = apriori_create(ISR_MAXIMAL, smin, 100.0, 100.0,
                           zmin, zmax, stat, APR_MAX, siglvl,
                           algo, mode);
  if (!apriori) error(E_NOMEM); /* create an Apriori miner */
  k = apriori_data(apriori, tabag, 0, sort);
  if (k) error(k);              /* prepare data for Apriori */
  report = isr_create(ibase);   /* create an item set reporter */
  if (!report) error(E_NOMEM);  /* and configure it */
  k = apriori_report(apriori, report);
  if (k) error(k);              /* prepare reporter for Apriori */
  if (setbdr(report, w, zmin, &border, bdrcnt) != 0)
    error(E_NOMEM);             /* set the support border (if any) */
  if (fn_psp && (isr_addpsp(report, NULL) < 0))
    error(E_NOMEM);             /* set a pattern spectrum if req. */
  if (isr_setfmt(report, scan, hdr, sep, NULL, info) != 0)
    error(E_NOMEM);             /* set the output format strings */
  k = isr_open(report, NULL, fn_out);
  if (k) error(k, isr_name(report));
  if (isr_setup(report) < 0)    /* open the output file and */
    error(E_NOMEM);             /* set up the item set reporter */
  k = apriori_mine(apriori, prune, 0.0, 0);
  if (k) error(k);              /* find frequent item sets */
  if (stats)                    /* print item set statistics */
    isr_prstats(report, stdout, 0);
  if (isr_close(report) != 0)   /* close the output file */
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

/*--------------------------------------------------------------------*/
#else                           /* standard Apriori algorithm */

int main (int argc, char *argv[])
{                               /* --- main function */
  int     i, k = 0;             /* loop variables, counters */
  char    *s;                   /* to traverse the options */
  CCHAR   **optarg = NULL;      /* option argument */
  CCHAR   *fn_inp  = NULL;      /* name of the input  file */
  CCHAR   *fn_out  = NULL;      /* name of the output file */
  CCHAR   *fn_sel  = NULL;      /* name of item selection file */
  CCHAR   *fn_psp  = NULL;      /* name of pattern spectrum file */
  CCHAR   *recseps = NULL;      /* record  separators */
  CCHAR   *fldseps = NULL;      /* field   separators */
  CCHAR   *blanks  = NULL;      /* blank   characters */
  CCHAR   *comment = NULL;      /* comment characters */
  CCHAR   *hdr     = "";        /* record header  for output */
  CCHAR   *sep     = " ";       /* item separator for output */
  CCHAR   *imp     = " <- ";    /* implication sign for ass. rules */
  CCHAR   *dflt    = " (%S)";   /* default format for check */
  CCHAR   *info    = dflt;      /* format for information output */
  int     target   = 's';       /* target type (e.g. closed/maximal) */
  ITEM    zmin     = 1;         /* minimum rule/item set size */
  ITEM    zmax     = ITEM_MAX;  /* maximum rule/item set size */
  double  smin     = 10;        /* minimum support of an item set */
  double  smax     = 100;       /* maximum support of an item set */
  double  conf     = 80;        /* minimum confidence (in percent) */
  int     eval     = 'x';       /* additional evaluation measure */
  int     eflgs    = 0;         /* evaluation measure flags */
  int     agg      = 'x';       /* aggregation mode for eval. measure */
  double  thresh   = 10;        /* threshold for evaluation measure */
  int     sort     = 2;         /* flag for item sorting and recoding */
  int     algo     = APR_BASIC; /* variant of apriori algorithm */
  int     mode     = APR_DEFAULT|APR_PREFMT;     /* search mode */
  ITEM    prune    = ITEM_MIN;  /* (min. size for) evaluation pruning */
  double  filter   = 0.01;      /* item usage filtering parameter */
  int     order    = 0;         /* size order item set/rule output */
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
  else {                        /* if no arguments are given */
    printf("usage: %s [options] infile [outfile]\n", argv[0]);
    printf("%s\n", DESCRIPTION);
    printf("%s\n", VERSION);
    printf("-t#      target type                              "
                    "(default: %c)\n", target);
    printf("         (s: frequent, c: closed, m: maximal item sets,\n");
    printf("          g: generators, r: association rules)\n");
    printf("-m#      minimum number of items per set/rule     "
                    "(default: %"ITEM_FMT")\n", zmin);
    printf("-n#      maximum number of items per set/rule     "
                    "(default: no limit)\n");
    printf("-s#      minimum support of an item set/rule      "
                    "(default: %g%%)\n", smin);
    printf("-S#      maximum support of an item set/rule      "
                    "(default: %g%%)\n", smax);
    printf("         (positive: percentage, "
                     "negative: absolute number)\n");
    printf("-o       use original rule support definition     "
                    "(body & head)\n");
    printf("-c#      minimum confidence of an assoc. rule     "
                    "(default: %g%%)\n", conf);
    printf("-e#      additional evaluation measure            "
                    "(default: none)\n");
    printf("-a#      aggregation mode for evaluation measure  "
                    "(default: none)\n");
    printf("-d#      threshold for add. evaluation measure    "
                    "(default: %g%%)\n", thresh);
    printf("-i       invalidate eval. below expected support  "
                    "(default: evaluate all)\n");
    printf("-p#      (min. size for) pruning with evaluation  "
                    "(default: no pruning)\n");
    printf("         (< 0: weak forward, > 0 strong forward, "
                     "= 0: backward pruning)\n");
    printf("-q#      sort items w.r.t. their frequency        "
                    "(default: %d)\n", sort);
    printf("         (1: ascending, -1: descending, 0: do not sort,\n"
           "          2: ascending, -2: descending w.r.t. "
                    "transaction size sum)\n");
    printf("-u#      filter unused items from transactions    "
                    "(default: %g)\n", filter);
    printf("         (0: do not filter items w.r.t. usage in sets,\n"
           "         <0: fraction of removed items for filtering,\n"
           "         >0: take execution times ratio into account)\n");
    printf("-x       do not prune with perfect extensions     "
                    "(default: prune)\n");
    printf("-y       a-posteriori pruning of infrequent item sets\n");
    printf("-T       do not organize transactions as a prefix tree\n");
    printf("-F#:#..  support border for filtering item sets   "
                    "(default: none)\n");
    printf("         (list of minimum support values, "
                    "one per item set size,\n");
    printf("         starting at the minimum size, "
                    "as given with option -m#)\n");
    printf("-R#      read item selection/appearance indicators\n");
    printf("-P#      write a pattern spectrum to a file\n");
    printf("-Z       print item set statistics "
                    "(number of item sets per size)\n");
    printf("-N       do not pre-format some integer numbers   "
                    "(default: do)\n");
    printf("-g       write item names in scanable form "
                    "(quote certain characters)\n");
    #ifdef USE_ZLIB             /* if optional output compression */
    printf("-z       compress output with zlib (deflate)      "
                    "(default: plain text)\n");
    #endif                      /* print compression option */
    printf("-h#      record header  for output                "
                    "(default: \"%s\")\n", hdr);
    printf("-k#      item separator for output                "
                    "(default: \"%s\")\n", sep);
    printf("-I#      implication sign for association rules   "
                    "(default: \"%s\")\n", imp);
    printf("-v#      output format for set/rule information   "
                    "(default: \"%s\")\n", info);
    printf("-j#      sort item sets in output by their size   "
                    "(default: no sorting)\n");
    printf("         (< 0: descending, > 0: ascending order)\n");
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
    printf("outfile  file to write item sets/assoc. rules to  "
                    "[optional]\n");
    return 0;                   /* print a usage message */
  }                             /* and abort the program */
  #endif  /* #ifndef QUIET */
  /* free option characters: l [A-Z]\[CFINPRSTZ] */

  /* --- evaluate arguments --- */
  for (i = 1; i < argc; i++) {  /* traverse the arguments */
    s = argv[i];                /* get an option argument */
    if (optarg) { *optarg = s; optarg = NULL; continue; }
    if ((*s == '-') && *++s) {  /* -- if argument is an option */
      while (*s) {              /* traverse the options */
        switch (*s++) {         /* evaluate the options */
          case '!': help();                          break;
          case 't': target = (*s) ? *s++ : 's';      break;
          case 'm': zmin   = (ITEM)strtol(s, &s, 0); break;
          case 'n': zmax   = (ITEM)strtol(s, &s, 0); break;
          case 's': smin   =       strtod(s, &s);    break;
          case 'S': smax   =       strtod(s, &s);    break;
          case 'o': mode  |= APR_ORIGSUPP;           break;
          case 'c': conf   =       strtod(s, &s);    break;
          case 'e': eval   = (*s) ? *s++ : 0;        break;
          case 'a': agg    = (*s) ? *s++ : 0;        break;
          case 'd': thresh =       strtod(s, &s);    break;
          case 'i': eflgs |= APR_INVBXS;             break;
          case 'p': prune  = (ITEM)strtol(s, &s, 0); break;
          case 'q': sort   = (int) strtol(s, &s, 0); break;
          case 'u': filter =       strtod(s, &s);    break;
          case 'x': mode  &= ~APR_PERFECT;           break;
          case 'y': mode  |=  APR_POST;              break;
          case 'T': mode  &= ~APR_TATREE;            break;
          case 'F': bdrcnt = getbdr(s, &s, &border); break;
          case 'R': optarg = &fn_sel;                break;
          case 'P': optarg = &fn_psp;                break;
          case 'Z': stats  = 1;                      break;
          case 'N': mode  &= ~APR_PREFMT;            break;
          case 'g': scan   = 1;                      break;
          #ifdef USE_ZLIB       /* if optional output compression */
          case 'z': mode  |= APR_ZLIB;               break;
          #endif                /* set the compression flag */
          case 'h': optarg = &hdr;                   break;
          case 'k': optarg = &sep;                   break;
          case 'I': optarg = &imp;                   break;
          case 'v': optarg = &info;                  break;
          case 'j': order  = (int) strtol(s, &s, 0); break;
          case 'w': mtar  |= TA_WEIGHT;              break;
          case 'r': optarg = &recseps;               break;
          case 'f': optarg = &fldseps;               break;
          case 'b': optarg = &blanks;                break;
          case 'C': optarg = &comment;               break;
          default : error(E_OPTION, *--s);           break;
        }                       /* set the option variables */
        if (optarg && *s) { *optarg = s; optarg = NULL; break; }
      } }                       /* get an option argument */
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
  if ((conf  < 0) || (conf > 100))
    error(E_CONF, conf);        /* check the minimum confidence */
  if ((!fn_inp || !*fn_inp) && (fn_sel && !*fn_sel))
    error(E_STDIN);             /* stdin must not be used twice */
  switch (target) {             /* check and translate target type */
    case 's': target = ISR_ALL;              break;
    case 'f': target = ISR_FREQUENT;         break;
    case 'c': target = ISR_CLOSED;           break;
    case 'm': target = ISR_MAXIMAL;          break;
    case 'g': target = ISR_GENERAS;          break;
    case 'r': target = ISR_RULES;            break;
    default : error(E_TARGET, (char)target); break;
  }                             /* (get target type code) */
  switch (eval) {               /* check and translate measure */
    case 'x': eval = RE_NONE;                break;
    case 'o': eval = RE_SUPP;                break;
    case 'c': eval = RE_CONF;                break;
    case 'd': eval = RE_CONFDIFF;            break;
    case 'l': eval = RE_LIFT;                break;
    case 'a': eval = RE_LIFTDIFF;            break;
    case 'q': eval = RE_LIFTQUOT;            break;
    case 'v': eval = RE_CVCT;                break;
    case 'e': eval = RE_CVCTDIFF;            break;
    case 'r': eval = RE_CVCTQUOT;            break;
    case 'k': eval = RE_CPROB;               break;
    case 'j': eval = RE_IMPORT;              break;
    case 'z': eval = RE_CERT;                break;
    case 'n': eval = RE_CHI2;                break;
    case 'p': eval = RE_CHI2PVAL;            break;
    case 'y': eval = RE_YATES;               break;
    case 't': eval = RE_YATESPVAL;           break;
    case 'i': eval = RE_INFO;                break;
    case 'g': eval = RE_INFOPVAL;            break;
    case 'f': eval = RE_FETPROB;             break;
    case 'h': eval = RE_FETCHI2;             break;
    case 'm': eval = RE_FETINFO;             break;
    case 's': eval = RE_FETSUPP;             break;
    case 'b': eval = APR_LDRATIO;            break;
    default : error(E_MEASURE, (char)eval);  break;
  }  /* free: u w */            /* (get evaluation measure code) */
  eval |= eflgs;                /* add evaluation measure flags */
  switch (agg) {                /* check and translate agg. mode */
    case 'x': agg  = APR_NONE;               break;
    case 'm': agg  = APR_MIN;                break;
    case 'n': agg  = APR_MAX;                break;
    case 'a': agg  = APR_AVG;                break;
    default : error(E_AGGMODE, (char)agg);   break;
  }                             /* (get aggregation mode code) */
  if ((filter <= -1) || (filter >= 1))
    filter = 0;                 /* check and adapt the filter option */
  if (target & ISR_RULES)       /* if to find association rules, */
    fn_psp = NULL;              /* no pattern spectrum possible */
  if (info == dflt) {           /* if default info. format is used, */
    if (target != ISR_RULES)    /* set default according to target */
         info = (smin < 0) ? " (%a)"     : " (%S)";
    else info = (smin < 0) ? " (%b, %C)" : " (%X, %C)";
  }                             /* select absolute/relative support */
  mode |= APR_VERBOSE|APR_NOCLEAN;
  MSG(stderr, "\n");            /* terminate the startup message */

  /* --- read item selection/appearance indicators --- */
  ibase = ib_create(0, 0);      /* create an item base */
  if (!ibase) error(E_NOMEM);   /* to manage the items */
  tread = trd_create();         /* create a transaction reader */
  if (!tread) error(E_NOMEM);   /* and configure the characters */
  trd_allchs(tread, recseps, fldseps, blanks, "", comment);
  if (fn_sel) {                 /* if an item selection is given */
    CLOCK(t);                   /* start timer, open input file */
    if (trd_open(tread, NULL, fn_sel) != 0)
      error(E_FOPEN, trd_name(tread));
    MSG(stderr, "reading %s ... ", trd_name(tread));
    m = (target == ISR_RULES)   /* depending on the target type */
      ? ib_readapp(ibase,tread) /* read the item appearances */
      : ib_readsel(ibase,tread);/* or a simple item selection */
    if (m < 0) error((int)-m, ib_errmsg(ibase, NULL, 0));
    trd_close(tread);           /* close the input file */
    MSG(stderr, "[%"ITEM_FMT" item(s)]", ib_cnt(ibase));
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

  /* --- find frequent item sets/association rules --- */
  apriori = apriori_create(target, smin, smax, conf, zmin, zmax,
                           eval, agg, thresh, algo, mode);
  if (!apriori) error(E_NOMEM); /* create an Apriori miner */
  k = apriori_data(apriori, tabag, 0, sort);
  if (k) error(k);              /* prepare data for Apriori */
  report = isr_create(ibase);   /* create an item set reporter */
  if (!report) error(E_NOMEM);  /* and configure it */
  k = apriori_report(apriori, report);
  if (k) error(k);              /* prepare reporter for Apriori */
  if (setbdr(report, w, zmin, &border, bdrcnt) != 0)
    error(E_NOMEM);             /* set the support border (if any) */
  if (fn_psp && (isr_addpsp(report, NULL) < 0))
    error(E_NOMEM);             /* set a pattern spectrum if req. */
  if (isr_setfmt(report, scan, hdr, sep, imp, info) != 0)
    error(E_NOMEM);             /* set the output format strings */
  k = isr_open(report, NULL, fn_out);
  if (k) error(k, isr_name(report));
  if (isr_setup(report) < 0)    /* open the output file and */
    error(E_NOMEM);             /* set up the item set reporter */
  k = apriori_mine(apriori, prune, filter, order);
  if (k) error(k);              /* find frequent item sets */
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

#endif  /* #ifdef APRIACC ... #else ... */
#endif  /* #ifdef APR_MAIN ... */
