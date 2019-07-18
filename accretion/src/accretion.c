/*----------------------------------------------------------------------
  File    : accretion.c
  Contents: accretion algorithm for identifying neural assemblies
  Author  : Christian Borgelt
  History : 2011.06.22 file created from file eclat.c
            2011.06.23 Fisher's exact test added (in various forms)
            2011.07.08 adapted to modified function tbg_recode()
            2011.07.22 adapted to new module ruleval (rule evaluation)
            2011.08.28 output of item set counters per size added
            2012.02.06 bug in selection of items to keep fixed
            2012.11.06 optional output of all and closed item sets added
            2013.01.31 option to invalidate statistic below expectation
            2013.03.07 direction parameter added to sorting functions
            2013.03.29 adapted to type changes in module tract
            2013.04.17 bug concerning support handling fixed
            2013.10.15 checks of return code of isr_report() added
            2013.10.18 optional pattern spectrum collection added
            2013.11.12 item selection file changed to option -R#
            2014.05.12 option -F# added (support border for filtering)
            2014.08.24 adapted to modified item set reporter interface
            2014.08.28 functions acc_data() and acc_report() added
            2014.10.24 changed from LGPL license to MIT license
            2016.11.15 accretion miner object and interface introduced
            2017.05.30 optional output compression with zlib added
------------------------------------------------------------------------
  Reference for the Accretion algorithm:
    G.L. Gerstein, D.H. Perkel and K.N. Subramanian.
    Identification of Functionally Related Neural Assemblies.
    Brain Research 140(1):43-62.
    Elsevier, Amsterdam, Netherlands 1978
  Reference for the Eclat algorithm:
    M.J. Zaki, S. Parthasarathy, M. Ogihara, and W. Li.
    New Algorithms for Fast Discovery of Association Rules.
    Proc. 3rd Int. Conf. on Knowledge Discovery and Data Mining
    (KDD'97, Newport Beach, CA), 283-296.
    AAAI Press, Menlo Park, CA, USA 1997
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
#ifdef ACC_MAIN
#ifndef PSP_REPORT
#define PSP_REPORT
#endif
#ifndef TA_READ
#define TA_READ
#endif
#endif
#ifdef ACC_ABORT
#include "sigint.h"
#endif
#include "accretion.h"
#ifdef ACC_MAIN
#include "error.h"
#endif
#ifdef STORAGE
#include "storage.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define PRGNAME     "accretion"
#define DESCRIPTION "accretion algorithm " \
                    "for identifying neural assemblies"
#define VERSION     "version 2.18 (2017.05.30)        " \
                    "(c) 2011-2017   Christian Borgelt"

/* --- error codes --- */
/* error codes   0 to  -4 defined in tract.h */
#define E_STDIN      (-5)       /* double assignment of stdin */
#define E_OPTION     (-6)       /* unknown option */
#define E_OPTARG     (-7)       /* missing option argument */
#define E_ARGCNT     (-8)       /* too few/many arguments */
#define E_TARGET     (-9)       /* invalid target type */
#define E_SIZE      (-10)       /* invalid item set size */
#define E_SUPPORT   (-11)       /* invalid minimum item set support */
#define E_STAT      (-12)       /* invalid test statistic */
#define E_SIGLVL    (-13)       /* invalid significance level */
/* error codes -15 to -25 defined in tract.h */

#define DIFFSIZE(p,q) ((size_t)((int*)(p)-(int*)(q)) *sizeof(int))

#ifndef QUIET                   /* if not quiet version, */
#define MSG         fprintf     /* print messages */
#define XMSG        if (accret->mode & ACC_VERBOSE) fprintf
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
typedef struct {                /* --- trans. identifier list --- */
  ITEM      item;               /* item identifier (last item in set) */
  SUPP      supp;               /* support of the item (or item set) */
  double    pval;               /* p-value of statistical test */
  TID       tids[1];            /* transaction identifiers */
} TIDLIST;                      /* (transaction identifier list) */

struct _accret {                /* --- accretion miner --- */
  int       target;             /* target type (e.g. closed/maximal) */
  double    smin;               /* minimum support of an item set */
  double    smax;               /* maximum support of an item set */
  SUPP      supp;               /* minimum support of an item set */
  ITEM      zmin;               /* minimum size of a rule/item set */
  ITEM      zmax;               /* maximum size of a rule/item set */
  int       stat;               /* evaluation statistic */
  int       invbxs;             /* invalidate stat. below expectation */
  RULEVALFN *statfn;            /* function for test statistic */
  double    siglvl;             /* significance level */
  int       mode;               /* operation/search mode */
  ITEM      maxext;             /* maximum number of extensions */
  TABAG     *tabag;             /* original transaction bag */
  ISREPORT  *report;            /* item set reporter */
  SUPP      ttw;                /* total transaction weight */
  TIDLIST   **lists;            /* transaction identifier lists */
  SUPP      *muls;              /* multiplicity of transactions */
  SUPP      *marks;             /* flags for tid occurrences */
} RECDATA;                      /* (recursion data) */

/*----------------------------------------------------------------------
  Constants
----------------------------------------------------------------------*/
#if !defined QUIET && defined ACC_MAIN
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
  /* E_STAT    -12 */  "invalid test statistic '%c'",
  /* E_SIGLVL  -13 */  "invalid significance level/p-value %g",
  /*           -14 */  NULL,
  /* E_NOITEMS -15 */  "no (frequent) items found",
  /*           -16 */  "unknown error"
};
#endif

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
#ifdef ACC_MAIN
#ifndef QUIET
static CCHAR    *prgname;       /* program name for error messages */
#endif
static TABREAD  *tread  = NULL; /* table/transaction reader */
static ITEMBASE *ibase  = NULL; /* item base */
static TABAG    *tabag  = NULL; /* transaction bag/multiset */
static ISREPORT *report = NULL; /* item set reporter */
static TABWRITE *twrite = NULL; /* table writer for pattern spectrum */
static double   *border = NULL; /* support border for filtering */
static ACCRET   *accret = NULL; /* accretion miner object */
#endif

/*----------------------------------------------------------------------
  Auxiliary Functions (for debugging)
----------------------------------------------------------------------*/
#if !defined NDEBUG && defined ACC_MAIN

static void indent (int k)
{ while (--k >= 0) printf("   "); }

/*--------------------------------------------------------------------*/

static void show (const char *text, TIDLIST **lists, int k, int ind)
{                               /* --- show a cond. trans. database */
  ITEM i, j;                    /* loop variable */
  TID  *s;                      /* to traverse the tids */

  if (text && *text) {          /* print the given text */
    indent(ind); printf("%s\n", text); }
  for (j = 0; j < k; j++) {     /* traverse the items/tid lists */
    indent(ind);                /* indent the output line */
    i = lists[j]->item;         /* print the item name and id */
    printf("%4s[%2"ITEM_FMT"]:", ib_name(ibase, i), i);
    for (s = lists[j]->tids; *s >= 0; s++)
      printf(" %"TID_FMT, *s);  /* print the item and the tids */
    printf(" (%"SUPP_FMT")\n", lists[i]->supp);
  }                             /* print the item support */
}  /* show() */

#endif
/*----------------------------------------------------------------------
  Accretion (with an Eclat-style scheme)
----------------------------------------------------------------------*/

static TID isect (TIDLIST *dst, TIDLIST *src1, TIDLIST *src2,SUPP *muls)
{                               /* --- intersect two tid lists */
  TID *s1, *s2, *d;             /* to traverse sources and dest. */

  assert(dst && src1 && src2    /* check the function arguments */
  &&    (src1->tids[0] >= 0) && (src2->tids[0] >= 0) && muls);
  dst->item = src1->item;       /* copy the first item and */
  dst->supp = 0;                /* initialize the support */
  if (src1->supp > src2->supp) { s2 = src1->tids; s1 = src2->tids; }
  else                         { s1 = src1->tids; s2 = src2->tids; }
  d = dst->tids;                /* get sources and destination */
  while (1) {                   /* tid list intersection loop */
    if      (*s1 < *s2) s2++;   /* if one transaction id is larger, */
    else if (*s1 > *s2) s1++;   /* simply skip this transaction id */
    else if (*s1 < 0) break;    /* check for the sentinel */
    else { dst->supp += muls[*d++ = *s1++]; s2++; }
  }                             /* copy equal elements to destination */
  *d++ = (TID)-1;               /* store a sentinel at the list end */
  return (TID)(d -dst->tids);   /* return the size of the new list */
}  /* isect() */

/*--------------------------------------------------------------------*/

static TID filter (TIDLIST *dst, TIDLIST *src, SUPP *muls)
{                               /* --- filter a tid list */
  SUPP m;                       /* multiplicity of transaction */
  TID  *s, *d;                  /* to traverse source and dest. */

  assert(dst && src && muls);   /* check function arguments */
  dst->item = src->item;        /* copy the first item and */
  dst->supp = 0;                /* initialize the support */
  for (d = dst->tids, s = src->tids; *s >= 0; s++)
    if ((m = muls[*s]) > 0) {   /* collect the marked trans. ids and */
      dst->supp += m; *d++ = *s; }    /* sum the transaction weights */
  *d++ = (TID)-1;               /* store a sentinel at the list end */
  return (TID)(d -dst->tids);   /* return the size of the new list */
}  /* filter() */

/*--------------------------------------------------------------------*/

static int cmp (const void *a, const void *b, void *data)
{                               /* --- compare tid list p-values */
  if (((TIDLIST*)a)->pval < ((TIDLIST*)b)->pval) return -1;
  if (((TIDLIST*)a)->pval > ((TIDLIST*)b)->pval) return +1;
  if (((TIDLIST*)a)->supp > ((TIDLIST*)b)->supp) return -1;
  if (((TIDLIST*)a)->supp < ((TIDLIST*)b)->supp) return +1;
  return 0;                     /* return sign of p-value difference */
}  /* cmp() */

/*--------------------------------------------------------------------*/

static SUPP recurse (ACCRET *accret, TIDLIST **lists, ITEM k, size_t x)
{                               /* --- eclat recursion with i.section */
  int     r;                    /* error status */
  ITEM    i, j, m, z;           /* loop variables */
  SUPP    s, smax;              /* (maximum) support of an item set */
  TIDLIST *l, *d;               /* to traverse the tid lists */
  TIDLIST **proj = NULL;        /* tid lists of projected database */
  TID     *p, *q;               /* to organize/traverse the tid lists */

  assert(accret && lists && (k > 0)); /* check the function arguments */
  #ifdef ACC_ABORT              /* if to check for interrupt */
  if (sig_aborted()) return -1; /* if execution was aborted, */
  #endif                        /* abort the recursion */
  if ((k > 1)                   /* if there is more than one item and */
  &&  isr_xable(accret->report, 2)) {    /* another item can be added */
    proj = (TIDLIST**)malloc((size_t)k *sizeof(TIDLIST*) +x);
    if (!proj) return -1;       /* allocate list and element arrays */
  }                             /* (memory for conditional databases) */
  smax = 0;                     /* clear the maximum item set support */
  ptr_qsort(lists, (size_t)k, +1, cmp, NULL);
  z = isr_cnt(accret->report);  /* sort by p-value of items/tid-lists */
  z = ((z <= 0) || (k < accret->maxext)) ? k : accret->maxext;
  for (i = r = 0; i < z; i++) { /* get the max. number of extensions */
    l = lists[i];               /* and traverse the items/tid lists */
    if (l->pval > accret->siglvl)  /* skip all extension items that */
      break;                    /* are not signficantly correlated */
    r = isr_add(accret->report, l->item, l->supp);
    if (r < 0) break;           /* add current item to the reporter */
    s = 0;                      /* default: no report in recursion */
    if (proj) {                 /* if another item can be added */
      proj[m = 0] = d = (TIDLIST*)(p = (TID*)(proj +k+1));
      if (k <= 2) {             /* if there is only one other item */
        /* Benchmark tests showed that this version is faster only */
        /* if there is only one other tid list to intersect with.  */
        for (j = 0; j < k; j++){/* intersect with other tid lists */
          if (j == i) continue; /* (skip the current item) */
          x = (size_t)isect(d, lists[j], l, accret->muls);
          if (d->supp < accret->supp)
            continue;           /* skip items that are infrequent */
          s = accret->lists[lists[j]->item]->supp;
          d->pval = accret->statfn(d->supp, l->supp, s, accret->ttw);
          proj[++m] = d = (TIDLIST*)(p = d->tids +x);
        } }                     /* collect tid lists of sign. items */
      else {                    /* if there are many items left */
        for (q = l->tids; *q >= 0; q++) /* mark transaction ids in */
          accret->marks[*q] = accret->muls[*q];    /* current list */
        for (j = 0; j < k; j++){/* intersect with other tid lists */
          if (j == i) continue; /* (skip the current item) */
          x = (size_t)filter(d, lists[j], accret->marks);
          if (d->supp < accret->supp)
            continue;           /* skip items that are infrequent */
          s = accret->lists[lists[j]->item]->supp;
          d->pval = (!accret->invbxs
                 || ((double)d->supp *(double)accret->ttw
                  >  (double)l->supp *(double)s))
                  ? accret->statfn(d->supp, l->supp, s, accret->ttw) :1;
          proj[++m] = d = (TIDLIST*)(p = d->tids +x);
        }                       /* collect tid lists of sign. items */
        for (q = l->tids; *q >= 0; q++)
          accret->marks[*q] = 0;/* unmark transaction ids */
      }                         /* in the current list */
      if (m > 0) {              /* if the projection is not empty */
        s = recurse(accret, proj, m, DIFFSIZE(p, proj[0]));
        if (s < 0) { r = (int)s; break; }
        if (s > smax) smax = s; /* recursively find freq. item sets */
      }                         /* in the created projection and */
    }                           /* update the maximum support */
    if (!(accret->target & (ISR_CLOSED|ISR_MAXIMAL))
    ||  ((accret->target & ISR_MAXIMAL) && (s < accret->supp))
    ||  ((accret->target & ISR_CLOSED)  && (s <      l->supp))) {
      if (l->supp > smax) smax = l->supp;
      r = isr_reportv(accret->report, l->pval);
      if (r < 0) break;         /* if current item set qualifies, */
    }                           /* report the current item set */
    isr_remove(accret->report, 1);
  }                             /* remove current item from reporter */
  if (proj) free(proj);         /* delete the list and element arrays */
  return (r < 0) ? (SUPP)r : smax;
}  /* recurse() */              /* return error status or max. supp */

/*--------------------------------------------------------------------*/

int accret_base (ACCRET *accret)
{                               /* --- search for frequent item sets */
  ITEM       i, k, m;           /* loop variable, number of items */
  TID        n;                 /* number transactions */
  size_t     x;                 /* extent (number of item instances) */
  SUPP       w;                 /* weight/support buffer */
  TRACT      *t;                /* to traverse the transactions */
  TIDLIST    **lists, *l;       /* to traverse the tid lists */
  TID        *tids, *p, **next; /* to traverse transaction ids */
  const ITEM *s;                /* to traverse transaction items */
  const TID  *c;                /* item occurrence counters */

  assert(accret);               /* check the function arguments */
  if (accret->supp > accret->ttw)
    return 0;                   /* check against minimum support */
  k = tbg_itemcnt(accret->tabag);
  if (k <= 0)                   /* get and check the number of items */
    return isr_reportv(accret->report, 1);
  n = tbg_cnt(accret->tabag);   /* get the number of transactions */
  c = tbg_icnts(accret->tabag, 0);
  if (!c) return -1;            /* get the number of containing */
  accret->lists =               /* transactions per item */
  lists = (TIDLIST**)malloc((size_t)(k+k) *sizeof(TIDLIST*)
                           +(size_t) k    *sizeof(TID*)
                           +(size_t)(n+n) *sizeof(SUPP));
  if (!lists) return -1;        /* create initial tid list array */
  next     = (TID**)(lists+k+k);/* and split off next position array, */
  accret->muls  = (SUPP*)(next+k);   /* transaction weight array, and */
  accret->marks = accret->muls +n;   /* transaction flags array */
  memset(accret->marks, 0, (size_t)n *sizeof(SUPP));
  x = tbg_extent(accret->tabag);/* get the number of item occurrences */
  p = tids = (TID*)malloc((size_t)k *sizeof(TIDLIST) +x *sizeof(TID));
  if (!p) { free(lists); return -1; } /* allocate tid list elements */
  for (i = 0; i < k; i++) {     /* traverse the items/tid lists */
    lists[i] = l = (TIDLIST*)p; /* get/create the next tid list */
    l->item  = i;               /* initialize the list item */
    l->supp  = 0;               /* and the support counter */
    l->pval  = 0;               /* clear the p-value (significant) */
    next[i]  = l->tids;         /* note position of next trans. id */
    p = l->tids +c[i] +1;       /* skip space for transaction ids */
  }                             /* and a sentinel at the end */
  while (n > 0) {               /* traverse the transactions */
    t = tbg_tract(accret->tabag, --n); /* get the next transaction */
    accret->muls[n] = w = ta_wgt(t);   /* and store its weight */
    for (s = ta_items(t); *s > TA_END; s++) {
      lists[*s]->supp += w;     /* traverse the transaction's items */
      *next[*s]++ = (TID)n;     /* sum the transaction weight and */
    }                           /* collect the transaction ids */
  }
  for (lists += k, i = m = 0; i < k; i++) {
    l = accret->lists[i];       /* traverse the items and eliminate */
    if (l->supp < accret->supp) continue;       /* infrequent items */
    *next[i]   = (TID)-1;       /* store a sentinel at the list end */
    lists[m++] = l;             /* collect lists for frequent items */
  }                             /* (eliminate infrequent items) */
  w = 0;                        /* init. return code/ext. support */
  if (m > 0)                    /* find freq. items sets recursively */
    w = recurse(accret, lists, m, DIFFSIZE(p, tids));
  if (!(accret->target & (ISR_CLOSED|ISR_MAXIMAL))
  ||  ((accret->target & ISR_MAXIMAL) && (w < accret->supp))
  ||  ((accret->target & ISR_CLOSED)  && (w < accret->ttw))) {
    if (isr_reportv(accret->report, 1) < 0) w = -1; }
  free(tids);                   /* finally report empty set and */
  free(accret->lists);          /* delete the allocated arrays */
  return (w < 0) ? (int)w : 0;  /* return the error status */
}  /* accret_base() */

/*----------------------------------------------------------------------
  Accretion (generic)
----------------------------------------------------------------------*/

ACCRET* accret_create (int target, double smin, double smax,
                       ITEM zmin, ITEM zmax,
                       int stat, double siglvl, int mode)
{                               /* --- create a accretion miner */
  ACCRET *accret;               /* created accretion miner */

  assert((stat & ~ACC_INVBXS) < RE_FNCNT);

  /* --- make parameters consistent --- */
  if      (target & ISR_MAXIMAL) target = ISR_MAXIMAL;
  else if (target & ISR_CLOSED)  target = ISR_CLOSED;
  else                           target = ISR_FREQUENT;

  /* --- create an accretion miner --- */
  accret = (ACCRET*)malloc(sizeof(ACCRET));
  if (!accret) return NULL;      /* create an accretion miner */
  accret->target = target;       /* and store all parameters */
  accret->smin   = smin;
  accret->smax   = smax;
  accret->supp   = 1;
  accret->zmin   = zmin;
  accret->zmax   = zmax;
  accret->stat   = stat & ~ACC_INVBXS;
  accret->invbxs = stat &  ACC_INVBXS;
  accret->statfn = re_function(accret->stat);
  accret->siglvl = (siglvl > 0) ? siglvl/100.0 : 0.01;
  accret->maxext = 2;
  accret->mode   = mode;
  accret->tabag  = NULL;
  accret->report = NULL;
  accret->ttw    = 0;
  accret->lists  = NULL;
  accret->muls   = NULL;
  accret->marks  = NULL;
  return accret;                 /* return created accretion miner */
}  /* accret_create() */

/*--------------------------------------------------------------------*/

void accret_delete (ACCRET *accret, int deldar)
{                               /* --- delete an accretion miner */
  if (deldar) {                 /* if to delete data and reporter */
    if (accret->report) isr_delete(accret->report, 0);
    if (accret->tabag)  tbg_delete(accret->tabag,  1);
  }                             /* delete if existing */
  free(accret);                 /* delete the base structure */
}  /* accret_delete() */

/*--------------------------------------------------------------------*/

int accret_data (ACCRET *accret, TABAG *tabag, int sort)
{                               /* --- prepare data for Accretion */
  ITEM    m;                    /* number of items */
  double  smin;                 /* absolute minimum support */
  #ifndef QUIET                 /* if to print messages */
  TID     n;                    /* number of transactions */
  SUPP    w;                    /* total transaction weight */
  clock_t t;                    /* timer for measurements */
  #endif                        /* (only needed for messages) */

  assert(accret && tabag);      /* check the function arguments */
  accret->tabag = tabag;        /* note the transaction bag */

  /* --- compute data-specific parameters --- */
  accret->ttw = tbg_wgt(tabag); /* compute absolute minimum support */
  smin = ceilsupp((accret->smin < 0) ? -accret->smin
                : (accret->smin/100.0) *(double)accret->ttw
                                       *(1-DBL_EPSILON));
  accret->supp = (SUPP)ceilsupp(smin);

  /* --- sort and recode items --- */
  CLOCK(t);                     /* start timer, print log message */
  XMSG(stderr, "filtering, sorting and recoding items ... ");
  m = tbg_recode(tabag, accret->supp, -1, -1, -sort);
  if (m <  0) return E_NOMEM;   /* recode items and transactions */
  if (m <= 0) return E_NOITEMS; /* and check the number of items */
  XMSG(stderr, "[%"ITEM_FMT" item(s)]", m);
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));

  /* --- sort and reduce transactions --- */
  CLOCK(t);                     /* start timer, print log message */
  XMSG(stderr, "sorting and reducing transactions ... ");
  tbg_itsort(tabag, -1, 0);     /* sort items in transactions and */
  tbg_sort  (tabag, -1, 0);     /* sort the trans. lexicographically */
  tbg_reduce(tabag, 0);         /* reduce transactions to unique ones */
  #ifndef QUIET                 /* if to print messages */
  n = tbg_cnt(tabag);           /* get the number of transactions */
  w = tbg_wgt(tabag);           /* and the new transaction weight */
  XMSG(stderr, "[%"TID_FMT, n); /* print number of transactions */
  if (w != (SUPP)n) { XMSG(stderr, "/%"SUPP_FMT, w); }
  XMSG(stderr, " transaction(s)] done [%.2fs].\n", SEC_SINCE(t));
  #endif
  return 0;                     /* return 'ok' */
}  /* acc_data() */

/*--------------------------------------------------------------------*/

int accret_report (ACCRET *accret, ISREPORT *report)
{                               /* --- prepare reporter for Accretion */
  TID    n;                     /* number of transactions */
  SUPP   w;                     /* total transaction weight */
  double smax;                  /* absolute maximum support */
  int    mrep;                  /* mode for item set reporter */

  assert(accret && report);     /* check the function arguments */
  accret->report = report;      /* note the item set reporter */

  /* --- get reporting mode --- */
  mrep = 0;                     /* set default reporting mode */
  #ifdef USE_ZLIB               /* if optional output compression */
  if (accret->mode & ACC_ZLIB)  /* if the compression flag is set, */
    mrep |= ISR_ZLIB;           /* transfer it to the report mode */
  #endif

  /* --- configure item set reporter --- */
  w = tbg_wgt(accret->tabag);   /* set support and size range */
  smax = (accret->smax < 0) ? -accret->smax
       : (accret->smax/100.0) *(double)w *(1-DBL_EPSILON);
  isr_setsupp(report, (RSUPP)accret->supp, (RSUPP)floorsupp(smax));
  isr_setsize(report, accret->zmin, accret->zmax);
  n = (accret->mode & ACC_PREFMT)/* get range of nums. to preformat */
    ? (TID)ib_maxfrq(tbg_base(accret->tabag)) : -1;
  if ((isr_prefmt(report, (TID)accret->supp, n) != 0)
  ||  (isr_settarg(report, ISR_ALL, mrep, -1)   != 0))
    return E_NOMEM;             /* set pre-format and target type */
  return 0;                     /* return 'ok' */
}  /* accret_report() */

/*--------------------------------------------------------------------*/

int accret_mine (ACCRET *accret, ITEM maxext)
{                               /* --- accretion algorithm */
  int     r;                    /* result of function call */
  #ifndef QUIET                 /* if to print messages */
  clock_t t;                    /* timer for measurements */
  #endif                        /* (only needed for messages) */

  assert(accret);               /* check the function arguments */
  CLOCK(t);                     /* start timer, print log message */
  XMSG(stderr, "writing %s ... ", isr_name(accret->report));
  accret->maxext = (maxext > 0) ? maxext : 1;
  r = accret_base(accret);      /* note maximum number of extensions */
  if (r < 0) return E_NOMEM;    /* and search for frequent item sets */
  XMSG(stderr, "[%"SIZE_FMT" set(s)]", isr_repcnt(accret->report));
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  return 0;                     /* return 'ok' */
}  /* accret_mine() */

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/
#ifdef ACC_MAIN

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

static int setbdr (ISREPORT *report, SUPP w, ITEM min,
                   double **border, ITEM n)
{                               /* --- set the support border */
  double s;                     /* to traverse the support values */

  assert(report                 /* check the function arguments */
  &&    (w > 0) && (n >= 0) && border && (*border || (n <= 0)));
  while (--n >= 0) {            /* traverse the support values */
    s = (*border)[n];           /* transform to absolute count */
    s = ceilsupp((s >= 0) ? s/100.0 *(double)w *(1-DBL_EPSILON) : -s);
    if (isr_setbdr(report, n+min, (RSUPP)s) < 0) return -1;
  }                             /* set support in item set reporter */
  if (*border) { free(*border); *border = NULL; }
  return 0;                     /* return 'ok' */
}  /* setbdr() */

/*--------------------------------------------------------------------*/

#ifndef NDEBUG                  /* if debug version */
  #undef  CLEANUP               /* clean up memory and close files */
  #define CLEANUP \
  if (accret) accret_delete(accret, 0); \
  if (twrite) twr_delete(twrite, 1);    \
  if (report) isr_delete(report, 0);    \
  if (tabag)  tbg_delete(tabag,  0);    \
  if (tread)  trd_delete(tread,  1);    \
  if (ibase)  ib_delete (ibase);        \
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
  CCHAR   *dflt    = " (%a,%4P)";    /* default format for check */
  CCHAR   *info    = dflt;      /* format for information output */
  int     target   = 'm';       /* target type (closed/maximal) */
  ITEM    zmin     = 2;         /* minimum size of an item set */
  ITEM    zmax     = ITEM_MAX;  /* maximum size of an item set */
  double  smin     = 1;         /* minimum support of an item set */
  int     stat     = 'p';       /* test statistic to use */
  int     sflgs    = 0;         /* test statistic flags */
  double  siglvl   = 1;         /* significance level (in percent) */
  ITEM    maxext   = 2;         /* maximum number of extension items */
  int     sort     = 2;         /* flag for item sorting and recoding */
  int     mode     = ACC_DEFAULT|ACC_PREFMT;    /* search mode */
  int     mtar     = 0;         /* mode for transaction reading */
  int     scan     = 0;         /* mode for item set reporting */
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
    printf("         (s: frequent, c: closed, m: maximal item sets)\n");
    printf("-m#      minimum number of items per item set     "
                    "(default: %"ITEM_FMT")\n", zmin);
    printf("-n#      maximum number of items per item set     "
                    "(default: no limit)\n");
    printf("-s#      minimum support of an item set           "
                    "(default: %g)\n", smin);
    printf("         (positive: percentage, "
                     "negative: absolute number)\n");
    printf("-e#      test statistic for item set evaluation   "
                    "(default: '%c')\n", stat);
    printf("-d#      significance level (maximum p-value)     "
                    "(default: %g%%)\n", siglvl);
    printf("-i       invalidate eval. below expected support  "
                    "(default: evaluate all)\n");
    printf("-x#      maximum number of extension items        "
                    "(default: %"ITEM_FMT")\n", maxext);
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
    printf("outfile  file to write found item sets to         "
                    "[optional]\n");
    return 0;                   /* print a usage message */
  }                             /* and abort the program */
  #endif  /* #ifndef QUIET */
  /* free option characters: acjlopuy [A-Z]\[CFNPRZ] */

  /* --- evaluate arguments --- */
  for (i = 1; i < argc; i++) {  /* traverse arguments */
    s = argv[i];                /* get option argument */
    if (optarg) { *optarg = s; optarg = NULL; continue; }
    if ((*s == '-') && *++s) {  /* -- if argument is an option */
      while (*s) {              /* traverse options */
        switch (*s++) {         /* evaluate switches */
          case '!': help();                          break;
          case 't': target = (*s) ? *s++ : 's';      break;
          case 'm': zmin   = (ITEM)strtol(s, &s, 0); break;
          case 'n': zmax   = (ITEM)strtol(s, &s, 0); break;
          case 's': smin   =       strtod(s, &s);    break;
          case 'e': stat   = (*s) ? *s++ : 'x';      break;
          case 'd': siglvl =       strtod(s, &s);    break;
          case 'i': sflgs  = ACC_INVBXS;             break;
          case 'x': maxext = (ITEM)strtol(s, &s, 0); break;
          case 'q': sort   = (int) strtol(s, &s, 0); break;
          case 'F': bdrcnt = getbdr(s, &s, &border); break;
          case 'R': optarg = &fn_sel;                break;
          case 'P': optarg = &fn_psp;                break;
          case 'Z': stats  = 1;                      break;
          case 'N': mode  &= ~ACC_PREFMT;            break;
          case 'g': scan   = 1;                      break;
          #ifdef USE_ZLIB       /* if optional output compression */
          case 'z': mode  |= ACC_ZLIB;               break;
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
  if (siglvl > 100) error(E_SIGLVL,  siglvl);
  if (bdrcnt < 0)   error(E_NOMEM);
  if ((!fn_inp || !*fn_inp) && (fn_sel && !*fn_sel))
    error(E_STDIN);             /* stdin must not be used twice */
  switch (target) {             /* check and translate target type */
    case 's': target = ISR_ALL;              break;
    case 'c': target = ISR_CLOSED;           break;
    case 'm': target = ISR_MAXIMAL;          break;
    default : error(E_TARGET, (char)target); break;
  }                             /* (get target type code) */
  switch (stat) {               /* check and translate target type */
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
  if (info == dflt)             /* adapt the default info. format */
    info = (smin < 0) ? " (%a,%4P)" : " (%3S,%4P)";
  if (maxext < 0)               /* a negative values means that */
    maxext = ITEM_MAX;          /* there is no limit on extensions */
  mode |= ACC_VERBOSE|ACC_NOCLEAN;
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

  /* --- find frequent item sets --- */
  accret = accret_create(target, smin, 100.0, zmin, zmax,
                         stat, siglvl, mode);
  if (!accret) error(E_NOMEM);  /* create an Accretion miner */
  k = accret_data(accret, tabag, sort);
  if (k) error(k);              /* prepare data for Accretion */
  report = isr_create(ibase);   /* create an item set reporter */
  if (!report) error(E_NOMEM);  /* and configure it */
  k = accret_report(accret, report);
  if (k) error(k);              /* prepare reporter for Accretion */
  if (setbdr(report, w, zmin, &border, bdrcnt) != 0)
    error(E_NOMEM);             /* set the support border */
  if (fn_psp && (isr_addpsp(report, NULL) < 0))
    error(E_NOMEM);             /* set a pattern spectrum if req. */
  if (isr_setfmt(report, scan, hdr, sep, NULL, info) != 0)
    error(E_NOMEM);             /* set the output format strings */
  k = isr_open(report, NULL, fn_out);
  if (k) error(k, isr_name(report)); /* open the item set file */
  if ((accret_report(accret, report) < 0)
  ||  (isr_setup(report) < 0))  /* prepare reporter for Accretion */
    error(E_NOMEM);             /* and set up the item set reporter */
  k = accret_mine(accret, maxext);
  if (k) error(k);              /* find frequent item sets */
  if (stats)                    /* print item set statistics */
    isr_prstats(report, stdout, 0);
  if (isr_close(report) != 0)   /* close the item set output file */
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

#endif
