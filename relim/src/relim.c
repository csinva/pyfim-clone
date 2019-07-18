/*----------------------------------------------------------------------
  File    : relim.c
  Contents: relim algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2004.11.05 file created from eclat.c
            2004.11.17 first reasonably fast version completed
            2004.11.18 start of loop over transactions lists improved
            2004.11.23 absolute/relative support output changed
            2004.12.09 filter added (binary logarithm of supp. quotient)
            2005.06.20 use of flag for "no item sorting" corrected
            2006.11.26 adapted to new structures ISFMTR and ISEVAL
            2007.02.13 adapted to modified module tabread
            2008.05.02 default limit for maximal number of items removed
            2008.10.13 adapted to modified module "tract", some redesign
            2008.10.15 simplification of the recursion parameters
            2008.10.29 reading item insertion penalties added
            2008.11.03 mining with insertion penalties added
            2008.11.11 adapted to insertion penalties in item base
            2008.11.13 adapted to changes in transaction management
            2008.11.19 sorting transaction list before processing added
            2008.12.05 perfect extension pruning added (optional)
            2008.12.11 special function for unlimited insertions added
            2009.05.28 adapted to modified function tbg_filter()
            2009.10.15 adapted to item set counter in reporter
            2009.10.16 closed and maximal item set mining added
            2010.03.04 sorting improved (appending the rest list)
            2010.03.10 projection memory combined into one block
            2010.03.15 bug in combined memory deallocation fixed
            2010.03.18 recording of actual number of occurrences added
            2010.04.07 threshold for both item set support and weight
            2010.07.14 output file made optional (for benchmarking)
            2010.08.19 item selection file added as optional input
            2010.08.22 adapted to modified modules tabread and tract
            2010.10.15 adapted to modified interface of module report
            2010.11.05 clearer interpretation of minimum support
            2010.11.24 adapted to modified error reporting (tract)
            2010.12.11 adapted to a generic error reporting function
            2011.03.16 closed/maximal item sets with item insertions
            2011.03.20 optional integer transaction weights added
            2011.05.30 item weight combination with t-norms added
            2011.06.02 initialize header table with memset()
            2011.07.08 adapted to modified function tbg_recode()
            2011.08.28 output of item set counters per size added
            2011.08.29 16 items machine added (without item insertions)
            2013.04.01 adapted to type changes in module tract
            2013.10.15 checks of return code of isr_report() added
            2013.10.18 optional pattern spectrum collection added
            2013.11.12 item insertion penalties changed to option -R#
            2014.05.12 option -F# added (support border for filtering)
            2014.08.02 option -c renamed to -i (min. supp. with insert.)
            2014.08.28 functions rel_data() and rel_report() added
            2014.10.24 changed from LGPL license to MIT license
            2016.02.19 added pre-formatting for some integer numbers
            2016.11.25 relim miner object and interface introduced
            2017.05.30 optional output compression with zlib added
            2017.07.04 tree-based algorithm variant added
            2017.07.06 use of 16-items machine in tree-based algorithm
------------------------------------------------------------------------
  Reference for the RElim algorithm:
    C. Borgelt.
    Keeping Things Simple:
    Finding Frequent Item Sets by Recursive Elimination.
    Workshop Open Source Data Mining Software
    (OSDM'05, Chicago, IL), 66--70.
    ACM Press, New York, NY, USA 2005
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
#ifdef RELIM_MAIN
#ifndef PSP_REPORT
#define PSP_REPORT
#endif
#ifndef TA_READ
#define TA_READ
#endif
#endif
#ifdef RELIM_ABORT
#include "sigint.h"
#endif
#include "memsys.h"
#include "relim.h"
#include "fim16.h"
#ifdef RELIM_MAIN
#include "error.h"
#endif
#ifdef STORAGE
#include "storage.h"
#endif

#ifndef INFINITY
#define INFINITY    (DBL_MAX+DBL_MAX)
#endif                          /* MSC still does not support C99 */

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define PRGNAME     "relim"
#define DESCRIPTION "find frequent item sets " \
                    "with a recursive elimination algorithm"
#define VERSION     "version 4.24 (2017.07.06)        " \
                    "(c) 2004-2017   Christian Borgelt"

/* --- error codes --- */
/* error codes   0 to  -4 defined in tract.h */
#define E_STDIN      (-5)       /* double assignment of stdin */
#define E_OPTION     (-6)       /* unknown option */
#define E_OPTARG     (-7)       /* missing option argument */
#define E_ARGCNT     (-8)       /* too few/many arguments */
#define E_TARGET     (-9)       /* invalid target type */
#define E_SIZE      (-10)       /* invalid item set size */
#define E_SUPPORT   (-11)       /* invalid minimum item set support */
#define E_WEIGHT    (-12)       /* invalid minimum transaction weight */
#define E_MEASURE   (-13)       /* invalid evaluation measure */
#define E_TNORM     (-14)       /* invalid triangular norm */
#define E_VARIANT   (-16)       /* invalid algorithm variant */
/* error codes -15 to -25 defined in tract.h */

#ifndef QUIET                   /* if not quiet version, */
#define MSG         fprintf     /* print messages */
#define XMSG        if (relim->mode & REL_VERBOSE) fprintf
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
typedef double TNORM (double a, double b);

typedef struct tsle {           /* --- trans. suffix list element --- */
  struct tsle   *succ;          /* successor element in list */
  const  ITEM   *items;         /* items in the transaction */
  SUPP          occ;            /* number of occurrences */
} TSLE;                         /* (transaction list element) */

typedef struct {                /* --- transaction suffix list --- */
  TSLE          *head;          /* head element of the list */
  SUPP          occ;            /* total number of occurrences */
} TSLIST;                       /* (transaction list) */

typedef struct txle {           /* --- trans. suffix list element --- */
  struct txle   *succ;          /* successor element in list */
  const  ITEM   *items;         /* items in the transaction */
  SUPP          occ;            /* number of actual occurrences */
  double        wgt;            /* weight of transactions */
} TXLE;                         /* (transaction list element) */

typedef struct {                /* --- transaction suffix list --- */
  TXLE          *head;          /* head element of the list */
  SUPP          occ;            /* number of occurrences */
  double        wgt;            /* total transaction weight */
} TXLIST;                       /* (transaction list) */

typedef struct tzle {           /* --- trans. suffix list element --- */
  struct tzle   *succ;          /* successor element in list */
  const  ITEM   *items;         /* items in the transaction */
  SUPP          occ;            /* number of actual occurrences */
  SUPP          cnt;            /* number of transactions */
  double        wgt;            /* weight per transaction */
} TZLE;                         /* (transaction list element) */

typedef struct {                /* --- transaction suffix list --- */
  TZLE          *head;          /* head element of the list */
  SUPP          occ;            /* number of actual occurrences */
  double        wgt;            /* total transaction weight */
} TZLIST;                       /* (transaction list) */

typedef struct renode {         /* --- recursive elimination node --- */
  ITEM          item;           /* associated item */
  SUPP          wgt;            /* weight of represented transactions */
  struct renode *sibling;       /* successor node in sibling list */
  struct renode *children;      /* list of child nodes */
} RENODE;                       /* (recursive elimination node) */

struct _relim {                 /* --- relim miner --- */
  int           target;         /* target type (e.g. closed/maximal) */
  double        smin;           /* minimum support of an item set */
  SUPP          supp;           /* minimum support of an item set */
  double        sins;           /* minimum support with insertions */
  double        isup;           /* minimum support with insertions */
  TNORM         *tnorm;         /* t-norm for comb. item penalties */
  double        twgt;           /* minimum transaction weight */
  ITEM          zmin;           /* minimum size of a rule/item set */
  ITEM          zmax;           /* maximum size of a rule/item set */
  int           eval;           /* additional evaluation measure */
  double        thresh;         /* threshold for evaluation measure */
  int           algo;           /* variant of relim algorithm */
  int           mode;           /* search mode (e.g. pruning) */
  ITEM          sort;           /* threshold for list sorting */
  TABAG         *tabag;         /* transaction bag/multiset */
  ISREPORT      *report;        /* item set reporter */
  MEMSYS        *mem;           /* memory system for tree version */
  FIM16         *fim16;         /* 16 items machine */
  void          *buf;           /* buffer for projection */
  ITEMBASE      *base;          /* underlying item base */
};                              /* (relim miner) */

/*----------------------------------------------------------------------
  Constants
----------------------------------------------------------------------*/
#if !defined QUIET && defined RELIM_MAIN
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
  /* E_WEIGHT  -12 */  "invalid minimum transaction weight %g",
  /* E_MEASURE -13 */  "invalid evaluation measure '%c'",
  /* E_TNORM   -14 */  "invalid triangular norm '%c'",
  /* E_NOITEMS -15 */  "no (frequent) items found",
  /* E_VARIANT -16 */  "invalid relim variant '%c'",
  /*           -17 */  "unknown error"
};
#endif

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
#ifdef RELIM_MAIN
#ifndef QUIET
static CCHAR    *prgname;       /* program name for error messages */
#endif
static TABREAD  *tread  = NULL; /* table/transaction reader */
static ITEMBASE *ibase  = NULL; /* item base */
static TABAG    *tabag  = NULL; /* transaction bag/multiset */
static ISREPORT *report = NULL; /* item set reporter */
static TABWRITE *twrite = NULL; /* table writer for pattern spectrum */
static double   *border = NULL; /* support border for filtering */
static RELIM    *relim  = NULL; /* relim miner object */
#endif

/*----------------------------------------------------------------------
  Auxiliary Functions (for debugging)
----------------------------------------------------------------------*/
#if !defined NDEBUG && defined RELIM_MAIN

static void indent (int k)
{ while (--k >= 0) printf("   "); }

/*--------------------------------------------------------------------*/

#define SHOW(show,type,weight,fmt,term) \
static void show (const char *text, type *list, int ind) \
{                               /* --- show a transaction list */      \
  ITEM       n;                 /* item counter */                     \
  const ITEM *p;                /* to traverse the item array */       \
  weight     w;                 /* total weight */                     \
                                                                       \
  indent(ind);                  /* indent the output line */           \
  printf("%s\n", text);         /* print the given text */             \
  for (w = 0, n = 0; list; list = list->succ){ /* traverse the list */ \
    indent(ind);                /* indent the output line */           \
    OCCUR;                      /* print the number of occurrences */  \
    printf(fmt" :", term);      /* print the transaction weight */     \
    for (p = list->items; *p >= 0; p++)                                \
      printf(" %s", ib_name(ibase, *p));                               \
    printf("\n");               /* print items in the transaction */   \
    w += term; n++;             /* sum the transaction weights */      \
  }                             /* and count the transactions */       \
  indent(ind);                  /* indent the output line */           \
  printf("total: %"ITEM_FMT"/"fmt"\n", n, w);                          \
}  /* show() */                 /* print total transaction weight */

/*--------------------------------------------------------------------*/

#define OCCUR
SHOW(show,     TSLE, SUPP, "%"SUPP_FMT, list->occ)
#undef  OCCUR
#define OCCUR  printf("%3"SUPP_FMT"/", list->occ)
SHOW(show_ins, TXLE, double, "%5.2f", list->wgt)
#undef  OCCUR
#define OCCUR  printf("%3"SUPP_FMT"/%3"TID_FMT"/%5.2f/", \
                      list->occ, list->cnt, list->wgt)
SHOW(show_lim, TZLE, double, "%5.2f", list->wgt *(double)list->cnt)

/*--------------------------------------------------------------------*/

static void show_tree (RENODE *node, int indent, ITEMBASE *ibase)
{                               /* --- show a recursive elim. tree */
  int i;                        /* loop variable */

  for ( ; node; node = node->sibling) {
    for (i = 0; i < indent; i++)/* traverse the sibling list and */
      printf("  ");             /* print indentation for each node */
    if (ibase) {                /* if an item base is given */
      if      (node->item == -1) printf("*/*");
      else if (node->item <   0) printf("%02x", node->item & 0xff);
      else printf("%s/%d", ib_name(ibase, node->item), node->item); }
    else {                      /* print item name and id, otherwise */
      if      (node->item == -1) printf("*");
      else if (node->item <   0) printf("%02x", node->item & 0xff);
      else printf("%d", node->item);
    }                           /* print only the item id */
    printf(":%d\n", node->wgt); /* print the transactin weight */
    if (node->children) show_tree(node->children, indent+1, ibase);
  }                             /* recursively print the subtrees */
}  /* show() */

#endif  /* #ifndef NDEBUG */
/*----------------------------------------------------------------------
  Triangular Norms (t-norms)
----------------------------------------------------------------------*/

static double t_min  (double a, double b)
{ return (a < b) ? a : b; }

static double t_nilp (double a, double b)
{ return (a+b <= 1) ? 0 : (a < b) ? a : b; }

static double t_prod (double a, double b)
{ return a*b; }

static double t_luka (double a, double b)
{ double x = a+b-1; return (x > 0) ? x : 0; }

static double t_hama (double a, double b)
{ double x = a+b-a*b; return (x > 0) ? (a*b)/x : 0; }

/*--------------------------------------------------------------------*/

static TNORM *tnorms[] = {      /* t-norms (triangular norms) */
  /* REL_MIN  0 */  t_min,      /* minimum */
  /* REL_NILP 1 */  t_nilp,     /* nil-potent minimum */
  /* REL_PROD 2 */  t_prod,     /* product */
  /* REL_LUKA 3 */  t_luka,     /* Lukasiewicz */
  /* REL_HAMA 4 */  t_hama,     /* Hamacher product */
};

/*----------------------------------------------------------------------
  Comparing and Sorting Transaction (Suffix) Lists
----------------------------------------------------------------------*/

static int cmp (const ITEM *a, const ITEM *b)
{                               /* --- compare two transactions */
  assert(a && b);               /* check the function arguments */
  for ( ; 1; a++, b++) {        /* lexicographic comparison loop */
    if (*a < *b) return -1;     /* compare corresponding items */
    if (*a > *b) return +1;     /* and if one is greater, abort */
    if (*a <= TA_END) return 0; /* otherwise check for the sentinel */
  }                             /* and abort if it is reached */
}  /* cmp() */

/*--------------------------------------------------------------------*/

#define SORT(sort,type) \
static type* sort (type *list) \
{                               /* --- sort a transaction list */      \
  type *b, *a = list;           /* to traverse the lists */            \
  type **e;                     /* end of the output list */           \
                                                                       \
  for (b = list->succ; b; ) {   /* traverse the list to sort */        \
    b = b->succ;                /* two steps on a, one step on list */ \
    if (b) { b = b->succ; list = list->succ; }                         \
  }                             /* (split list into two halves) */     \
  b = list->succ;               /* get the second list and */          \
  list->succ = NULL;            /* terminate the first list */         \
  if (a->succ) a = sort(a);     /* if longer than one element, */      \
  if (b->succ) b = sort(b);     /* sort the lists recursively */       \
  for (e = &list; 1; ) {        /* transaction list merge loop */      \
    register int c = cmp(a->items, b->items);                          \
    if      (c < 0) { *e = a; e = &a->succ; if (!(a = *e)) break; }    \
    else if (c > 0) { *e = b; e = &b->succ; if (!(b = *e)) break; }    \
    COMBINE                     /* copy unique transactions and */     \
  }                             /* combine equal transactions */       \
  *e = (a) ? a : b;             /* append the non-empty list */        \
  return list;                  /* return the sorted list */           \
}  /* sort() */

/*--------------------------------------------------------------------*/

#define COMBINE \
    else                      { a->occ += b->occ; b = b->succ; \
                                *e = a; e = &a->succ; a = *e;  \
                                if (!a || !b) break; }
SORT(sort,     TSLE)

#undef  COMBINE
#define COMBINE \
    else                      { a->occ += b->occ;              \
                                a->wgt += b->wgt; b = b->succ; \
                                *e = a; e = &a->succ; a = *e;  \
                                if (!a || !b) break; }
SORT(sort_ext, TXLE)

#undef  COMBINE
#define COMBINE \
    else if (a->wgt < b->wgt) { *e = a; e = &a->succ;          \
                                if (!(a = *e)) break; }        \
    else if (a->wgt > b->wgt) { *e = b; e = &b->succ;          \
                                if (!(b = *e)) break; }        \
    else                      { a->occ += b->occ;              \
                                a->cnt += b->cnt; b = b->succ; \
                                *e = a; e = &a->succ; a = *e;  \
                                if (!a || !b) break; }
SORT(sort_wgt, TZLE)

/*----------------------------------------------------------------------
  Recursive Elimination: Basic Version
----------------------------------------------------------------------*/

static int recurse (RELIM *relim, TSLIST *lists, ITEM k, TID n)
{                               /* --- recursive elimination (arrays) */
  int    r;                     /* status indicator */
  TSLIST *proj  = NULL;         /* list(s) of (projected) database */
  TSLE   *elems = NULL;         /* to traverse the transaction lists */
  TSLIST *cur, *tal;            /* current/projected transaction list */
  TSLE   *src, *dst;            /* to traverse the transaction lists */
  SUPP   pex;                   /* minimum support for perfect exts. */

  assert(relim && lists && (k > 0)); /* check the function arguments */
  #ifdef RELIM_ABORT            /* if to check for interrupt */
  if (sig_aborted()) return -1; /* if execution was aborted, */
  #endif                        /* abort the recursion */
  if ((k > 1)                   /* if there is more than one item and */
  &&  isr_xable(relim->report, 2)) {     /* another item can be added */
    proj = (TSLIST*)malloc((size_t)(k-1) *sizeof(TSLIST)
                          +(size_t) n    *sizeof(TSLE));
    if (!proj) return -1;       /* allocate list and element arrays */
    elems = (TSLE*)(proj +k-1); /* and organize the memory */
    memset(proj, 0, (size_t)(k-1) *sizeof(TSLIST));
  }                             /* initialize the projection header */
  pex = (relim->mode & REL_PERFECT)
      ? isr_supp(relim->report) : SUPP_MAX;
  for (r = 0; --k >= 0; ) {     /* traverse the transaction lists */
    cur = lists +k;             /* get the next transaction list */
    if      (cur->occ >= pex)   /* if item is a perfect extension, */
      isr_addpex(relim->report, k);/* add it to the item set reporter */
    else if (cur->occ >= relim->supp) {  /* if support is high enough */
      r = isr_add(relim->report, k, cur->occ);
      if (r < 0) break;         /* add current item to the reporter */
      if (r > 0) {              /* if the item needs processing */
        if (cur->head && proj){ /* if another item can be added */
          if (cur->head->succ   /* if list has more than one element */
          && (k <= relim->sort))/* and there are few enough items, */
            cur->head = sort(cur->head);   /* sort the trans. list */
          dst = elems;          /* traverse list for current item */
          for (src = cur->head; src; src = src->succ) {
            tal = proj +*src->items;  /* get the first item and */
            tal->occ  += src->occ;    /* sum the transaction weight */
            if (src->items[1] < 0) continue;
            dst->items = src->items+1;
            dst->occ   = src->occ;    /* copy the item array and */
            dst->succ  = tal->head;   /* the number of occurrences */
            tal->head  = dst++; /* add the new element at the head */
          }                     /* of the corresponding list */
          r = recurse(relim, proj, k, (TID)(dst-elems));
          if (r < 0) break;     /* find frequent item sets */
        }                       /* recursively in the projection */
        r = isr_report(relim->report);
        if (r < 0) break;       /* report the current item set */
        isr_remove(relim->report, 1);
      }                         /* remove the current item */
    }                           /* from the item set reporter */
    cur->occ = 0;               /* clear the number of occurrences */
    while (cur->head) {         /* while the list is not empty, */
      src       = cur->head;    /* remove the first element and */
      cur->head = src->succ;    /* get the list of the first item */
      tal = lists +*src->items++;
      tal->occ += src->occ;     /* sum the number of occurrences */
      if (*src->items < 0) continue;
      src->succ = tal->head;    /* reassign the transactions */
      tal->head = src;          /* based on their first items */
    }                           /* and skip this first item */
  }
  if (proj) free(proj);         /* delete the list and element arrays */
  return r;                     /* return the error status */
}  /* recurse() */

/*--------------------------------------------------------------------*/

int relim_base (RELIM *relim)
{                               /* --- recursive elimination (arrays) */
  int    r;                     /* result of recursion */
  ITEM   i, k;                  /* loop variable, number of items */
  TID    n;                     /* number of transactions */
  TRACT  *t;                    /* to traverse the transactions */
  TSLIST *lists, *tal;          /* (array of) transaction list(s) */
  TSLE   *elems, *dst;          /* (array of) trans. list element(s) */

  assert(relim);                /* check the function arguments */
  if (tbg_wgt(relim->tabag) < relim->supp)
    return 0;                   /* check the total transaction weight */
  k = tbg_itemcnt(relim->tabag);/* get and check the number of items */
  if (k <= 0) return isr_report(relim->report);
  n = tbg_cnt(relim->tabag);    /* get the number of transactions */
  lists = (TSLIST*)malloc((size_t)k *sizeof(TSLIST)
                         +(size_t)n *sizeof(TSLE));
  if (!lists) return -1;        /* allocate lists and element arrays */
  dst = elems = (TSLE*)(lists +k);       /* and initialize the lists */
  memset(lists, 0, (size_t)k *sizeof(TSLIST));
  while (--n >= 0) {            /* traverse the transactions */
    t = tbg_tract(relim->tabag, n); /* get the current transaction */
    dst->items = ta_items(t);       /* and its item array */
    i = *dst->items++;          /* get the first item and skip it */
    if (i < 0) continue;        /* skip empty transactions */
    tal = lists +i;             /* otherwise sum transaction weight */
    tal->occ += dst->occ = ta_wgt(t);
    if (*dst->items < 0) continue;
    dst->succ = tal->head;      /* skip one element transactions */
    tal->head = dst++;          /* add the new element to the */
  }                             /* list for the first item */
  r = recurse(relim, lists, k, (TID)(dst-elems));
  free(lists);                  /* execute recursive elimination */
  if (r >= 0)                   /* finally report the empty item set */
    r = isr_report(relim->report);
  return r;                     /* return the error status */
}  /* relim_base() */

/*----------------------------------------------------------------------
  Recursive Elimination with 16 Items Machine
----------------------------------------------------------------------*/

static int rec_m16 (RELIM *relim, TSLIST *lists, ITEM k, TID n)
{                               /* --- recursive elimination (arrays) */
  int    r;                     /* status indicator */
  ITEM   i;                     /* item buffer */
  TSLIST *proj  = NULL;         /* list(s) of (projected) database */
  TSLE   *elems = NULL;         /* to traverse the transaction lists */
  TSLIST *cur, *tal;            /* current/projected transaction list */
  TSLE   *src, *dst;            /* to traverse the transaction lists */
  SUPP   pex;                   /* minimum support for perfect exts. */

  assert(relim && lists && (k > 0)); /* check the function arguments */
  #ifdef RELIM_ABORT            /* if to check for interrupt */
  if (sig_aborted()) return -1; /* if execution was aborted, */
  #endif                        /* abort the recursion */
  if ((k > 1)                   /* if there is more than one item and */
  &&  isr_xable(relim->report, 2)) {     /* another item can be added */
    proj = (TSLIST*)malloc((size_t)(k-1) *sizeof(TSLIST)
                          +(size_t)(n+1) *sizeof(TSLE));
    if (!proj) return -1;       /* allocate list and element arrays */
    elems = (TSLE*)(proj +k-1); /* and organize the memory */
    memset(proj, 0, (size_t)(k-1) *sizeof(TSLIST));
  }                             /* initialize the projection header */
  pex = (relim->mode & REL_PERFECT)
      ? isr_supp(relim->report) : SUPP_MAX;
  for (r = 0; --k >= 16; ) {    /* traverse the transaction lists */
    cur = lists +k;             /* get the next transaction list */
    if      (cur->occ >= pex)   /* if item is a perfect extension, */
      isr_addpex(relim->report, k);/* add it to the item set reporter */
    else if (cur->occ >= relim->supp) {    /* if the item is frequent */
      r = isr_add(relim->report, k, cur->occ);
      if (r < 0) break;         /* add current item to the reporter */
      if (r > 0) {              /* if the item needs processing */
        if (cur->head && proj){ /* if another item can be added */
          if (cur->head->succ   /* if list has more than one element */
          && (k <= relim->sort))/* and there are few enough items, */
            cur->head = sort(cur->head);   /* sort the trans. list */
          dst = elems;          /* traverse list for current item */
          for (src = cur->head; src; src = src->succ) {
            i = *src->items;    /* get the first item in suffix */
            if (i < 0) {        /* if packed items (bit represent.) */
              proj->occ += dst->occ = src->occ;
              dst->items = src->items;
              dst->succ  = proj->head;
              proj->head = dst++;
              continue;         /* transactions with packed items */
            }                   /* are stored in the first list */
            tal = proj +i;      /* get the destination list and */
            tal->occ += src->occ; /* sum the transaction weight */
            if (src->items[1] <= TA_END) continue;
            dst->items = src->items+1;
            dst->occ   = src->occ;    /* store trans. suffix and */
            dst->succ  = tal->head;   /* the number of occurrences */
            tal->head  = dst++; /* add the new element at the head */
          }                     /* of the corresponding list */
          r = rec_m16(relim, proj, k, (TID)(dst-elems));
          if (r < 0) break;     /* find frequent item sets */
        }                       /* recursively in the projection */
        r = isr_report(relim->report);
        if (r < 0) break;       /* report the current item set */
        isr_remove(relim->report, 1);
      }                         /* remove the current item */
    }                           /* from the item set reporter */
    cur->occ = 0;               /* clear the number of occurrences */
    while (cur->head) {         /* while the list is not empty, */
      src       = cur->head;    /* remove the first element */
      cur->head = src->succ;    /* and get the first item */
      i = *src->items;          /* get first item of trans. suffix */
      if (i < 0) {              /* if packed items (bit represent.) */
        lists->occ += src->occ; /* sum the transaction weights */
        src->succ   = lists->head;
        lists->head = src;      /* transactions with packed items */
        continue;               /* are stored in the first list */
      }                         /* (catchment basin for 16 items) */
      tal = lists +i;           /* get transaction list for item */
      tal->occ += src->occ;     /* sum the number of occurrences */
      if (*++src->items <= TA_END) continue;
      src->succ = tal->head;    /* reassign the transactions */
      tal->head = src;          /* based on their first items */
    }                           /* and skip this first item */
  }
  if ((r >= 0) && (lists->occ >= relim->supp)) {
    for (src = lists->head; src; src = src->succ)
      m16_add(relim->fim16, (BITTA)(src->items[0] & ~TA_END), src->occ);
    r = m16_mine(relim->fim16); /* traverse and add packed items */
  }                             /* and mine with 16 items machine */
  lists->head = NULL; lists->occ = 0;
  if (proj) free(proj);         /* delete the list and element arrays */
  return r;                     /* return the error status */
}  /* rec_m16() */

/*--------------------------------------------------------------------*/

int relim_m16 (RELIM *relim)
{                               /* --- recursive elimination (arrays) */
  int    r;                     /* result of recursion */
  ITEM   i, k;                  /* loop variable, number of items */
  TID    n;                     /* number of transactions */
  TRACT  *t;                    /* to traverse the transactions */
  TSLIST *lists, *tal;          /* (array of) transaction list(s) */
  TSLE   *elems, *dst;          /* (array of) trans. list element(s) */

  assert(relim);                /* check the function arguments */
  if (tbg_wgt(relim->tabag) < relim->supp)
    return 0;                   /* check the total transaction weight */
  k = tbg_itemcnt(relim->tabag);/* get and check the number of items */
  if (k <= 0) return isr_report(relim->report);
  n = tbg_cnt(relim->tabag);    /* get the number of transactions */
  lists = (TSLIST*)malloc((size_t)k *sizeof(TSLIST)
                         +(size_t)n *sizeof(TSLE));
  if (!lists) return -1;        /* allocate lists and element arrays */
  dst = elems = (TSLE*)(lists +k);       /* and initialize the lists */
  memset(lists, 0, (size_t)k *sizeof(TSLIST));
  relim->fim16 = m16_create(-1, relim->supp, relim->report);
  if (!relim->fim16) { free(lists); return -1; }
                                /* create a 16 items machine */
  while (--n >= 0) {            /* traverse the transactions */
    t = tbg_tract(relim->tabag, n); /* get the current transaction */
    dst->items = ta_items(t);       /* and its item array */
    i = *dst->items;            /* get the first item */
    if (i <= TA_END) continue;  /* skip empty transactions */
    if (i <  0) {               /* if packed items (bit represent.) */
      lists->occ += dst->occ = ta_wgt(t);
      dst->succ   = lists->head;
      lists->head = dst++;      /* transactions with packed items */
      continue;                 /* are stored in the first list */
    }                           /* (catchment basin for 16 items) */
    tal = lists +i;             /* otherwise sum transaction weight */
    tal->occ += dst->occ = ta_wgt(t);
    if (*++dst->items <= TA_END) continue;
    dst->succ = tal->head;      /* skip one element transactions */
    tal->head = dst++;          /* add the new element to the */
  }                             /* list for the first item */
  r = rec_m16(relim, lists, k, (TID)(dst -elems));
  m16_delete(relim->fim16);     /* execute recursive elimination */
  free(lists);                  /* and deallocate working memory */
  if (r >= 0)                   /* finally report the empty item set */
    r = isr_report(relim->report);
  return r;                     /* return the error status */
}  /* relim_m16() */

/*----------------------------------------------------------------------
  Recursive Elimination: Tree-based Version
----------------------------------------------------------------------*/

static int add (RENODE *node, MEMSYS *mem, const ITEM *items, SUPP wgt)
{                               /* --- add a transaction to the tree */
  ITEM   i;                     /* buffer for an item */
  RENODE **p;                   /* to traverse a child node list */

  assert(node                   /* check the function arguments */
  &&     mem && items && (wgt >= 0));
  while (1) {                   /* traverse the items of the set */
    node->wgt += wgt;           /* update the item set support */
    i = *items++;               /* get the next item in the set */
    if (i == TA_END) return 0;  /* if all items are processed, abort */
    p = &node->children;        /* traverse the list of children */
    while (*p && ((*p)->item > i)) p = &(*p)->sibling;
    if (!(node = *p) || (node->item != i)) break;
  }                             /* find the item/insertion position */
  node = (RENODE*)ms_alloc(mem);
  if (!node) return -1;         /* create a new prefix tree node */
  node->item    = i;            /* store the current item and */
  node->wgt     = wgt;          /* the support of the item set */
  node->sibling = *p;           /* insert the created node */
  *p = node;                    /* into the sibling list */
  while ((i = *items++) != TA_END) { /* traverse the remaining items */
    node = node->children = (RENODE*)ms_alloc(mem);
    if (!node) return -1;       /* create a new prefix tree node */
    node->item    = i;          /* store the next item and */
    node->wgt     = wgt;        /* the support of the item set */
    node->sibling = NULL;       /* there are no siblings yet */
  }                             /* insert node into the item list */
  node->children = NULL;        /* last created node is a leaf */
  return 1;                     /* return that nodes were added */
}  /* add() */

/*--------------------------------------------------------------------*/

static RENODE* copy (const RENODE *src, MEMSYS *mem)
{                               /* --- copy a c/m subtree */
  RENODE *dst;                  /* created copy of the source node */
  RENODE *child;                /* to traverse the child nodes */
  RENODE **end;                 /* end of the created child node list */

  assert(src && mem);           /* check the function arguments */
  dst = ms_alloc(mem);          /* create a copy of */
  if (!dst) return NULL;        /* the given source node */
  dst->item    = src->item;     /* transfer the associated item */
  dst->wgt     = src->wgt;      /* and the transaction weight */
  dst->sibling = NULL;          /* there are no sibling nodes */
  end = &dst->children;         /* get the end of the child node list */
  for (child = src->children; child; child = child->sibling) {
    *end = copy(child, mem);    /* traverse and copy the children, */
    if (!*end) return NULL;     /* append them to the sibling list, */
    end = &(*end)->sibling;     /* and advance the end of the list */
  }                             /* of copied child nodes */
  *end = NULL;                  /* terminate copied child node list */
  return dst;                   /* return the created copy */
}  /* copy() */

/*--------------------------------------------------------------------*/

static RENODE* merge (RENODE *s1, RENODE *s2, MEMSYS *mem)
{                               /* --- merge two node list */
  RENODE *out, **end;           /* output node list and end pointer */

  assert(mem);                  /* check the function arguments */
  if (!s1) return s2;           /* if there is only one node list, */
  if (!s2) return s1;           /* simply return the other list */
  end = &out;                   /* start the output list */
  while (1) {                   /* node list merge loop */
    if      (s1->item > s2->item) {
      *end = s1; end = &s1->sibling; s1 = *end; if (!s1) break; }
    else if (s1->item < s2->item) {
      *end = s2; end = &s2->sibling; s2 = *end; if (!s2) break; }
    else {                      /* copy nodes with singular items */
      s1->children = merge(s1->children, s2->children, mem);
      s1->wgt     += s2->wgt;   /* merge the children recursively */
      s2 = s2->sibling; *end = s1; end = &s1->sibling; s1 = *end;
      if (!s1 || !s2) break;    /* move node from the first source */
    }                           /* to the output and delete the one */
  }                             /* from the second source */
  *end = (s1) ? s1 : s2;        /* append the remaining nodes */
  return out;                   /* return the merged prefix tree */
}  /* merge() */

/*--------------------------------------------------------------------*/

static int rec_tree (RELIM *relim, RENODE *tree)
{                               /* --- recursive elimination (arrays) */
  int    r = 0;                 /* status/error indicator */
  RENODE *proj;                 /* tree of projected trans. database */
  SUPP   pex;                   /* minimum support for perfect exts. */

  assert(relim && tree);        /* check the function arguments */
  #ifdef RELIM_ABORT            /* if to check for interrupt */
  if (sig_aborted()) return -1; /* if execution was aborted, */
  #endif                        /* abort the recursion */
  pex = (relim->mode & REL_PERFECT)
      ? isr_supp(relim->report) : SUPP_MAX;
  while (tree                   /* while the tree is not empty and */
  &&    (tree->item >= 0)) {    /* the packed items are not reached */
    ms_push(relim->mem);        /* store the current memory state */
    proj = copy(tree, relim->mem);/* copy and merge first subtree */
    tree = merge(tree->children, tree->sibling, relim->mem);
    if (proj->wgt >= pex)       /* collect perfect extensions */
      isr_addpex(relim->report, proj->item);
    else if (proj->wgt >= relim->supp) {  /* process frequent items */
      r = isr_add(relim->report, proj->item, proj->wgt);
      if (r < 0) return -1;     /* add current item to the reporter */
      if (r > 0) {              /* if the item needs processing */
        if (proj->children      /* recursively find freq. item sets */
        &&  isr_xable(relim->report, 1))
          r = rec_tree(relim, proj->children);
        if (r >= 0) r = isr_report(relim->report);
        isr_remove(relim->report, 1);
        if (r <  0) return -1;  /* report the current item set */
      }                         /* and remove the current item */
    }                           /* from the item set reporter */
    ms_pop(relim->mem);         /* restore the memory mgmt. state */
  }                             /* (delete nodes of projection) */
  if (tree) {                   /* if there are packed items left */
    for ( ; tree; tree = tree->sibling)
      m16_add(relim->fim16, (BITTA)(tree->item & ~TA_END), tree->wgt);
    r = m16_mine(relim->fim16); /* traverse and add packed items */
  }                             /* and mine with 16 items machine */
  return r;                     /* return the error status */
}  /* rec_tree() */

/*--------------------------------------------------------------------*/

int relim_tree (RELIM *relim)
{                               /* --- recursive elimination (tree) */
  int    r = 0;                 /* result of recursion */
  ITEM   k;                     /* number of items */
  TID    n;                     /* number of transactions */
  TRACT  *t;                    /* to traverse the transactions */
  RENODE *tree;                 /* recursive elimination tree */

  assert(relim);                /* check the function arguments */
  if (tbg_wgt(relim->tabag) < relim->supp)
    return 0;                   /* check the total transaction weight */
  k = tbg_itemcnt(relim->tabag);/* get and check the number of items */
  if (k <= 0) return isr_report(relim->report);
  relim->mem = ms_create(sizeof(RENODE), 65535);
  if (!relim->mem) return -1;   /* create memory system for the nodes */
  tree = ms_alloc(relim->mem);  /* and with it a root node */
  tree->item = -1; tree->wgt = 0;
  tree->sibling = tree->children = NULL;
  relim->fim16 = NULL;          /* default: no 16-items machine */
  if (relim->mode & REL_FIM16){ /* if to use a 16-items machine */
    relim->fim16 = m16_create(-1, relim->supp, relim->report);
    if (!relim->fim16) { ms_delete(relim->mem); return -1; }
  }                             /* create a 16-items machine */
  for (n = tbg_cnt(relim->tabag); --n >= 0; ) {
    t = tbg_tract(relim->tabag, n);
    r = add(tree, relim->mem, ta_items(t), ta_wgt(t));
    if (r < 0) break;           /* traverse the transactions */
  }                             /* and add them to the tree */
  if (r >= 0)                   /* find frequent item sets */
    r = rec_tree(relim, tree->children);
  if (r >= 0)                   /* finally report the empty item set */
    r = isr_report(relim->report);
  if (relim->fim16)             /* if a 16-items machine was used, */
    m16_delete(relim->fim16);   /* delete the 16-items machine */
  ms_delete(relim->mem);        /* delete the memory mgmt. system */
  return r;                     /* return the error status */
}  /* relim_tree() */

/*----------------------------------------------------------------------
  Recursive Elimination: Unlimited Item Insertions
----------------------------------------------------------------------*/

static int rec_ins (RELIM *relim, TXLIST *lists, ITEM k, TID n)
{                               /* --- recursive elimination */
  int    r;                     /* error status */
  ITEM   i;                     /* current item */
  TXLIST *proj  = NULL;         /* list(s) of (projected) database */
  TXLE   *elems = NULL;         /* to traverse the transaction lists */
  TXLIST *cur, *tal;            /* current/projected transaction list */
  TXLE   *src, *dst;            /* to traverse the transaction lists */
  double pen, wgt;              /* insertion penalty, trans. weight */
  double pex;                   /* minimum weight for perfect exts. */

  assert(relim && lists && (k > 0)); /* check the function arguments */
  #ifdef RELIM_ABORT            /* if to check for interrupt */
  if (sig_aborted()) return -1; /* if execution was aborted, */
  #endif                        /* abort the recursion */
  if ((k > 1)                   /* if there is more than one item and */
  &&  isr_xable(relim->report, 2)) {     /* another item can be added */
    proj = (TXLIST*)malloc((size_t)k *sizeof(TXLIST)
                          +(size_t)n *sizeof(TXLE));
    if (!proj) return -1;       /* allocate list and element arrays */
    elems = (TXLE*)(proj +k);   /* and organize the memory */
    memset(proj, 0, (size_t)k *sizeof(TXLIST));
  }                             /* initialize the projection header */
  pex = (relim->mode & REL_PERFECT)
      ? isr_wgt(relim->report) : INFINITY;
  for (r = 0; --k >= 0; ) {     /* traverse the transaction lists */
    if (proj) {                 /* clear list of empty transactions */
      proj->head = NULL; proj->wgt = 0; proj->occ = 0; }
    dst = elems;                /* init. the projected database */
    cur = lists +k+1;           /* get the current transaction list */
    pen = ib_getpen(relim->base, k);
    if (pen > 0) {              /* if insertion penalty is positive */
      for (i = k+1; --i >= 0;){ /* traverse the preceding lists */
        for (src = lists[i].head; src; src = src->succ) {
          cur->wgt  += wgt = relim->tnorm(src->wgt, pen);
          if (!dst) continue;   /* sum the transaction weight */
          tal        = proj +i; /* for current and projected list */
          tal->wgt  += dst->wgt = wgt;
          dst->occ   = 0;       /* no occurrences (item is missing) */
          dst->items = src->items;
          dst->succ  = tal->head;
          tal->head  = dst++;   /* keep the item array unchanged */
        }                       /* and add the new list element to */
      }                         /* the corresp. transaction list */
    } /* The above loop inserts the current item into transactions  */
      /* that do not contain it. It is the only difference to the   */
      /* standard version of the algorithm apart from the fact that */
      /* empty transactions are processed in an additional list.    */
    i = -1;                     /* default: no recursion */
    if (cur->wgt >= pex) {      /* if item is a perfect extension, */
      isr_addpex(relim->report, k);      /* add it to the reporter */
      i = 0; }                  /* and clear the projection flag */
    else if ((cur->occ >= relim->supp)    /* if both support values */
    &&       (cur->wgt >= relim->isup)) { /* are large enough */
      r = isr_addwgt(relim->report, k, cur->occ, cur->wgt);
      if (r < 0) break;         /* add current item to the reporter */
      if (r > 0) {              /* if item needs processing */
        if ((k > 0) && proj) {  /* if another item can be added */
          if (cur->head         /* if list has more than one element */
          &&  cur->head->succ   /* and there are few enough items, */
          &&  (k <= relim->sort))     /* sort the transaction list */
            cur->head = sort_ext(cur->head);
          for (src = cur->head; src; src = src->succ) {
            i = *src->items +1; /* get and check the first item */
            if (i <= 0) i = 0;  /* get first list for sentinels */
            tal = proj +i;      /* sum the transaction weights */
            tal->occ  += dst->occ = src->occ;
            tal->wgt  += dst->wgt = src->wgt;
            dst->items = src->items +((i > 0) ? 1 : 0);
            dst->succ  = tal->head;
            tal->head  = dst++; /* add a new list element to */
          }                     /* the corresp. transaction list */
          r = rec_ins(relim, proj, k, (TID)(dst-elems));
          if (r < 0) break;     /* find frequent item sets */
        }                       /* recursively in the projection */
        r = isr_report(relim->report);
        if (r < 0) break;       /* report the current item set */
        isr_remove(relim->report, 1);
      }                         /* remove the current item */
    }                           /* from the item set reporter */
    if ((i < 0) && proj)        /* clear the projected database */
      memset(proj, 0, (size_t)k *sizeof(TXLIST));
    cur->wgt = 0; cur->occ = 0; /* clear the transaction weight */
    while (cur->head) {         /* while the list is not empty */
      src       = cur->head;    /* remove the first element */
      cur->head = src->succ;    /* from the transaction list */
      i = *src->items +1;       /* get and check the first item */
      if (i <= 0) i = 0;        /* if sentinel, get the first list, */
      else src->items++;        /* otherwise skip the item */
      tal       = lists +i;     /* get the destination list */
      tal->occ += src->occ;     /* sum the number of occurrences */
      tal->wgt += src->wgt;     /* and the transaction weight */
      src->succ = tal->head;    /* reassign the transactions */
      tal->head = src;          /* based on their first items */
    }                           /* and skip this first item */
  }
  if (proj) free(proj);         /* delete the list and element arrays */
  return r;                     /* return the error status */
}  /* rec_ins() */

/*--------------------------------------------------------------------*/

int relim_ins (RELIM *relim)
{                               /* --- recursive elimination  */
  int    r;                     /* result of recursion */
  ITEM   i, k;                  /* loop variable, number of items */
  TID    n;                     /* number of transactions */
  TRACT  *t;                    /* to traverse the transactions */
  TXLIST *lists, *tal;          /* (array of) transaction list(s) */
  TXLE   *elems, *dst;          /* (array of) trans. list element(s) */

  assert(relim);               /* check the function arguments */
  if (tbg_wgt(relim->tabag) < relim->supp)
    return 0;                   /* check the total transaction weight */
  k = tbg_itemcnt(relim->tabag);/* get and check the number of items */
  if (k <= 0) return isr_report(relim->report);
  n = tbg_cnt(relim->tabag);    /* get the number of transactions */
  lists = (TXLIST*)malloc((size_t)(k+1) *sizeof(TXLIST)
                         +(size_t) n    *sizeof(TXLE));
  if (!lists) return -1;        /* allocate lists and element arrays */
  dst = elems = (TXLE*)(lists +k+1);     /* and initialize the lists */
  memset(lists, 0, (size_t)(k+1) *sizeof(TXLIST));
  while (--n >= 0) {            /* traverse the transactions */
    t = tbg_tract(relim->tabag, n); /* initialize a new list element */
    dst->items = ta_items(t);       /* for each transaction */
    i = *dst->items +1;         /* get and check the first item */
    if (i <= 0) i = 0;          /* if sentinel, get the first list, */
    else dst->items++;          /* otherwise skip the item */
    tal       = lists +i;       /* sum transaction weights per item */
    tal->occ += dst->occ = ta_wgt(t);
    tal->wgt += dst->wgt = (double)dst->occ;
    dst->succ = tal->head;      /* add the new element to the */
    tal->head = dst++;          /* list for the first item */
  }                             /* and skip this item */
  r = rec_ins(relim, lists, k, (TID)(dst-elems));
  free(lists);                  /* execute recursive elimination */
  if (r >= 0)                   /* finally report the empty item set */
    r = isr_report(relim->report);
  return r;                     /* return the error status */
}  /* relim_ins() */

/*----------------------------------------------------------------------
  Recursive Elimination: Limited Item Insertions
----------------------------------------------------------------------*/

static int rec_lim (RELIM *relim, TZLIST *lists, ITEM k, TID n)
{                               /* --- recursive elimination */
  int    r;                     /* error status */
  ITEM   i;                     /* current item */
  TZLIST *proj  = NULL;         /* list(s) of (projected) database */
  TZLE   *elems = NULL;         /* to traverse the transaction lists */
  TZLIST *cur, *tal;            /* current/projected transaction list */
  TZLE   *src, *dst;            /* to traverse the transaction lists */
  double pen, wgt, w;           /* insertion penalty and tra. weight */
  double pex;                   /* minimum weight for perfect exts. */

  assert(relim && lists && (k > 0)); /* check the function arguments */
  #ifdef RELIM_ABORT            /* if to check for interrupt */
  if (sig_aborted()) return -1; /* if execution was aborted, */
  #endif                        /* abort the recursion */
  if ((k > 1)                   /* if there is more than one item and */
  &&  isr_xable(relim->report, 2)) {     /* another item can be added */
    proj = (TZLIST*)malloc((size_t)k *sizeof(TZLIST)
                          +(size_t)n *sizeof(TZLE));
    if (!proj) return -1;       /* allocate list and element arrays */
    elems = (TZLE*)(proj +k);   /* and organize the memory */
    memset(proj, 0, (size_t)k *sizeof(TZLIST));
  }                             /* initialize the projection header */
  pex = (relim->mode & REL_PERFECT)
      ? isr_wgt(relim->report) : INFINITY;
  for (r = 0; --k >= 0; ) {     /* traverse the transaction lists */
    if (proj) {                 /* clear list of empty transactions */
      proj->head = NULL; proj->wgt = 0; proj->occ = 0; }
    dst = elems;                /* init. the projected database */
    cur = lists +k+1;           /* get the current transaction list */
    pen = ib_getpen(relim->base, k);
    if (pen > 0) {              /* if insertion penalty is positive */
      for (i = k; i >= 0; i--){ /* traverse the preceding lists */
        for (src = lists[i].head; src; src = src->succ) {
          wgt = relim->tnorm(src->wgt, pen);
          if (wgt < relim->twgt)/* compute the penalized weight and */
            continue;           /* skip trans. with insuff. weight */
          cur->wgt  += w = wgt *(double)src->cnt;
          if (!dst) continue;   /* sum the transaction weight */
          tal        = proj +i; /* for current and projected list */
          tal->wgt  += w;       /* initialize a new list element */
          dst->occ   = 0;       /* with the penalized weight, but */
          dst->wgt   = wgt;     /* no occurrences (item is missing) */
          dst->cnt   = src->cnt;
          dst->items = src->items;
          dst->succ  = tal->head;
          tal->head  = dst++;   /* keep the item array unchanged */
        }                       /* and add the new list element to */
      }                         /* the corresp. transaction list */
    } /* The above loop inserts the current item into transactions  */
      /* that do not contain it. It is the only difference to the   */
      /* standard version of the algorithm apart from the fact that */
      /* empty transactions are processed in an additional list.    */
    i = -1;                     /* default: no recursion */
    if (cur->wgt >= pex) {      /* if item is a perfect extension, */
      isr_addpex(relim->report, k);      /* add it to the reporter */
      i = 0; }                  /* and clear the projection flag */
    else if ((cur->occ >= relim->supp)    /* if both support values */
    &&       (cur->wgt >= relim->isup)) { /* are large enough */
      r = isr_addwgt(relim->report, k, cur->occ, cur->wgt);
      if (r < 0) break;         /* add current item to the reporter */
      if (r > 0) {              /* if the item needs processing */
        if ((k > 0) && proj) {  /* if another item can be added */
          if (cur->head         /* if list has more than one element */
          &&  cur->head->succ   /* and there are few enough items, */
          &&  (k <= relim->sort))     /* sort the transaction list */
            cur->head = sort_wgt(cur->head);
          for (src = cur->head; src; src = src->succ) {
            i = *src->items +1; /* get and check the first item */
            if (i <= 0) i = 0;  /* get first list for sentinels */
            tal = proj +i;      /* sum the transaction weights */
            tal->occ  += dst->occ = src->occ;
            dst->cnt   = src->cnt;  /* copy and sum number of occs. */
            tal->wgt  += src->wgt *(double)src->cnt;
            dst->wgt   = src->wgt;  /* copy and sum trans. weight */
            dst->items = src->items +((i > 0) ? 1 : 0);
            dst->succ  = tal->head;
            tal->head  = dst++; /* add the new list element to */
          }                     /* the corresp. transaction list */
          r = rec_lim(relim, proj, k, (TID)(dst-elems));
          if (r < 0) break;     /* find frequent item sets */
        }                       /* recursively in the projection */
        r = isr_report(relim->report);
        if (r < 0) break;       /* report the current item set */
        isr_remove(relim->report, 1);
      }                         /* remove the current item */
    }                           /* from the item set reporter */
    if ((i < 0) && proj)        /* clear the projected database */
      memset(proj, 0, (size_t)k *sizeof(TZLIST));
    cur->wgt = 0; cur->occ = 0; /* clear the transaction weight */
    while (cur->head) {         /* while the list is not empty */
      src       = cur->head;    /* remove the first element */
      cur->head = src->succ;    /* from the transaction list */
      i = *src->items +1;       /* get and check the first item */
      if (i <= 0) i = 0;        /* if sentinel, get the first list, */
      else src->items++;        /* otherwise skip the item */
      tal       = lists +i;     /* sum the number of occurrences */
      tal->occ += src->occ;     /* and the transaction weights */
      tal->wgt += (double)src->cnt *src->wgt;
      src->succ = tal->head;    /* reassign the transactions */
      tal->head = src;          /* based on their first items */
    }                           /* and skip this first item */
  }
  if (proj) free(proj);         /* delete list and element arrays */
  return r;                     /* return the error status */
}  /* rec_lim() */

/*--------------------------------------------------------------------*/

int relim_lim (RELIM *relim)
{                               /* --- recursive elimination  */
  int    r;                     /* result of recursion */
  ITEM   i, k;                  /* loop variable, number of items */
  TID    n;                     /* number of transactions */
  TRACT  *t;                    /* to traverse the transactions */
  TZLIST *lists, *tal;          /* (array of) transaction list(s) */
  TZLE   *elems, *dst;          /* (array of) trans. list element(s) */

  assert(relim);                /* check the function arguments */
  if (tbg_wgt(relim->tabag) < relim->supp)
    return 0;                   /* check the total transaction weight */
  k = tbg_itemcnt(relim->tabag);/* get the number of items */
  if (k <= 0) return isr_report(relim->report);
  n = tbg_cnt(relim->tabag);    /* get the number of transactions */
  lists = (TZLIST*)malloc((size_t)(k+1) *sizeof(TZLIST)
                         +(size_t) n    *sizeof(TZLE));
  if (!lists) return -1;        /* allocate lists and element arrays */
  dst = elems = (TZLE*)(lists +k+1);     /* and initialize the lists */
  memset(lists, 0, (size_t)(k+1) *sizeof(TZLIST));
  while (--n >= 0) {            /* traverse the transactions */
    t = tbg_tract(relim->tabag, n); /* initialize a new list element */
    dst->items = ta_items(t);       /* for each transaction */
    i = *dst->items +1;         /* get the first item and */
    if (i <= 0) i = 0;          /* if sentinel, get the first list, */
    else dst->items++;          /* otherwise skip the item */
    tal       = lists +i;       /* sum the number of occurrences */
    tal->occ += dst->occ = dst->cnt = ta_wgt(t);
    dst->wgt  = 1.0;            /* sum transaction weights per item */
    tal->wgt += (double)dst->cnt;
    dst->succ = tal->head;      /* add the new element to the */
    tal->head = dst++;          /* list for the first item */
  }                             /* (or the list of empty trans.) */
  n = (TID)(dst -elems);        /* shrink memory to needed size */
  lists = (TZLIST*)realloc(lists, (size_t)(k+1) *sizeof(TZLIST)
                                 +(size_t) n    *sizeof(TZLE));
  r = rec_lim(relim, lists, k, n);
  free(lists);                  /* recursively find freq. item sets */
  if (r >= 0)                   /* finally report the empty item set */
    r = isr_report(relim->report);
  return r;                     /* return the error status */
}  /* relim_lim() */

/*----------------------------------------------------------------------
  RElim (generic)
----------------------------------------------------------------------*/

RELIM* relim_create (int target, double smin, double sins,
                     ITEM zmin, ITEM zmax, int tnorm, double twgt,
                     int eval, double thresh, int algo, int mode)
{                               /* create a relim miner */
  RELIM *relim;                 /* created relim miner */

  /* --- make parameters consistent --- */
  if      (target & REL_MAXIMAL) target = ISR_MAXIMAL;
  else if (target & REL_CLOSED)  target = ISR_CLOSED;
  else                           target = ISR_FREQUENT;
  if ((tnorm < 0) || (tnorm >= (int)(sizeof(tnorms)/sizeof(*tnorms))))
    tnorm = REL_MIN;            /* check and adapt the t-norm */

  /* --- create a relim miner --- */
  relim = (RELIM*)malloc(sizeof(RELIM));
  if (!relim) return NULL;      /* create a relim miner */
  relim->target = target;       /* and store all parameters */
  relim->smin   = smin;
  relim->supp   = 1;
  relim->sins   = sins;
  relim->isup   = DBL_MIN;
  relim->tnorm  = tnorms[tnorm];
  relim->twgt   = twgt;
  relim->zmin   = zmin;
  relim->zmax   = zmax;
  relim->eval   = eval;
  relim->thresh = thresh/100.0;
  relim->algo   = algo;
  relim->mode   = mode;
  relim->sort   = 32;
  relim->tabag  = NULL;
  relim->report = NULL;
  relim->mem    = NULL;
  relim->fim16  = NULL;
  relim->buf    = NULL;
  relim->base   = NULL;
  return relim;                 /* return created relim miner */
}  /* relim_create() */

/*--------------------------------------------------------------------*/

void relim_delete (RELIM *relim, int deldar)
{                               /* --- delete an relim miner */
  if (deldar) {                 /* if to delete data and reporter */
    if (relim->report) isr_delete(relim->report, 0);
    if (relim->tabag)  tbg_delete(relim->tabag,  1);
  }                             /* delete if existing */
  free(relim);                  /* delete the base structure */
}  /* relim_delete() */

/*--------------------------------------------------------------------*/

int relim_data (RELIM *relim, TABAG *tabag, int sort)
{                               /* --- prepare data for relim */
  ITEM    m;                    /* number of items */
  int     pack;                 /* number of items to pack */
  double  smin;                 /* absolute minimum support */
  double  sins;                 /* absolute minimum support */
  SUPP    w;                    /* total transaction weight */
  #ifndef QUIET                 /* if to print messages */
  TID     n;                    /* number of transactions */
  clock_t t;                    /* timer for measurements */
  #endif                        /* (only needed for messages) */

  assert(tabag);                /* check the function arguments */
  relim->tabag = tabag;         /* note the transaction bag and */
  relim->base  = tbg_base(tabag);   /* the underlying item base */
  pack = relim->mode & REL_FIM16;
  if (pack > 16) pack = 16;     /* get number of items to pack */
  w = tbg_wgt(tabag);           /* compute absolute minimum support */
  smin = ceilsupp((relim->smin < 0) ? -relim->smin
                : (relim->smin/100.0) *(double)w *(1-DBL_EPSILON));
  relim->supp = (SUPP)ceilsupp(smin);
  sins =          (relim->sins < 0) ? -relim->sins
                : (relim->sins/100.0) *(double)w *(1-DBL_EPSILON);
  relim->isup = (sins > 0) ? sins : DBL_MIN;

  /* --- sort and recode items --- */
  CLOCK(t);                     /* start timer, print log message */
  XMSG(stderr, "filtering, sorting and recoding items ... ");
  m = tbg_recode(tabag, relim->supp, -1, -1, -sort);
  if (m <  0) return E_NOMEM;   /* recode items and transactions */
  if (m <= 0) return E_NOITEMS; /* and check the number of items */
  XMSG(stderr, "[%"ITEM_FMT" item(s)]", m);
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));

  /* --- sort and reduce transactions --- */
  CLOCK(t);                     /* start timer, print log message */
  XMSG(stderr, "sorting and reducing transactions ... ");
  tbg_filter(tabag, (relim->twgt >= 0) ? 0 : relim->zmin, NULL, 0);
  tbg_itsort(tabag, -1, 0);     /* sort items in transactions and */
  tbg_sort  (tabag, -1, 0);     /* sort the trans. lexicographically */
  tbg_reduce(tabag, 0);         /* reduce transactions to unique ones */
  if ((relim->twgt < 0)         /* if insertions are not possible */
  &&  (pack        > 0))        /* and 16-items machine is used, */
    tbg_pack(tabag, pack);      /* pack the most frequent items */
  #ifndef QUIET                 /* if to print messages */
  n = tbg_cnt(tabag);           /* get the number of transactions */
  w = tbg_wgt(tabag);           /* and the new transaction weight */
  XMSG(stderr, "[%"TID_FMT, n); /* print number of transactions */
  if (w != (SUPP)n) { XMSG(stderr, "/%"SUPP_FMT, w); }
  XMSG(stderr, " transaction(s)] done [%.2fs].\n", SEC_SINCE(t));
  #endif
  return 0;                     /* return 'ok' */
}  /* relim_data() */

/*--------------------------------------------------------------------*/

int relim_report (RELIM *relim, ISREPORT *report)
{                               /* --- prepare reporter for relim */
  TID n;                        /* number of transactions */
  int mrep;                     /* mode for item set reporter */

  assert(report);               /* check the function arguments */
  relim->report = report;       /* note the item set reporter */

  /* --- get reporting mode --- */
  mrep = 0;                     /* init. the reporting mode */
  #ifdef USE_ZLIB               /* if optional output compression */
  if (relim->mode & REL_ZLIB)   /* if the compression flag is set, */
    mrep |= ISR_ZLIB;           /* transfer it to the report mode */
  #endif

  /* --- configure item set reporter --- */
  isr_setsupp(report, (RSUPP)relim->supp, RSUPP_MAX);
  isr_setsize(report, relim->zmin, relim->zmax);
  if (relim->eval == REL_LDRATIO)  /* set add. evaluation measure */
    isr_seteval(report, isr_logrto, NULL, +1, relim->thresh);
  n = (relim->mode & REL_PREFMT)/* get range of numbers to preformat */
    ? (TID)ib_maxfrq(tbg_base(relim->tabag)) : -1;
  if ((isr_prefmt(report, (TID)relim->supp, n)      != 0)
  ||  (isr_settarg(report, relim->target, mrep, -1) != 0))
    return E_NOMEM;             /* set pre-format and target type */
  return 0;                     /* return 'ok' */
}  /* relim_report() */

/*--------------------------------------------------------------------*/

int relim_mine (RELIM *relim, ITEM sort)
{                               /* --- RElim algorithm (generic) */
  int     r;                    /* result of function call */
  #ifndef QUIET                 /* if to print messages */
  clock_t t;                    /* timer for measurements */
  #endif                        /* (only needed for messages) */

  assert(relim);                /* check the function arguments */
  CLOCK(t);                     /* start timer, print log message */
  XMSG(stderr, "writing %s ... ", isr_name(relim->report));
  relim->sort = sort;           /* note the sorting threshold */
  if      (relim->twgt >  0)
    r = relim_lim (relim);      /* limited   item insertions */
  else if (relim->twgt >= 0)
    r = relim_ins (relim);      /* unlimited item insertions */
  else if (relim->algo == REL_TREE)
    r = relim_tree(relim);      /* tree-based algorithm */
  else if (relim->mode & REL_FIM16)
    r = relim_m16 (relim);      /* using a 16 items machine */
  else
    r = relim_base(relim);      /* standard list-based algorithm */
  if (r < 0) return E_NOMEM;    /* search for frequent item sets */
  XMSG(stderr, "[%"SIZE_FMT" set(s)]", isr_repcnt(relim->report));
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  return 0;                     /* return 'ok' */
}  /* relim_mine() */

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/
#ifdef RELIM_MAIN

static void help (void)
{                               /* --- print add. option information */
  #ifndef QUIET
  fprintf(stderr, "\n");        /* terminate startup message */
  printf("t-norms (triangular norms) for combining item penalties"
         " (option -N#)\n");
  printf("  m   minimum              T(a,b) = min(a,b)\n");
  printf("  n   nil-potent minimum   T(a,b) = min(a,b) "
                                    "if a+b > 1 else 0\n");
  printf("  p   product              T(a,b) = a*b\n");
  printf("  l   Lukasiewicz          T(a,b) = max(0,a+b-1)\n");
  printf("  h   Hamacher product     T(a,b) = 0 if a = b = 0 "
                                    "else a*b/(a+b-a*b)\n");
  printf("\n");
  printf("additional evaluation measures (option -e#)\n");
  printf("  x   no measure (default)\n");
  printf("  b   binary logarithm of support quotient\n");
  printf("\n");
  printf("relim algorithm variants (option -A#)\n");
  printf("  l   basic list-based algorithm (default)\n");
  printf("  t   tree-based algorithm\n");
  printf("\n");
  printf("information output format characters (option -v#)\n");
  printf("  %%%%  a percent sign\n");
  printf("  %%i  number of items (item set size)\n");
  printf("  %%a  absolute item set support\n");
  printf("  %%s  relative item set support as a fraction\n");
  printf("  %%S  relative item set support as a percentage\n");
  printf("  %%w  absolute support with insertions\n");
  printf("  %%r  relative support with insertions as a fraction\n");
  printf("  %%R  relative support with insertions as a percentage\n");
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

  assert(s && end);             /* check the function arguments */
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
  if (relim)  relim_delete(relim, 0); \
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
  int     target   = 's';       /* target type (closed/maximal) */
  ITEM    zmin     =  1;        /* minimum size of an item set */
  ITEM    zmax     = ITEM_MAX;  /* maximum size of an item set */
  double  supp     = 10;        /* minimum support of an item set */
  double  sins     = 10;        /* minimum support with insertions */
  int     tnorm    = 'p';       /* t-norm for combining item weights */
  double  twgt     = -1;        /* minimum transaction weight */
  int     eval     = 'x';       /* additional evaluation measure */
  double  thresh   = 10;        /* threshold for evaluation measure */
  int     sort     =  2;        /* flag for item sorting and recoding */
  int     algo     = 'l';       /* variant of RElim algorithm */
  int     mode     = REL_DEFAULT|REL_PREFMT;   /* search mode */
  int     pack     = 16;        /* number of bit-packed items */
  ITEM    slist    = 32;        /* transaction list sorting threshold */
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
    printf("         (s: frequent, c: closed, m: maximal item sets)\n");
    printf("-m#      minimum number of items per item set     "
                    "(default: %"ITEM_FMT")\n", zmin);
    printf("-n#      maximum number of items per item set     "
                    "(default: no limit)\n");
    printf("-s#      minimum support of an item set           "
                    "(default: %g%%)\n", supp);
    printf("         (positive: percentage, "
                     "negative: absolute number)\n");
    printf("-i#      minimum support with item insertions     "
                    "(default: %g%%)\n", sins);
    printf("         (only with item insertions, option -u)\n");
    printf("-T#      t-norm for combining item penalties      "
                    "(default: %c)\n",   tnorm);
    printf("-u#      minimum weight of a transaction          "
                    "(default: %g)\n",   twgt);
    printf("         (a value >= 0 selects item insertions)\n");
    printf("-e#      additional evaluation measure            "
                    "(default: none)\n");
    printf("-d#      threshold for add. evaluation measure    "
                    "(default: %g%%)\n", thresh);
    printf("-q#      sort items w.r.t. their frequency        "
                    "(default: %d)\n", sort);
    printf("         (1: ascending, -1: descending, 0: do not sort,\n"
           "          2: ascending, -2: descending w.r.t. "
                    "transaction size sum)\n");
    printf("-A#      variant of the relim algorithm to use    "
                    "(default: %c)\n", algo);
    printf("-x       do not prune with perfect extensions     "
                    "(default: prune)\n");
    printf("-l#      number of items for k-items machine      "
                    "(default: %d)\n", pack);
    printf("-y#      threshold for transaction list sorting   "
                    "(default: %"ITEM_FMT")\n", slist);
    printf("-F#:#..  support border for filtering item sets   "
                    "(default: none)\n");
    printf("         (list of minimum support values, "
                    "one per item set size,\n");
    printf("         starting at the minimum size, "
                    "as given with option -m#)\n");
    printf("-R#      read item selection/insertion penalties\n");
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
  /* free option characters: acjop [A-Z]\[CFNPRTZ] */

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
          case 's': supp   =       strtod(s, &s);    break;
          case 'i': sins   =       strtod(s, &s);    break;
          case 'T': tnorm  = (*s) ? *s++ : 'p';      break;
          case 'u': twgt   =       strtod(s, &s);    break;
          case 'e': eval   = (*s) ? *s++ : 0;        break;
          case 'd': thresh =       strtod(s, &s);    break;
          case 'q': sort   = (int) strtol(s, &s, 0); break;
          case 'A': algo   = (*s) ? *s++ : 0;        break;
          case 'x': mode  &= ~REL_PERFECT;           break;
          case 'l': pack   = (int) strtol(s, &s, 0); break;
          case 'y': slist  = (ITEM)strtol(s, &s, 0); break;
          case 'F': bdrcnt = getbdr(s, &s, &border); break;
          case 'R': optarg = &fn_sel;                break;
          case 'P': optarg = &fn_psp;                break;
          case 'Z': stats  = 1;                      break;
          case 'N': mode  |= REL_PREFMT;             break;
          case 'g': scan   = 1;                      break;
          #ifdef USE_ZLIB       /* if optional output compression */
          case 'z': mode  |= REL_ZLIB;               break;
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
  if (supp   > 100) error(E_SUPPORT, supp);
  if (sins   > 100) error(E_SUPPORT, sins);
  if (twgt   > 1)   error(E_WEIGHT,  twgt);
  if (bdrcnt < 0)   error(E_NOMEM);
  if ((!fn_inp || !*fn_inp) && (fn_sel && !*fn_sel))
    error(E_STDIN);             /* stdin must not be used twice */
  switch (target) {             /* check and translate target type */
    case 's': target = ISR_ALL;              break;
    case 'c': target = ISR_CLOSED;           break;
    case 'm': target = ISR_MAXIMAL;          break;
    default : error(E_TARGET, (char)target); break;
  }                             /* (get the target type code) */
  switch (tnorm) {              /* check and translate t-norm */
    case 'm': tnorm = REL_MIN;               break;
    case 'n': tnorm = REL_NILP;              break;
    case 'p': tnorm = REL_PROD;              break;
    case 'l': tnorm = REL_LUKA;              break;
    case 'h': tnorm = REL_HAMA;              break;
    default : error(E_TNORM, (char)tnorm);   break;
  }                             /* (get triangular norm) */
  switch (eval) {               /* check and translate measure */
    case 'x': eval = REL_NONE;               break;
    case 'b': eval = REL_LDRATIO;            break;
    default : error(E_MEASURE, (char)eval);  break;
  }                             /* (get evaluation measure code) */
  switch (algo) {               /* check and translate alg. variant */
    case 'l': algo = REL_BASIC;              break;
    case 't': algo = REL_TREE;               break;
    default : error(E_VARIANT, (char)algo);  break;
  }                             /* (get fpgrowth algorithm code) */
  mode = (mode & ~REL_FIM16)    /* add packed items to search mode */
       | ((pack <= 0) ? 0 : (pack < 16) ? pack : 16);
  if (slist < 0) slist = ITEM_MAX;
  if (info == dflt)             /* adapt the default info. format */
    info = (supp < 0) ? " (%a)" : " (%S)";
  MSG(stderr, "\n");            /* terminate the startup message */
  mode |= REL_VERBOSE|REL_NOCLEAN;

  /* --- read item selection/insertion penalties --- */
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
    m = (twgt >= 0)             /* depending on the target type */
      ? ib_readpen(ibase,tread) /* read the insertion penalties */
      : ib_readsel(ibase,tread);/* or the given item selection */
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
  relim = relim_create(target, supp, sins, zmin, zmax, tnorm, twgt,
                       eval, thresh, algo, mode);
  if (!relim) error(E_NOMEM);   /* create a split and merge miner */
  k = relim_data(relim, tabag, sort);
  if (k) error(k);              /* prepare data for split and merge */
  report = isr_create(ibase);   /* create an item set reporter */
  if (!report) error(E_NOMEM);  /* and configure it */
  k = relim_report(relim, report);
  if (k) error(k);              /* prepare reporter for relim */
  if (setbdr(report, w, zmin, &border, bdrcnt) != 0)
    error(E_NOMEM);             /* set the support border */
  if (fn_psp && (isr_addpsp(report, NULL) < 0))
    error(E_NOMEM);             /* set a pattern spectrum if req. */
  if (isr_setfmt(report, scan, hdr, sep, NULL, info) != 0)
    error(E_NOMEM);             /* set the output format strings */
  k = isr_open(report, NULL, fn_out);
  if (k) error(k, isr_name(report));
  if (isr_setup(report) < 0)     /* open the item set file and */
    error(E_NOMEM);              /* set up the item set reporter */
  k = relim_mine(relim, slist);  /* find frequent item sets */
  if (k) error(k);
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

#endif
