/*----------------------------------------------------------------------
  File    : fpgpsp.c
  Contents: generate or estimate a pattern spectrum (FP-growth)
  Author  : Christian Borgelt
  History : 2015.08.28 file created
            2016.11.20 fpgrowth miner object and interface introduced
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <assert.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <pthread.h>
#endif
#ifdef FPG_ABORT
#include "sigint.h"
#endif
#include "random.h"
#ifndef TA_SURR
#define TA_SURR
#endif
#ifdef FPGPSP_MAIN
#ifndef TA_READ
#define TA_READ
#endif
#endif
#include "tract.h"
#ifndef PSP_ESTIM
#define PSP_ESTIM
#endif
#ifdef FPGPSP_MAIN
#ifndef PSP_REPORT
#define PSP_REPORT
#endif
#endif
#ifndef ISR_PATSPEC
#define ISR_PATSPEC
#endif
#include "patspec.h"
#include "fpgpsp.h"
#ifdef FPGPSP_MAIN
#include "error.h"
#endif
#ifdef STORAGE
#include "storage.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define PRGNAME     "fpgpsp"
#define DESCRIPTION "generate or estimate a pattern spectrum " \
                    "(FP-growth)"
#define VERSION     "version 1.2 (2016.11.21)         " \
                    "(c) 2015-2016   Christian Borgelt"

/* --- thread definitions --- */
#ifdef _WIN32                   /* if Microsoft Windows system */
#define THREAD       HANDLE     /* threads identified by handles */
#define THREAD_OK    0          /* return value is DWORD */
#define WORKERDEF(n,p)  DWORD WINAPI n (LPVOID p)
#else                           /* if Linux/Unix system */
#define THREAD       pthread_t  /* use the POSIX thread type */
#define THREAD_OK    NULL       /* return value is void* */
#define WORKERDEF(n,p)  void*        n (void* p)
#endif                          /* definition of a worker function */

/* --- error codes --- */
/* error codes   0 to  -4 defined in tract.h */
#define E_STDIN      (-5)       /* double assignment of stdin */
#define E_OPTION     (-6)       /* unknown option */
#define E_OPTARG     (-7)       /* missing option argument */
#define E_ARGCNT     (-8)       /* too few/many arguments */
#define E_TARGET     (-9)       /* invalid target type */
#define E_SIZE      (-10)       /* invalid item set size */
#define E_SUPPORT   (-11)       /* invalid item set support */
#define E_VARIANT   (-12)       /* invalid algorithm variant */
#define E_SURRCNT   (-13)       /* invalid number of surrogates */
#define E_SURR      (-14)       /* invalid surrogate data gen. method */
/* error codes -15 to -26 defined in tract.h and train.h */
#define E_NOTABLE   (-27)       /* transaction bag is not a table */
#define E_ABORTED   (-28)       /* processing aborted by user */

#ifndef QUIET                   /* if not quiet version, */
#define MSG         fprintf     /* print messages */
#define XMSG        if (mode & FPG_VERBOSE) fprintf
#define CLOCK(t)    ((t) = clock())
#else                           /* if quiet version, */
#define MSG(...)    ((void)0)   /* suppress messages */
#define XMSG(...)   ((void)0)
#define CLOCK(t)    ((void)0)
#endif

#define SEC_SINCE(t)  ((double)(clock()-(t)) /(double)CLOCKS_PER_SEC)

/*----------------------------------------------------------------------
  Constants
----------------------------------------------------------------------*/
#if !defined QUIET && defined FPGPSP_MAIN
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
  /* E_SUPPORT -11 */  "invalid minimum support %"SUPP_FMT,
  /* E_VARIANT -12 */  "invalid fpgrowth variant '%c'",
  /* E_SURRCNT -13 */  "invalid number of surrogates %ld",
  /* E_SURR    -14 */  "invalid surrogate data generation method '%c'",
  /* E_NOITEMS -15 */  "no (frequent) items found",
  /*    -16 to -21 */  NULL, NULL, NULL, NULL, NULL, NULL,
  /*    -22 to -26 */  NULL, NULL, NULL, NULL, NULL,
  /* E_NOTABLE -27 */  "transactions are not table-derived",
  /* E_ABORTED -28 */  "processing aborted by user",
  /*           -28 */  "unknown error"
};
#endif

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- thread worker data --- */
  int       id;                 /* id of corresponding thread */
  FPGROWTH  *fpgrowth;          /* fpgrowth miner */
  TABAG     *tabag;             /* transaction bag to analyze */
  TABAG     *tasur;             /* buffer for surrogate data set */
  TBGSURRFN *surrfn;            /* surrogate data generator function */
  long      cnt;                /* number of surrogate data sets */
  RNG       *rng;               /* random number generator */
  ISREPORT  *report;            /* item set reporter (one per thread) */
  int       err;                /* error indicator */
  volatile long *done;          /* number of completed data sets */
  PRGREPFN  *repfn;             /* progress reporting function */
  void      *data;              /* progress reporting function data */
} WORKDATA;                     /* (thread worker data) */

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
static TBGSURRFN *sur_tab[] = { /* surrogate generation functions */
  /* IDENT   0 */  tbg_ident,   /* identity (keep original data) */
  /* RANDOM  1 */  tbg_random,  /* transaction randomization */
  /* SWAP    2 */  tbg_swap,    /* permutation by pair swaps */
  /* SHUFFLE 3 */  tbg_shuffle, /* shuffle table-derived data */
};

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
#ifdef FPGPSP_MAIN
#ifndef QUIET
static CCHAR    *prgname;       /* program name for error messages */
#endif
static TABREAD  *tread  = NULL; /* table/train reader */
static ITEMBASE *ibase  = NULL; /* item base */
static TABAG    *tabag  = NULL; /* transaction bag */
static PATSPEC  *psp    = NULL; /* pattern spectrum */
static TABWRITE *twrite = NULL; /* table writer */
static RSUPP    *border = NULL; /* decision border */
#endif

/*----------------------------------------------------------------------
  CPUinfo Function
----------------------------------------------------------------------*/

static int cpucnt (void)
{                               /* --- get the number of processors */
  #ifdef _WIN32                 /* if Microsoft Windows system */
  SYSTEM_INFO sysinfo;          /* system information structure */
  GetSystemInfo(&sysinfo);      /* get system information */
  return (int)sysinfo.dwNumberOfProcessors;
  #elif defined _SC_NPROCESSORS_ONLN
  return (int)sysconf(_SC_NPROCESSORS_ONLN);
  #else                         /* if no direct function available */
  static int cpuinfo[5] = { 0, 0, 0, 0, 0 };
  if (!cpuinfo[4]) {            /* if cpu info. not yet retrieved */
    __asm__ __volatile__ (      /* execute cpuid instruction */
      "cpuid": "=a" (cpuinfo[0]), "=b" (cpuinfo[1]),
               "=c" (cpuinfo[2]), "=d" (cpuinfo[3]) : "a" (1) );
    cpuinfo[4] = -1;            /* execute cpuid instruction */
  }                             /* and set initialization flag */
  return (cpuinfo[1] >> 16) & 0xff;  /* get the number of */
  #endif                        /* addressable processors */
}  /* cpucnt() */

/*----------------------------------------------------------------------
  Pattern Spectrum Generation
----------------------------------------------------------------------*/

static WORKERDEF(worker, p)
{                               /* --- worker function for a thread */
  WORKDATA *w = p;              /* type the argument pointer */
  long     i;                   /* loop variable for data sets */

  assert(p);                    /* check the function argument */
  for (i = 1; i <= w->cnt; i++){/* generate cnt surrogate data sets */
    w->tasur = w->surrfn(w->tabag, w->rng, w->tasur);
    #ifdef FPG_ABORT            /* if a signal handler is present */
    if (sig_aborted()) break;   /* check for an abort interrupt */
    #endif
    w->err |= fpg_data(w->fpgrowth, w->tasur, FPG_SURR, 0);
    if (w->err < 0)    break;   /* prepare the transactions */
    #ifdef FPG_ABORT            /* if a signal handler is present */
    if (sig_aborted()) break;   /* check for an abort interrupt */
    #endif
    w->err |= fpg_mine(w->fpgrowth, ITEM_MIN, 0);
    if (w->err < 0)    break;   /* execute the fpgrowth algorithm */
    #ifdef FPG_ABORT            /* if a signal handler is present */
    if (sig_aborted()) break;   /* check for an abort interrupt */
    #endif
    w->done[0]++;               /* count the surrogate data set */
    if (w->repfn) w->repfn(w->done[0], w->data);
  }                             /* report the progress */
  return THREAD_OK;             /* return a dummy result */
}  /* worker() */

/*--------------------------------------------------------------------*/

PATSPEC* fpg_genpsp (TABAG *tabag, int target, double supp,
                     ITEM zmin, ITEM zmax, int algo, int mode,
                     size_t cnt, int surr, long seed,
                     int cpus, PRGREPFN *rep, void* data)
{                               /* --- generate a pattern spectrum */
  PATSPEC   *psp = NULL;        /* created pattern spectrum */
  int       r;                  /* result of function call */
  TABAG     *tasur, *s;         /* surrogate data set */
  TBGSURRFN *surrfn;            /* surrogate data generation function */
  RNG       *rng;               /* random number generator */
  ISREPORT  *report;            /* item set reporter */
  FPGROWTH  *fpgrowth;          /* fpgrowth miner */
  THREAD    *threads;           /* thread handles */
  WORKDATA  *w;                 /* data for worker thread */
  long      c, x;               /* number of data sets per thread */
  long      i;                  /* loop variable for data sets */
  int       k, n;               /* loop variables for threads */
  #ifdef _WIN32                 /* if Microsoft Windows system */
  DWORD     thid;               /* dummy for storing the thread id */
  #endif                        /* (not really needed here) */
  volatile long done = 0;       /* number of completed surrogates */

  assert(tabag                  /* check the function arguments */
  &&    (algo >= FPG_SIMPLE) && (algo <= FPG_TOPDOWN));
  if (seed == 0) seed = (long)time(NULL);

  /* --- prepare data --- */
  fpgrowth = fpg_create(target, supp, 100.0, 100.0, zmin, zmax,
                        RE_NONE, FPG_NONE, 0.0, algo, mode);
  if (!fpgrowth) return NULL;   /* create an fpgrowth miner */
  r = (surr == FPG_SHUFFLE) ? FPG_NORECODE|FPG_NOSORT : 0;
  r = fpg_data(fpgrowth, tabag, r|FPG_NOPACK, +2);  /* prepare the */
  if (r) { fpg_delete(fpgrowth, 0); return NULL; }  /* transactions */

  /* --- generate pattern spectrum --- */
  if (cpus <= 0) cpus = cpucnt();
  if ((cpus > 1) && (cnt > 1)){ /* if to use multi-threading */
    threads = calloc((size_t)cpus, sizeof(THREAD));
    if (!threads) return NULL;  /* create array of thread handles */
    w = calloc((size_t)cpus, sizeof(WORKDATA));
    if (!w) { free(threads); return NULL; }
    c = ((long)cnt+cpus-1)/cpus;/* number of data sets per thread */
    for (n = 0; n < cpus; n++){ /* traverse the cpus/threads */
      x = (long)cnt -n*c;       /* get the number of data sets and */
      if (x <= 0) break;        /* check whether thread is needed */
      if (!fpgrowth) {          /* if no miner for this thread */
        fpgrowth = fpg_create(target, supp, 100.0, 100.0, zmin, zmax,
                              RE_NONE, FPG_NONE, 0.0, algo, mode);
      }                         /* create an fpgrowth miner */
      w[n].id       = n;        /* note the thread identifier */
      w[n].fpgrowth = fpgrowth; /* create an fpgrowth miner and */
      w[n].tabag    = tabag;    /* note transactions and a clone */
      w[n].tasur    = tbg_clone(tabag);
      w[n].surrfn   = sur_tab[surr];
      w[n].cnt      = (x < c) ? x : c;
      w[n].rng      = rng_create((unsigned int)(seed+n));
      w[n].report   = isr_create(tbg_base(tabag));
      w[n].err      = 0;        /* create random number generator */
      w[n].done     = &done;    /* and an item set reporter */
      w[n].repfn    = rep;
      w[n].data     = data;
      if (!w[n].fpgrowth || !w[n].tasur || !w[n].rng || !w[n].report) {
        w[n].err = -1; break; } /* check for successful creation */
      if ((fpg_data  (w[n].fpgrowth, w[n].tasur,
                      FPG_NORECODE|FPG_NOREDUCE, 0) != 0)
      ||  (fpg_report(w[n].fpgrowth, w[n].report) != 0)
      ||  (isr_addpsp(w[n].report, NULL) < 0)
      ||  (isr_setup (w[n].report) != 0)) {
        w[n].err = -1; break; } /* set up the item set reporter */
      #ifdef _WIN32             /* if Microsoft Windows system */
      threads[n] = CreateThread(NULL, 0, worker, w+n, 0, &thid);
      if (!threads[n]) { w[n].err = -1; break; }
      #else                     /* if Linux/Unix system */
      if (pthread_create(threads+n, NULL, worker, w+n) != 0) {
        w[n].err = -1; break; } /* create a thread for each data set */
      #endif                    /* to compute surrogates in parallel */
      fpgrowth = NULL;          /* fpgrowth miner has been used */
    }
    #ifdef _WIN32               /* if Microsoft Windows system */
    WaitForMultipleObjects((DWORD)n, threads, TRUE, INFINITE);
    for (k = n; --k >= 0; )     /* wait for threads to finish, */
      CloseHandle(threads[k]);  /* then close all thread handles */
    #else                       /* if Linux/Unix system */
    for (k = n; --k >= 0; )     /* wait for threads to finish */
      pthread_join(threads[k], NULL);
    #endif                      /* (join threads with this one) */
    for (k = r = 0; k < n; k++)
      r |= w[k].err;            /* join the error indicators */
    if (r >= 0) {               /* if processing was successful */
      psp = isr_rempsp(w[0].report, 0);
      for (k = 1; k < n; k++) { /* traverse the threads */
        r = psp_addpsp(psp, isr_getpsp(w[k].report));
        if (r < 0) break;       /* traverse and sum pattern spectrums */
      }                         /* in the one of the first thread */
    }
    for (k = n; --k >= 0; ) {   /* traverse the worker data */
      if (w[k].fpgrowth) fpg_delete(w[k].fpgrowth, 0);
      if (w[k].tasur)    tbg_delete(w[k].tasur, 0);
      if (w[k].rng)      rng_delete(w[k].rng);
      if (w[k].report)   isr_delete(w[k].report, 0);
    }                           /* delete the data structures */
    free(w);                    /* delete array of worker data */
    free(threads); }            /* and array of thread handles */
  else {                        /* if to use only one thread */
    report = isr_create(tbg_base(tabag));
    if (!report) return NULL;   /* create an item set reporter */
    if ((fpg_report(fpgrowth, report) < 0)
    ||  (isr_addpsp(report, NULL)     < 0)
    ||  (isr_setup (report) != 0)) {
      isr_delete(report, 0); fpg_delete(fpgrowth, 0); return NULL; }
    rng = rng_create((unsigned int)seed);
    if (!rng) {                 /* create a random number generator */
      isr_delete(report, 0); return NULL; }
    surrfn = sur_tab[surr];     /* get the surrogate data function */
    tasur  = NULL; r = 0;       /* init. surrogate and return code */
    for (i = 1; i <= (long)cnt; i++) {
      s = surrfn(tabag, rng, tasur);
      if (!s) { r = -1; break;} /* generate a surrogate data set */
      r = fpg_data(fpgrowth, tasur = s, FPG_SURR, 0);
      if (r < 0) break;         /* prepare the data set */
      r = fpg_mine(fpgrowth, 0, 0);
      if (r < 0) break;         /* execute the FP-growth algorithm */
      if (rep) rep(i, data);    /* report the progress and */
      #ifdef CCN_ABORT          /* if a signal handler is present */
      if (sig_aborted()) break; /* check for an interrupt */
      #endif
    }
    psp = isr_rempsp(report,0); /* get the created pattern spectrum */
    if (tasur) tbg_delete(tasur, 0);
    fpg_delete(fpgrowth, 0);    /* delete the data objects */
    rng_delete(rng);            /* (fpgrowth miner, */
    isr_delete(report, 0);      /* random number generator, */
  }                             /* and item set reporter) */

  /* --- clean up --- */
  #ifdef FPG_ABORT              /* if a signal handler is present */
  if (sig_aborted()) { sig_remove(); return NULL; }
  #endif
  if (r < 0) { if (psp) psp_delete(psp); return NULL; }
  #ifdef FPG_ABORT              /* if a signal handler is present */
  sig_remove();                 /* remove the signal handler */
  #endif
  return psp;                   /* return the created Java object */
}  /* fpg_genpsp() */

/*----------------------------------------------------------------------
  Pattern Spectrum Estimation
----------------------------------------------------------------------*/

PATSPEC* fpg_estpsp (TABAG *tabag, int target, double supp,
                     ITEM zmin, ITEM zmax, size_t equiv,
                     double alpha, size_t smpls, long seed)
{                               /* --- estimate a pattern spectrum */
  PATSPEC  *psp;                /* created pattern spectrum */
  FPGROWTH *fpgrowth;           /* fpgrowth miner */
  int       r;                  /* result of function call */

  assert(tabag);                /* check the function arguments */
  if (seed <= 0) seed = (long)time(NULL);
  rseed((unsigned)seed);        /* seed random number generator */
  fpgrowth = fpg_create(target, supp, 100.0, 100.0, zmin, zmax,
                        RE_NONE, FPG_NONE, 0.0, FPG_AUTO, FPG_DEFAULT);
  if (!fpgrowth) return NULL;   /* create an fpgrowth miner */
  r = fpg_data(fpgrowth, tabag, FPG_NOPACK, +2);
  fpg_delete(fpgrowth, 0);      /* prepare the transactions */
  if (r) return NULL;           /* and check for an error */
  supp = (supp < 0) ? -supp     /* compute absolute minimum support */
       : (supp/100.0) *(double)tbg_wgt(tabag) *(1-DBL_EPSILON);
  psp = psp_create(zmin, zmax, (SUPP)ceilsupp(supp), tbg_cnt(tabag));
  if (!psp) return NULL;        /* create a pattern spectrum */
  if (psp_tbgest(tabag, psp, (size_t)equiv, alpha, (size_t)smpls) == 0)
    return psp;                 /* estimate pattern spectrum */
  psp_delete(psp); return NULL; /* on failure delete pattern spectrum */
}  /* fpg_estpsp() */

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/
#ifdef FPGPSP_MAIN

static void help (void)
{                               /* --- print add. option information */
  #ifndef QUIET
  fprintf(stderr, "\n");        /* terminate startup message */
  printf("fpgrowth algorithm variants (option -A#)\n");
  printf("  s   simple  tree nodes with only successor and parent\n");
  printf("  c   complex tree nodes with children and siblings "
               "(default)\n");
  printf("  d   top-down processing on a single prefix tree\n");
  printf("  t   top-down processing of the prefix trees\n");
  printf("Variant 'd' does not support mining closed/maximal item ");
  printf("sets,\nvariant 't' does not support the use of a k-items ");
  printf("machine, and\nonly variant 'c' supports item reordering ");
  printf("w.r.t. conditional support,\nbut closed/maximal item sets ");
  printf("can only be mined without reordering.\nThese restrictions ");
  printf("may be removed in future versions of this program.\n");
  printf("surrogate data generation methods "
         "(option -g#, default: -ge)\n");
  printf("  e    estimate a pattern spectrum (no surrogates)\n");
  printf("  i    identity (keep original data)\n");
  printf("  r    random transaction generation\n");
  printf("  p    permutation by pair swaps\n");
  printf("  s    shuffle table-derived data\n");
#endif                        /* print help information */
  exit(0);                      /* abort the program */
}  /* help() */

/*--------------------------------------------------------------------*/

static void repfn (long cnt, void *data)
{                               /* --- progress reporting function */
  if ((cnt > *(long*)data) && ((cnt % 20) == 0))
    fprintf(stderr, "%10ld\b\b\b\b\b\b\b\b\b\b", *(long*)data = cnt);
}  /* repfn() */

/*--------------------------------------------------------------------*/

#ifndef NDEBUG                  /* if debug version */
  #undef  CLEANUP               /* clean up memory and close files */
  #define CLEANUP \
  if (twrite) twr_delete(twrite, 1); \
  if (psp)    psp_delete(psp);       \
  if (tabag)  tbg_delete(tabag,  0); \
  if (tread)  trd_delete(tread,  1); \
  if (ibase)  ib_delete (ibase);
#endif

GENERROR(error, exit)           /* generic error reporting function */

/*--------------------------------------------------------------------*/

int main (int argc, char *argv[])
{                               /* --- main function */
  int     i, k = 0;             /* loop variables, counters */
  char    *s;                   /* to traverse the options */
  CCHAR   **optarg = NULL;      /* option argument */
  CCHAR   *fn_inp  = NULL;      /* name of input  file */
  CCHAR   *fn_sel  = NULL;      /* name of item selection file */
  CCHAR   *fn_psp  = NULL;      /* name of pattern spectrum file */
  CCHAR   *recseps = NULL;      /* record  separators */
  CCHAR   *fldseps = NULL;      /* field   separators */
  CCHAR   *blanks  = NULL;      /* blank   characters */
  CCHAR   *comment = NULL;      /* comment characters */
  int     target   = 's';       /* target type (closed/maximal) */
  double  supp     = 10;        /* minimum support of an item set */
  SUPP    smin     = 1;         /* minimum support of an item set */
  ITEM    zmin     = 1;         /* minimum size of an item set */
  ITEM    zmax     = ITEM_MAX;  /* maximum size of an item set */
  int     algo     = 'c';       /* variant of fpgrowth algorithm */
  int     mode     = FPG_DEFAULT;  /* search mode (e.g. pruning) */
  int     pack     = 16;        /* number of bit-packed items */
  int     mtar     = 0;         /* mode for transaction reading */
  long    cnt      = 1000;      /* number of surrogate data sets */
  int     surr     = 'e';       /* surrogate data generation method */
  double  alpha    = 0.5;       /* probability dispersion factor */
  long    smpls    = 1000;      /* number of item set samples */
  long    seed     = 0;         /* seed for random numbers */
  int     cpus     = 0;         /* number of cpus to use */
  long    done     = 0;         /* number of completed data sets */
  int     decbdr   = 0;         /* flag for printing the border */
  ITEM    m;                    /* number of items/trains */
  TID     n;                    /* number of points/spikes */
  SUPP    w;                    /* total transaction weight */
  size_t  z;                    /* number of signatures */
  RSUPP   bmin;                 /* for printing the decision border */
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
    printf("         (s: frequent, c: closed, m: maximal item sets,\n");
    printf("          g: generators)\n");
    printf("-m#      minimum number of items per item set     "
                    "(default: %"ITEM_FMT")\n", zmin);
    printf("-n#      maximum number of items per item set     "
                    "(default: no limit)\n");
    printf("-s#      minimum support of an item set           "
                    "(default: %g)\n", supp);
    printf("         (positive: percentage, "
                     "negative: absolute number)\n");
    printf("-A#      variant of the fpgrowth algorithm to use "
                    "(default: %c)\n", algo);
    printf("-x       do not prune with perfect extensions     "
                    "(default: prune)\n");
    printf("-l#      number of items for k-items machine      "
                    "(default: %d)\n", pack);
    printf("         (only for variants s and d, "
                    "options -As or -Ad)\n");
    printf("-i       do not sort items w.r.t. cond. support   "
                    "(default: sort)\n");
    printf("         (only for algorithm variant c, option -Ac)\n");
    printf("-u       do not use head union tail (hut) pruning "
                    "(default: use hut)\n");
    printf("         (only for maximal item sets, option -tm)\n");
    printf("-g#      surrogate generation method              "
                    "(default: %c)\n",  surr);
    printf("-c#      number of surrogate data sets            "
                    "(default: %ld)\n", cnt);
    printf("-e#      probability dispersion factor            "
                    "(default: %g)\n",  alpha);
    printf("-z#      number of item set samples per size      "
                    "(default: %ld)\n", smpls);
    printf("-S#      seed for random numbers                  "
                    "(default: time)\n");
    printf("-Z#      number of cpus/processor cores to use    "
                    "(default: %d)\n", cpus);
    printf("         (a value <= 0 means all cpus "
                    "reported as available)\n");
    printf("-D       print induced decision border\n");
    printf("-R#      read an item selection from a file       "
                    "(default: use all items)\n");
    printf("-r#      record/transaction separators            "
                    "(default: \"\\n\")\n");
    printf("-f#      field /item        separators            "
                    "(default: \" \\t,\")\n");
    printf("-b#      blank   characters                       "
                    "(default: \" \\t\\r\")\n");
    printf("-C#      comment characters                       "
                    "(default: \"#\")\n");
    printf("-!       print additional option information\n");
    printf("infile   file to read trains from                 "
                    "[required]\n");
    printf("outfile  file to write pattern spectrum to        "
                    "[optional]\n");
    return 0;                   /* print a usage message */
  }                             /* and abort the program */
  #endif  /* #ifndef QUIET */

  /* --- evaluate arguments --- */
  seed = (long)time(NULL);      /* and get a default seed value */
  for (i = 1; i < argc; i++) {  /* traverse arguments */
    s = argv[i];                /* get option argument */
    if (optarg) { *optarg = s; optarg = NULL; continue; }
    if ((*s == '-') && *++s) {  /* -- if argument is an option */
      while (*s) {              /* traverse options */
        switch (*s++) {         /* evaluate switches */
          case '!': help();                          break;
          case 't': target = (*s) ? *s++ : 0;        break;
          case 'm': zmin   = (ITEM)strtol(s, &s, 0); break;
          case 'n': zmax   = (ITEM)strtol(s, &s, 0); break;
          case 's': supp   =       strtod(s, &s);    break;
          case 'A': algo   = (*s) ? *s++ : 0;        break;
          case 'x': mode  &= ~FPG_PERFECT;           break;
          case 'l': pack   = (int) strtol(s, &s, 0); break;
          case 'i': mode  &= ~FPG_REORDER;           break;
          case 'u': mode  &= ~FPG_TAIL;              break;
          case 'g': surr   = (*s) ? *s++ : 0;        break;
          case 'c': cnt    =       strtol(s, &s, 0); break;
          case 'e': alpha  =       strtod(s, &s);    break;
          case 'z': smpls  =       strtol(s, &s, 0); break;
          case 'S': seed   =       strtol(s, &s, 0); break;
          case 'Z': cpus   =  (int)strtol(s, &s, 0); break;
          case 'D': decbdr = -1;                     break;
          case 'R': optarg = &fn_sel;                break;
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
        case  1: fn_psp = s;      break;
        default: error(E_ARGCNT); break;
      }                         /* note filenames */
    }
  }
  if (optarg) error(E_OPTARG);  /* check (option) arguments */
  if (k < 2)  error(E_ARGCNT);  /* and number of arguments */
  if (zmin < 0)   error(E_SIZE, zmin); /* check the size limits */
  if (zmax < 0)   error(E_SIZE, zmax); /* and the minimum support */
  if (supp > 100) error(E_SUPPORT, supp);
  if (cnt  <= 0)  error(E_SURRCNT, cnt);
  switch (target) {             /* check and translate target type */
    case 's': target = ISR_ALL;              break;
    case 'c': target = ISR_CLOSED;           break;
    case 'm': target = ISR_MAXIMAL;          break;
    case 'g': target = ISR_GENERAS;          break;
    default : error(E_TARGET, (char)target); break;
  }                             /* (get target type code) */
  switch (algo) {               /* check and translate alg. variant */
    case 's': algo = FPG_SIMPLE;             break;
    case 'c': algo = FPG_COMPLEX;            break;
    case 'd': algo = FPG_SINGLE;             break;
    case 't': algo = FPG_TOPDOWN;            break;
    default : error(E_VARIANT, (char)algo);  break;
  }                             /* (get fpgrowth algorithm code) */
  mode = (mode & ~FPG_FIM16)    /* add packed items to search mode */
       | ((pack <= 0) ? 0 : (pack < 16) ? pack : 16);
  switch (surr) {               /* check and translate surrogate */
    case 'e': surr = -1;                     break;
    case 'i': surr = FPG_IDENTITY;           break;
    case 'r': surr = FPG_RANDOM;             break;
    case 'p': surr = FPG_SWAP;               break;
    case 's': surr = FPG_SHUFFLE;            break;
    default : error(E_SURR, (char)surr);     break;
  }                             /* (get surrogate method code) */
  if (surr == FPG_IDENTITY) cnt = 1;
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
  if ((surr == FPG_SHUFFLE) && !tbg_istab(tabag))
    error(E_NOTABLE);           /* check for tabular data */
  MSG(stderr, "\n");            /* terminate the log message */

  /* --- estimate/generate pattern spectrum --- */
  CLOCK(t);                     /* start timer, print log message */
  if (surr < 0) {               /* if to estimate a pattern spectrum */
    MSG(stderr, "estimating pattern spectrum ... ");
    psp = fpg_estpsp(tabag, target, supp, zmin, zmax,
                     (size_t)cnt, alpha, (size_t)smpls, seed); }
  else {                        /* if to generate a pattern spectrum */
    MSG(stderr, "generating pattern spectrum ... ");
    psp = fpg_genpsp(tabag, target, supp, zmin, zmax,
                     algo, mode, (size_t)cnt, surr, seed,
                     cpus, repfn, &done);
  }                             /* (create a pattern spectrum) */
  #ifdef FPG_ABORT              /* if a signal handler is present */
  if (!psp) error((sig_aborted()) ? E_ABORTED : E_NOMEM);
  #else                         /* check for abort or error */
  if (!psp) error(E_NOMEM);     /* if no signal handler is present */
  #endif                        /* it must be out of memory */
  z = psp_sigcnt(psp);          /* get the number of signatures */
  MSG(stderr, "[%"SIZE_FMT" signature(s)]", z);
  MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));

  /* --- write pattern spectrum --- */
  CLOCK(t);                     /* start timer, print log message */
  twrite = twr_create();        /* create a table writer and */
  if (!twrite) error(E_NOMEM);  /* open the output file */
  if (twr_open(twrite, NULL, fn_psp) != 0)
    error(E_FOPEN,  twr_name(twrite));
  MSG(stderr, "writing %s ... ", twr_name(twrite));
  if (psp_report(psp, twrite, 1.0/(double)cnt) != 0)
    error(E_FWRITE, twr_name(twrite));
  twr_delete(twrite, 1);        /* write the pattern spectrum */
  twrite = NULL;                /* and delete the table writer */
  MSG(stderr, "[%"SIZE_FMT" signature(s)]", z);
  MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));

  /* --- print decision border --- */
  if (decbdr) {                 /* if to print the decision border */
    for (zmax = psp_max(psp); zmax >= zmin; zmax--)
      if (psp_max4sz(psp, zmax) >= psp_min4sz(psp, zmin))
        break;                  /* find the largest pattern size */
    bmin   = (RSUPP)smin;       /* get the minimum support */
    border = (RSUPP*)malloc(((size_t)zmax+1) *sizeof(RSUPP));
    if (!border) error(E_NOMEM);/* create a decision border */
    for (m = zmax+1; --m >= zmin; ) { /* traverse the pattern sizes */
      border[m] = psp_max4sz(psp, m) +SUPP_EPS;
      if (bmin > border[m]) border[m] = bmin;
      bmin = border[m];         /* get min. support for significance */
    }                           /* and ensure a monotone border */
    for (m = zmin; m < zmax; m++) printf("%"RSUPP_FMT":", border[m]);
    if (zmax > zmin) { printf("%"RSUPP_FMT"\n", border[zmax]); }
  }                             /* print the decision border */

  /* --- clean up --- */
  CLEANUP;                      /* clean up memory and close files */
  SHOWMEM;                      /* show (final) memory usage */
  return 0;                     /* return 'ok' */
}  /* main() */

#endif  /* #ifdef CCNPSP_MAIN */
