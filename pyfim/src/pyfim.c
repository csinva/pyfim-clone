/*----------------------------------------------------------------------
  File    : pyfim.c
  Contents: Frequent Item set Mining for Python
  Author  : Christian Borgelt
  History : 2011.07.13 file created
            2011.07.19 interface for apriori implementation added
            2011.07.22 adapted to new module ruleval (rule evaluation)
            2011.07.25 adapted to modified apriori() interface
            2011.07.28 translation of test statistic strings etc. added
            2011.08.03 subset filtering options added to apriacc
            2012.04.17 functions py_eclat() and py_fpgrowth() added
            2012.08.02 missing Py_DECREFs with Append/BuildValue added
            2013.01.16 bug in function apriori() fixed (isr_create())
            2013.01.24 transactions generalized to any sequence type
            2013.02.10 module initialization adapted to Python 3
            2013.02.11 items generalized to any hashable object
            2013.02.25 transactions generalized to any iterable object
            2013.03.06 dictionary allowed for transaction database
            2013.04.01 adapted to type changes in module tract
            2013.05.08 memory leak in tbg_fromPyObj() fixed (dict)
            2013.10.18 optional pattern spectrum collection added
            2013.10.31 carpenter and supporting functions added
            2013.11.21 sam and relim and supporting functions added
            2013.12.07 returned pattern spectrum as list or dictionary
            2014.05.06 pattern spectrum generation and estimation added
            2014.05.12 support border added as fim function argument
            2014.05.15 item set evaluations "cprob" and "import" added
            2014.05.27 bugs concerning new argument 'border' fixed
            2014.06.12 bug in function py_patspec() fixed (NULL init.)
            2014.08.01 minimum improvement of evaluation measure removed
            2014.08.24 adapted to modified item set reporter interface
            2014.08.25 ista algorithm and supporting functions added
            2014.08.28 adapted to new FIM function interfaces
            2014.09.19 function py_arules() and mode parameters added
            2014.10.02 rules as target added to apriori/eclat/fpgrowth
            2014.10.09 functions repinit() and repterm() added
            2014.10.15 bug in function fim() fixed (call to eclat fn.)
            2014.10.17 bug in function py_patspec() fixed (~FPG_FIM16)
            2014.10.24 changed from LGPL license to MIT license
            2015.02.27 item appearances param. added to rule functions
            2015.08.07 ensured signal handler removal on error
            2015.08.13 function py_patspec() renamed to py_genpsp()
            2015.08.15 function py_psp2bdr() added (decision border)
            2015.08.17 function py_patred() added (pat. set reduction)
            2015.08.19 adapted to generic pattern set reduction
            2015.09.04 reporting of associated values made more flexible
            2016.03.30 some cleanup of option strings (app, eval etc.)
            2016.05.02 excess Py_DECREF() in ib_appPyObj() deleted
            2016.10.06 optional floating point support made possible
            2016.10.07 conditional help compilation for float support
            2016.11.05 adapted to modified apriori interface
            2016.11.11 adapted to modified eclat interface
            2016.11.15 adapted to modified accretion interface
            2016.11.21 adapted to modified fpgrowth interface
            2017.03.24 adapted to modified carpenter/ista interfaces
----------------------------------------------------------------------*/
#include <assert.h>
#include <float.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <pthread.h>
#endif
#if defined _WIN32 && !defined HAVE_ROUND
#define HAVE_ROUND 1
#endif
#include <Python.h>
#ifndef TATREEFN
#define TATREEFN
#endif
#ifndef TA_SURR
#define TA_SURR
#endif
#ifndef PSP_ESTIM
#define PSP_ESTIM
#endif
#ifndef ISR_PATSPEC
#define ISR_PATSPEC
#endif
#ifndef ISR_CLOMAX
#define ISR_CLOMAX
#endif
#include "sigint.h"
#include "report.h"
#include "apriori.h"
#include "eclat.h"
#include "fpgrowth.h"
#include "sam.h"
#include "relim.h"
#include "carpenter.h"
#include "ista.h"
#include "accretion.h"
#include "fpgpsp.h"
#include "patred.h"

/*--------------------------------------------------------------------*/

#define int         1           /* to check definition of RSUPP */
#define long        2           /* for certain types */
#define ptrdiff_t   3
#define double      4

#if SUPP==double
#define FPSUPP      1           /* if support type is floating point */
#else
#define FPSUPP      0           /* if support type is integer */
#endif

#undef int                      /* remove preprocessor definitions */
#undef long                     /* needed for the type checking */
#undef ptrdiff_t
#undef double

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#if PY_MAJOR_VERSION >= 3
#define PyInt_Check     PyLong_Check
#define PyInt_AsLong    PyLong_AsLong
#define PyInt_FromLong  PyLong_FromLong
#else
#define Py_hash_t       long    /* type was introduced with Python 3 */
#endif

#if FPSUPP                      /* if floating point support */
#define PyObj_FromSupp(s)   PyFloat_FromDouble((double)s)
#else                           /* create a floating point object */
#define PyObj_FromSupp(s)   PyInt_FromLong((long)s)
#endif                          /* create an integer object */

/* --- error handling --- */
#define ERR_RTN()       error(0, NULL)
#define ERR_VALUE(s)    error(PyExc_ValueError,   s)
#define ERR_TYPE(s)     error(PyExc_TypeError,    s)
#define ERR_MEM()       error(PyExc_MemoryError,  "not enough memory")
#define ERR_ABORT()     error(PyExc_RuntimeError, "user abort")

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- item set report data --- */
  PyObject *res;                /* constructed result */
  int      mode;                /* value reporting mode */
  int      cnt;                 /* number of value indicators */
  CCHAR    *rep;                /* indicators of values to report */
  int      err;                 /* error flag */
} REPDATA;                      /* (item set report data) */

/*----------------------------------------------------------------------
  Item Functions
----------------------------------------------------------------------*/

static size_t hashitem (const void *a, int type)
{                               /* --- compute hash code */
  return (size_t)PyObject_Hash(*(PyObject**)a);
}  /* hashitem() */

/*--------------------------------------------------------------------*/

static int cmpitems (const void *a, const void *b, void *data)
{                               /* --- compare two Python objects */
  return PyObject_RichCompareBool(*(PyObject**)a,*(PyObject**)b,Py_NE);
}  /* cmpitems() */             /* return 0 if objects are equal */

/*----------------------------------------------------------------------
  Parameter Functions
----------------------------------------------------------------------*/

static int get_app (const char *s)
{                               /* --- get item appearance indicator */
  assert(s);                    /* check the function argument */
  if (s[0] && !s[1]) {          /* translate single character */
    switch (s[0]) {             /* evaluate single character */
      case 'n': s = "-"; break;
      case 'i': s = "a"; break;
      case 'b': s = "a"; break;
      case 'o': s = "c"; break;
      case 'h': s = "c"; break;
    } }
  else if (s[0] && s[1]) {      /* evaluate textual identifier */
    if      (strcmp(s, "none")       == 0) s = "-";
    else if (strcmp(s, "neither")    == 0) s = "-";
    else if (strcmp(s, "ign")        == 0) s = "-";
    else if (strcmp(s, "ignore")     == 0) s = "-";
    else if (strcmp(s, "in")         == 0) s = "a";
    else if (strcmp(s, "inp")        == 0) s = "a";
    else if (strcmp(s, "input")      == 0) s = "a";
    else if (strcmp(s, "out")        == 0) s = "c";
    else if (strcmp(s, "output")     == 0) s = "c";
    else if (strcmp(s, "ante")       == 0) s = "a";
    else if (strcmp(s, "antecedent") == 0) s = "a";
    else if (strcmp(s, "cons")       == 0) s = "c";
    else if (strcmp(s, "consequent") == 0) s = "c";
    else if (strcmp(s, "body")       == 0) s = "a";
    else if (strcmp(s, "head")       == 0) s = "c";
    else if (strcmp(s, "io")         == 0) s = "x";
    else if (strcmp(s, "i&o")        == 0) s = "x";
    else if (strcmp(s, "o&i")        == 0) s = "x";
    else if (strcmp(s, "inout")      == 0) s = "x";
    else if (strcmp(s, "in&out")     == 0) s = "x";
    else if (strcmp(s, "ac")         == 0) s = "x";
    else if (strcmp(s, "a&c")        == 0) s = "x";
    else if (strcmp(s, "c&a")        == 0) s = "x";
    else if (strcmp(s, "canda")      == 0) s = "x";
    else if (strcmp(s, "bh")         == 0) s = "x";
    else if (strcmp(s, "b&h")        == 0) s = "x";
    else if (strcmp(s, "h&b")        == 0) s = "x";
    else if (strcmp(s, "both")       == 0) s = "x";
  }
  if (s[0] && !s[1]) {          /* check for a valid string */
    switch (s[0]) {             /* evaluate the appearance code */
      case '-': return APP_NONE;
      case 'a': return APP_BODY;
      case 'c': return APP_HEAD;
      case 'x': return APP_BOTH;
    }
  }
  PyErr_SetString(PyExc_ValueError,"invalid item appearance indicator");
  return -1;                    /* return an error code */
}  /* get_app() */

/*--------------------------------------------------------------------*/

static int get_target (const char *s, const char *targets)
{                               /* --- get target */
  assert(s);                    /* check the function argument */
  if (s[0] && s[1]) {           /* evaluate textual identifier */
    if      (strcmp(s, "set")        == 0) s = "s";
    else if (strcmp(s, "sets")       == 0) s = "s";
    else if (strcmp(s, "all")        == 0) s = "s";
    else if (strcmp(s, "allset")     == 0) s = "s";
    else if (strcmp(s, "allsets")    == 0) s = "s";
    else if (strcmp(s, "frq")        == 0) s = "s";
    else if (strcmp(s, "freq")       == 0) s = "s";
    else if (strcmp(s, "frequent")   == 0) s = "s";
    else if (strcmp(s, "frqset")     == 0) s = "s";
    else if (strcmp(s, "frqsets")    == 0) s = "s";
    else if (strcmp(s, "freqset")    == 0) s = "s";
    else if (strcmp(s, "freqsets")   == 0) s = "s";
    else if (strcmp(s, "cls")        == 0) s = "c";
    else if (strcmp(s, "clsd")       == 0) s = "c";
    else if (strcmp(s, "closed")     == 0) s = "c";
    else if (strcmp(s, "max")        == 0) s = "m";
    else if (strcmp(s, "maxi")       == 0) s = "m";
    else if (strcmp(s, "maximal")    == 0) s = "m";
    else if (strcmp(s, "gen")        == 0) s = "g";
    else if (strcmp(s, "gens")       == 0) s = "g";
    else if (strcmp(s, "generas")    == 0) s = "g";
    else if (strcmp(s, "generators") == 0) s = "g";
    else if (strcmp(s, "rule")       == 0) s = "r";
    else if (strcmp(s, "rules")      == 0) s = "r";
    else if (strcmp(s, "arule")      == 0) s = "r";
    else if (strcmp(s, "arules")     == 0) s = "r";
  }
  if (s[0] && !s[1]             /* check for a valid string */
  && (strchr(targets, s[0]) != NULL)) {
    switch (s[0]) {             /* evaluate the target code */
      case 'a': return ISR_ALL;
      case 's': return ISR_SETS;
      case 'f': return ISR_FREQUENT;
      case 'c': return ISR_CLOSED;
      case 'm': return ISR_MAXIMAL;
      case 'g': return ISR_GENERAS;
      case 'r': return ISR_RULES;
    }
  }
  PyErr_SetString(PyExc_ValueError, "invalid target type");
  return -1;                    /* return an error code */
}  /* get_target() */

/*--------------------------------------------------------------------*/

static int get_stat (const char *s)
{                               /* --- get statistic code */
  assert(s);                    /* check the function argument */
  if (s[0] && s[1]) {           /* evaluate textual identifier */
    if      (strcmp(s, "none")      == 0) s = "x";
    else if (strcmp(s, "X2")        == 0) s = "p";
    else if (strcmp(s, "chi2")      == 0) s = "p";
    else if (strcmp(s, "X2pval")    == 0) s = "p";
    else if (strcmp(s, "chi2pval")  == 0) s = "p";
    else if (strcmp(s, "yates")     == 0) s = "t";
    else if (strcmp(s, "yatespval") == 0) s = "t";
    else if (strcmp(s, "info")      == 0) s = "g";
    else if (strcmp(s, "infopval")  == 0) s = "g";
    else if (strcmp(s, "fetprob")   == 0) s = "f";
    else if (strcmp(s, "fetchi2")   == 0) s = "h";
    else if (strcmp(s, "fetX2")     == 0) s = "h";
    else if (strcmp(s, "fetinfo")   == 0) s = "m";
    else if (strcmp(s, "fetsupp")   == 0) s = "s";
  }
  if (s[0] && !s[1]) {          /* check for a valid string */
    switch (s[0]) {             /* evaluate the statistic code */
      case 'x': return RE_NONE;
      case 'c': return RE_CHI2PVAL;
      case 'p': return RE_CHI2PVAL;
      case 'n': return RE_CHI2PVAL;
      case 'y': return RE_YATESPVAL;
      case 't': return RE_YATESPVAL;
      case 'i': return RE_INFOPVAL;
      case 'g': return RE_INFOPVAL;
      case 'f': return RE_FETPROB;
      case 'h': return RE_FETCHI2;
      case 'm': return RE_FETINFO;
      case 's': return RE_FETSUPP;
    }
  }
  PyErr_SetString(PyExc_ValueError, "invalid statistic");
  return -1;                    /* return an error code */
}  /* get_stat() */

/*--------------------------------------------------------------------*/

static int get_eval (const char *s)
{                               /* --- get evaluation measure code */
  assert(s);                    /* check the function argument */
  if (s[0] && s[1]) {
    if (strcmp(s, "none")    == 0) return 'x';
    if (strcmp(s, "ldratio") == 0) return 'b';
  }
  if (strchr("xb", s[0]) != NULL) return s[0];
  PyErr_SetString(PyExc_ValueError, "invalid evaluation measure");
  return -1;                    /* return an error code */
}  /* get_eval() */

/*--------------------------------------------------------------------*/

static int get_evalx (const char *s)
{                               /* --- get evaluation measure code */
  assert(s);                    /* check the function argument */
  if (s[0] && s[1]) {           /* evaluate textual identifier */
    if (strcmp(s, "none")       == 0) s = "x";
    if (strcmp(s, "supp")       == 0) s = "o";
    if (strcmp(s, "support")    == 0) s = "o";
    if (strcmp(s, "conf")       == 0) s = "c";
    if (strcmp(s, "confidence") == 0) s = "c";
    if (strcmp(s, "confdiff")   == 0) s = "d";
    if (strcmp(s, "lift")       == 0) s = "l";
    if (strcmp(s, "liftdiff")   == 0) s = "a";
    if (strcmp(s, "liftquot")   == 0) s = "q";
    if (strcmp(s, "cvct")       == 0) s = "v";
    if (strcmp(s, "conviction") == 0) s = "v";
    if (strcmp(s, "cvctdiff")   == 0) s = "e";
    if (strcmp(s, "cvctquot")   == 0) s = "r";
    if (strcmp(s, "cprob")      == 0) s = "k";
    if (strcmp(s, "import")     == 0) s = "j";
    if (strcmp(s, "importance") == 0) s = "j";
    if (strcmp(s, "cert")       == 0) s = "z";
    if (strcmp(s, "chi2")       == 0) s = "n";
    if (strcmp(s, "X2")         == 0) s = "n";
    if (strcmp(s, "chi2pval")   == 0) s = "p";
    if (strcmp(s, "X2pval")     == 0) s = "p";
    if (strcmp(s, "yates")      == 0) s = "y";
    if (strcmp(s, "yatespval")  == 0) s = "t";
    if (strcmp(s, "info")       == 0) s = "i";
    if (strcmp(s, "infopval")   == 0) s = "g";
    if (strcmp(s, "gpval")      == 0) s = "g";
    if (strcmp(s, "fetprob")    == 0) s = "f";
    if (strcmp(s, "fetchi2")    == 0) s = "h";
    if (strcmp(s, "fetX2")      == 0) s = "h";
    if (strcmp(s, "fetinfo")    == 0) s = "m";
    if (strcmp(s, "fetsupp")    == 0) s = "s";
    if (strcmp(s, "ldratio")    == 0) s = "b";
  }
  if (s[0] && !s[1]) {          /* check for a valid string */
    switch (s[0]) {             /* evaluate the measure code */
      case 'x': return RE_NONE;
      case 'o': return RE_SUPP;
      case 'c': return RE_CONF;
      case 'd': return RE_CONFDIFF;
      case 'l': return RE_LIFT;
      case 'a': return RE_LIFTDIFF;
      case 'q': return RE_LIFTQUOT;
      case 'v': return RE_CVCT;
      case 'e': return RE_CVCTDIFF;
      case 'r': return RE_CVCTQUOT;
      case 'k': return RE_CPROB;
      case 'j': return RE_IMPORT;
      case 'z': return RE_CERT;
      case 'n': return RE_CHI2;
      case 'p': return RE_CHI2PVAL;
      case 'y': return RE_YATES;
      case 't': return RE_YATESPVAL;
      case 'i': return RE_INFO;
      case 'g': return RE_INFOPVAL;
      case 'f': return RE_FETPROB;
      case 'h': return RE_FETCHI2;
      case 'm': return RE_FETINFO;
      case 's': return RE_FETSUPP;
      case 'b': return RE_FNCNT;
    }
  }
  PyErr_SetString(PyExc_ValueError, "invalid evaluation measure");
  return -1;                    /* return an error code */
}  /* get_evalx() */

/*--------------------------------------------------------------------*/

static int get_agg (const char *s)
{                               /* --- get aggregation mode */
  assert(s);                    /* check the function argument */
  if (s[0] && s[1]) {           /* evaluate textual identifier */
    if      (strcmp(s, "none")    == 0) s = "x";
    else if (strcmp(s, "min")     == 0) s = "m";
    else if (strcmp(s, "minimum") == 0) s = "m";
    else if (strcmp(s, "max")     == 0) s = "n";
    else if (strcmp(s, "maximum") == 0) s = "n";
    else if (strcmp(s, "avg")     == 0) s = "a";
    else if (strcmp(s, "average") == 0) s = "a";
  }
  if (s[0] && !s[1]) {          /* check for avalid string */
    switch (s[0]) {             /* evaluate the aggregation code */
      case 'x': return IST_NONE;
      case 'm': return IST_MIN;
      case 'n': return IST_MAX;
      case 'a': return IST_AVG;
    }
  }
  PyErr_SetString(PyExc_ValueError, "invalid aggregation mode");
  return -1;                    /* return an error code */
}  /* get_agg() */

/*--------------------------------------------------------------------*/

static int get_surr (const char *s)
{                               /* --- get surrogate function code */
  assert(s);                    /* check the function argument */
  if (s[0] && s[1]) {           /* evaluate textual identifier */
    if      (strcmp(s, "ident")     == 0) s = "i";
    else if (strcmp(s, "identity")  == 0) s = "i";
    else if (strcmp(s, "random")    == 0) s = "r";
    else if (strcmp(s, "randomize") == 0) s = "r";
    else if (strcmp(s, "swap")      == 0) s = "p";
    else if (strcmp(s, "perm")      == 0) s = "p";
    else if (strcmp(s, "permute")   == 0) s = "p";
    else if (strcmp(s, "shuffle")   == 0) s = "s";
  }
  if (s[0] && !s[1]) {          /* check for avalid string */
    switch (s[0]) {             /* evaluate the surrogate method code */
      case 'i': return FPG_IDENTITY;
      case 'r': return FPG_RANDOM;
      case 'p': return FPG_SWAP;
      case 'w': return FPG_SWAP;
      case 's': return FPG_SHUFFLE;
    }
  }
  PyErr_SetString(PyExc_ValueError,
                  "invalid surrogate generation method");
  return -1;                    /* return an error code */
}  /* get_surr() */

/*--------------------------------------------------------------------*/

static int get_red (const char *s)
{                               /* --- get random function code */
  assert(s);                    /* check the function argument */
  if (s[0] && s[1]) {           /* evaluate textual identifier */
    if      (strcmp(s, "none")        == 0) s = "x";
    else if (strcmp(s, "coins")       == 0) s = "c";
    else if (strcmp(s, "coins0")      == 0) s = "c";
    else if (strcmp(s, "coins1")      == 0) s = "C";
    else if (strcmp(s, "coins+1")     == 0) s = "C";
    else if (strcmp(s, "items")       == 0) s = "i";
    else if (strcmp(s, "items2")      == 0) s = "i";
    else if (strcmp(s, "neurons")     == 0) s = "i";
    else if (strcmp(s, "cover")       == 0) s = "s";
    else if (strcmp(s, "cover0")      == 0) s = "s";
    else if (strcmp(s, "covered")     == 0) s = "s";
    else if (strcmp(s, "covered0")    == 0) s = "s";
    else if (strcmp(s, "cover1")      == 0) s = "S";
    else if (strcmp(s, "covered1")    == 0) s = "S";
    else if (strcmp(s, "leni")        == 0) s = "l";
    else if (strcmp(s, "leni0")       == 0) s = "l";
    else if (strcmp(s, "lenient")     == 0) s = "l";
    else if (strcmp(s, "lenient0")    == 0) s = "l";
    else if (strcmp(s, "leni1")       == 0) s = "L";
    else if (strcmp(s, "lenient1")    == 0) s = "L";
    else if (strcmp(s, "strict")      == 0) s = "t";
    else if (strcmp(s, "strict0")     == 0) s = "t";
    else if (strcmp(s, "strict1")     == 0) s = "T";
  }
  if (s[0] && !s[1]) {          /* translate surrogate method string */
    switch (s[0]) {             /* evaluate the surrogate method code */
      case 'x': return PSR_NONE;
      case 'c': return PSR_COINS0;
      case 'C': return PSR_COINS1;
      case 'i': return PSR_ITEMS2;
      case 's': return PSR_COVER0;
      case 'S': return PSR_COVER1;
      case 'l': return PSR_LENIENT0;
      case 'L': return PSR_LENIENT1;
      case 't': return PSR_STRICT0;
      case 'T': return PSR_STRICT1;
    }
  }
  PyErr_SetString(PyExc_ValueError,
                  "invalid pattern set reduction method");
  return -1;                    /* return an error code */
}  /* get_red() */

/*----------------------------------------------------------------------
  Auxiliary Functions
----------------------------------------------------------------------*/

static void clean1 (PyObject *a)
{ if (a) Py_DECREF(a); }

/*--------------------------------------------------------------------*/

static void clean2 (PyObject *a, PyObject *b)
{ if (a) Py_DECREF(a); if (b) Py_DECREF(b); }

/*--------------------------------------------------------------------*/

static void clean3 (PyObject *a, PyObject *b, PyObject *c)
{ if (a) Py_DECREF(a); if (b) Py_DECREF(b); if (c) Py_DECREF(c); }

/*--------------------------------------------------------------------*/

static void clean4 (PyObject *a, PyObject *b, PyObject *c, PyObject *d)
{ if (a) Py_DECREF(a); if (b) Py_DECREF(b);
  if (c) Py_DECREF(c); if (d) Py_DECREF(d); }

/*--------------------------------------------------------------------*/

static void* error (PyObject *type, const char *msg)
{                               /* --- report an error */
  sig_remove();                 /* remove the signal handler */
  if (msg) PyErr_SetString(type, msg); /* set the error message */
  return NULL;                  /* return NULL as an error indicator */
}  /* error() */

/*--------------------------------------------------------------------*/

static ITEMBASE* ib_appPyObj (ITEMBASE *ibase, PyObject *appear)
{                               /* --- get item appearances */
  ITEM      k;                  /* item identifier */
  int       app;                /* item appearance indicator */
  PyObject  *ii;                /* item iterator (dictionary) */
  PyObject  *item;              /* to traverse the items */
  PyObject  *pyapp;             /* item appearance indicator */
  Py_hash_t h;                  /* hash value of item */

  assert(ibase);                /* check the function arguments */
  if (!appear) return ibase;    /* if no item apps. given, abort */
  if (!PyDict_Check(appear)) {  /* check for a dictionary */
    return ERR_TYPE("item appearances must be a dictionary"); }
  ii = PyObject_GetIter(appear);/* get an iterator for the dictionary */
  if (!ii) { clean1(ii);        /* and check for success */
    return ERR_TYPE("item appearances must be iterable"); }
  while ((item = PyIter_Next(ii)) != NULL) { /* traverse the items */
    if (item == Py_None)        /* if the item is 'None', */
      k = -1;                   /* it is the default appearance */
    else {                      /* if the item is a normal object */
      h = PyObject_Hash(item);  /* check whether item is hashable */
      if (h == -1) { clean2(item, ii);
        return ERR_TYPE("items must be hashable"); }
      k = ib_add(ibase, &item); /* add the item to the item base */
      if (k <   0) { clean2(item, ii); return ERR_MEM(); }
    }                           /* check for successful addition */
    pyapp = PyDict_GetItem(appear, item);
    Py_DECREF(item);            /* get the appearance indicator */
    if (PyUnicode_Check(pyapp)) {   /* if unicode string */
      pyapp = PyUnicode_AsUTF8String(pyapp);
      if (!pyapp) { clean1(ii); ERR_MEM(); }
      app = get_app(PyBytes_AS_STRING(pyapp)); }
    #if PY_MAJOR_VERSION < 3    /* if Python 2.x */
    else if (PyString_Check(pyapp)) /* if standard string */
      app = get_app(PyString_AS_STRING(pyapp));
    #endif                      /* get the string directly */
    else { clean1(ii);          /* if neither unicode nor string */
      return ERR_TYPE("item appearance indicators must be strings"); }
    if (app < 0) { clean1(ii); return NULL; }
    ib_setapp(ibase, k, app);   /* set the item appearance */
  }                             /* (default or per item) */
  Py_DECREF(ii);                /* drop the item iterator */
  return ibase;                 /* return the modified item base */
}  /* ib_appPyObj() */

/*--------------------------------------------------------------------*/

static TABAG* tbg_fromPyObj (PyObject *tracts, PyObject *appear)
{                               /* --- create a transaction bag */
  PyObject  *ti, *ii;           /* transaction and item iterator */
  PyObject  *trans;             /* to traverse the transactions */
  PyObject  *item;              /* to traverse the items */
  PyObject  *mul;               /* to get transaction multiplicity */
  Py_hash_t h;                  /* hash value of item */
  ITEM      k;                  /* number of items, buffers */
  SUPP      w;                  /* weight/support buffer */
  int       dict;               /* flag for transaction dictionary */
  TABAG     *tabag;             /* created transaction bag */
  ITEMBASE  *ibase;             /* underlying item base */

  assert(tracts);               /* check the function argument */
  ibase = ib_create(IB_OBJNAMES, 0, hashitem, cmpitems, NULL, NULL);
  if (!ibase) return ERR_MEM(); /* create an item base and add */
  if (!ib_appPyObj(ibase, appear)) {   /* the item appearances */
    ib_delete(ibase); return NULL; }
  tabag = tbg_create(ibase);    /* create a transaction bag */
  if (!tabag) { ib_delete(ibase); return ERR_MEM(); }
  dict = PyDict_Check(tracts);  /* check whether it is a dictionary */
  ti = PyObject_GetIter(tracts);/* get an iterator for transactions */
  if (!ti) { tbg_delete(tabag, 1);
    return ERR_TYPE("transaction database must be iterable"); }
  while ((trans = PyIter_Next(ti)) != NULL) {
    ib_clear(ibase);            /* traverse the transactions */
    ii = PyObject_GetIter(trans);
    if (!ii) { clean2(trans, ti); tbg_delete(tabag, 1);
      return ERR_TYPE("transactions must be iterable"); }
    w = 1;                      /* default: unit transaction weight */
    if (dict) {                 /* if trans. multiplicities given */
      mul = PyDict_GetItem(tracts, trans);
      if      (PyInt_Check (mul))  w = (SUPP)PyInt_AsLong (mul);
      else if (PyLong_Check(mul))  w = (SUPP)PyLong_AsLong(mul);
      else if (PyFloat_Check(mul)) w = (SUPP)PyFloat_AsDouble(mul);
      else { clean3(ii, trans, ti); tbg_delete(tabag, 1);
        return ERR_TYPE("transaction multiplicities must be numbers"); }
    }                           /* get transaction multiplicity */
    Py_DECREF(trans);           /* drop the transaction reference */
    while ((item = PyIter_Next(ii)) != NULL) {
      h = PyObject_Hash(item);  /* check whether item is hashable */
      if (h == -1) { clean3(item, ii, ti); tbg_delete(tabag, 1);
        return ERR_TYPE("items must be hashable"); }
      k = ib_add2ta(ibase, &item);
      Py_DECREF(item);          /* add item to transaction */
      if (k < 0) { clean2(ii, ti); tbg_delete(tabag, 1);
        return ERR_MEM(); }     /* check for an error */
    }
    Py_DECREF(ii);              /* drop the item iterator and */
    ib_finta(ibase, w);         /* set the transaction weight */
    if (PyErr_Occurred()) {     /* check for an iteration error */
      clean1(ti); tbg_delete(tabag, 1); return NULL; }
    if (tbg_addib(tabag) < 0) { /* add the transaction to the bag */
      clean1(ti); tbg_delete(tabag, 1); return ERR_MEM(); }
  }
  Py_DECREF(ti);                /* drop the transaction iterator */
  if (PyErr_Occurred()) { tbg_delete(tabag, 1); return NULL; }
  return tabag;                 /* return the created trans. bag */
}  /* tbg_fromPyObj() */

/*--------------------------------------------------------------------*/

static void* isr_pyborder (ISREPORT *rep, PyObject *border)
{                               /* --- set reporter filtering border */
  Py_ssize_t n;                 /* loop variable, sequence length */
  PyObject   *o;                /* to traverse the sequence elements */
  double     s;                 /* support threshold as a double */
  RSUPP      supp;              /* support threshold (scaled) */

  assert(rep && border);        /* check the function arguments */
  if (!border) return (void*)1; /* check for a border */
  if (!PySequence_Check(border))/* check for a sequence */
    return ERR_TYPE("border must be a list or tuple of numbers");
  n = PySequence_Length(border);/* get the sequence length */
  if (n <= 0)  return (void*)1; /* empty sequences need no processing */
  while (n > 0) { --n;          /* traverse the sequence elements */
    o = PySequence_GetItem(border, n);
    if      (PyLong_Check(o))   /* if element is a long integer */
      supp = (RSUPP)PyLong_AsLong(o);
    else if (PyInt_Check(o))    /* if element is an integer */
      supp = (RSUPP)PyInt_AsLong(o);
    else if (PyFloat_Check(o)){ /* if element is a float */
      s = PyFloat_AsDouble(o);
      supp = (s >= (double)SUPP_MAX) ? RSUPP_MAX : (RSUPP)s; }
    else { clean1(o);           /* get the element value as support */
      return ERR_TYPE("border elements must be numbers"); }
    Py_DECREF(o);               /* drop the element reference */
    if (isr_setbdr(rep, (ITEM)n, supp) < 0) return ERR_MEM();
  }                             /* set the border support */
  return (void*)1;              /* return 'ok' */
}  /* isr_pyborder() */

/*--------------------------------------------------------------------*/

static void isr_iset2PyObj (ISREPORT *rep, void *data)
{                               /* --- report an item set */
  int      v;                   /* loop variable, index offset */
  ITEM     i, n;                /* loop variables, number of items */
  RSUPP    supp, base;          /* item set support and base support */
  SUPP     s;                   /* support value for reporting */
  double   e, x;                /* evaluation and scaling factor */
  PyObject *pat;                /* pattern with item set and values */
  PyObject *iset;               /* found item set (as a tuple) */
  PyObject *vals;               /* values associated to item set */
  PyObject *obj;                /* to traverse the items, values */
  REPDATA  *rd = data;          /* report data structure */

  assert(rep && data);          /* check the function arguments */
  n    = isr_cnt(rep);          /* get the size of the item set */
  iset = PyTuple_New(n);        /* create an item set tuple */
  if (!iset) { rd->err = -1; return; }
  for (i = 0; i < n; i++) {     /* traverse the items */
    obj = (PyObject*)isr_itemobj(rep, isr_itemx(rep, i));
    Py_INCREF(obj);             /* get the corresp. Python object */
    PyTuple_SET_ITEM(iset, i, obj);
  }                             /* store the item in the set */
  if      (rd->mode == '[') vals = PyList_New (rd->cnt);
  else if (rd->mode == '(') vals = PyTuple_New(rd->cnt);
  else                      vals = PyTuple_New(rd->cnt+1);
  if (!vals) { clean1(iset); rd->err = -1; return; }
  supp = isr_supp(rep);         /* get the item set support */
  base = isr_suppx(rep, 0);     /* and the total transaction weight */
  s = 0; e = 0;                 /* initialize the report variables */
  for (v = 0; v < rd->cnt; v++){/* traverse the values to store */
    switch (rd->rep[v]) {       /* evaluate the value indicator */
      case 'a': s = (SUPP)supp;                 x =   0.0; break;
      case 's': e = (double)supp /(double)base; x =   1.0; break;
      case 'S': e = (double)supp /(double)base; x = 100.0; break;
      case 'p': e = isr_eval(rep);              x =   1.0; break;
      case 'P': e = isr_eval(rep);              x = 100.0; break;
      case 'e': e = isr_eval(rep);              x =   1.0; break;
      case 'E': e = isr_eval(rep);              x = 100.0; break;
      case 'Q': s = (SUPP)base;                 x =   0.0; break;
      default : s = 0;                          x =  -1.0; break;
    }                           /* get the requested value */
    if      (x <  0) obj = PyInt_FromLong(0L);
    else if (x <= 0) obj = PyObj_FromSupp(s);
    else             obj = PyFloat_FromDouble(e*x);
    if (!obj) { clean2(iset, vals); rd->err = -1; return; }
    if      (rd->mode == '[') PyList_SET_ITEM (vals, v,   obj);
    else if (rd->mode == '(') PyTuple_SET_ITEM(vals, v,   obj);
    else                      PyTuple_SET_ITEM(vals, v+1, obj);
  }                             /* store the created value */
  if (rd->mode == 0)            /* if to store in pattern object */
    pat = vals;                 /* just get the value tuple */
  else {                        /* if to store in separate object */
    pat = PyTuple_New(2);       /* create a pair of set and values */
    if (!pat) { clean2(iset, vals); rd->err = -1; return; }
    PyTuple_SET_ITEM(pat, 1, vals);
  }                             /* set the pair elements */
  PyTuple_SET_ITEM(pat, 0, iset);
  if (PyList_Append(rd->res, pat) != 0)
    rd->err = -1;               /* append pattern to the result list */
  Py_DECREF(pat);               /* remove reference to pattern */
}  /* isr_iset2PyObj() */

/*--------------------------------------------------------------------*/

static double lift (RSUPP supp, RSUPP body, RSUPP head, RSUPP base)
{                               /* --- compute lift value of a rule */
  return ((body <= 0) || (head <= 0)) ? 0
       : ((double)supp*(double)base) /((double)body*(double)head);
}  /* lift() */

/*--------------------------------------------------------------------*/

static void isr_rule2PyObj (ISREPORT *rep, void *data,
                            ITEM item, RSUPP body, RSUPP head)
{                               /* --- report an association rule */
  int      v;                   /* loop variable for values */
  ITEM     i, k, n;             /* loop variables, number of items */
  ITEM     z;                   /* to traverse the items */
  RSUPP    supp, base;          /* item set support and base support */
  SUPP     s;                   /* support value for reporting */
  double   e, x;                /* evaluation and scaling factor */
  PyObject *rule;               /* triplet of head, body and values */
  PyObject *ante;               /* antecedent of association rule */
  PyObject *cons;               /* consequent of association rule */
  PyObject *obj;                /* current item, to create objects */
  PyObject *vals;               /* values associated to rule */
  REPDATA  *rd = data;          /* report data structure */

  assert(rep && data            /* check the function arguments */
  &&    (body > 0) && (head > 0));
  assert(isr_uses(rep, item));  /* head item must be in item set */
  n    = isr_cnt(rep);          /* get the size of the item set */
  ante = PyTuple_New(n-1);      /* create a tuple for the rule body */
  if (!ante) { rd->err = -1; return; }
  for (i = k = 0; i < n; i++) { /* traverse the items */
    z = isr_itemx(rep, i);      /* get the next item and skip it */
    if (z == item) continue;    /* if it is the head of the rule */
    obj = (PyObject*)isr_itemobj(rep, z);
    Py_INCREF(obj);             /* get the corresp. Python object */
    PyTuple_SET_ITEM(ante, k, obj);
    k++;                        /* store the item in the rule body */
  }                             /* and advance the item position */
  if      (rd->mode == '[') vals = PyList_New (rd->cnt);
  else if (rd->mode == '(') vals = PyTuple_New(rd->cnt);
  else                      vals = PyTuple_New(rd->cnt+2);
  if (!vals) { clean1(ante); rd->err = -1; return; }
  supp = isr_supp(rep);         /* get the item set support */
  base = isr_suppx(rep, 0);     /* and the total transaction weight */
  s = 0; e = 0;                 /* initialize the report variables */
  for (v = 0; v < rd->cnt; v++){/* traverse the values to store */
    switch (rd->rep[v]) {       /* evaluate the value indicator */
      case 'a': s = (SUPP)supp;                   x =   0.0; break;
      case 'b': s = (SUPP)body;                   x =   0.0; break;
      case 'h': s = (SUPP)head;                   x =   0.0; break;
      case 's': e = (double)supp /(double)base;   x =   1.0; break;
      case 'S': e = (double)supp /(double)base;   x = 100.0; break;
      case 'x': e = (double)body /(double)base;   x =   1.0; break;
      case 'X': e = (double)body /(double)base;   x = 100.0; break;
      case 'y': e = (double)head /(double)base;   x =   1.0; break;
      case 'Y': e = (double)head /(double)base;   x = 100.0; break;
      case 'c': e = (double)supp /(double)body;   x =   1.0; break;
      case 'C': e = (double)supp /(double)body;   x = 100.0; break;
      case 'l': e = lift(supp, body, head, base); x =   1.0; break;
      case 'L': e = lift(supp, body, head, base); x = 100.0; break;
      case 'e': e = isr_eval(rep);                x =   1.0; break;
      case 'E': e = isr_eval(rep);                x = 100.0; break;
      case 'Q': s = (SUPP)base;                   x =   0.0; break;
      default : s = 0;                            x =  -1.0; break;
    }                           /* get the requested value */
    if      (x <  0) obj = PyInt_FromLong(0L);
    else if (x <= 0) obj = PyObj_FromSupp(s);
    else             obj = PyFloat_FromDouble(e*x);
    if (!obj) { clean2(ante, vals); rd->err = -1; return; }
    if      (rd->mode == '[') PyList_SET_ITEM (vals, v,   obj);
    else if (rd->mode == '(') PyTuple_SET_ITEM(vals, v,   obj);
    else                      PyTuple_SET_ITEM(vals, v+2, obj);
  }                             /* store the created value */
  cons = (PyObject*)isr_itemobj(rep, item);
  Py_INCREF(cons);              /* get rule head as a Python object */
  if (rd->mode == 0)            /* if to store in pattern object */
    rule = vals;                /* just get the value tuple */
  else {                        /* if to store in separate object */
    rule = PyTuple_New(3);      /* create triplet for rule and values */
    if (!rule) { clean3(cons, ante, vals); rd->err = -1; return; }
    PyTuple_SET_ITEM(rule, 2, vals);
  }                             /* set the triplet elements */
  PyTuple_SET_ITEM(rule, 0, cons);
  PyTuple_SET_ITEM(rule, 1, ante);
  if (PyList_Append(rd->res, rule) != 0)
    rd->err = -1;               /* append the pair to the result list */
  Py_DECREF(rule);              /* remove internal reference to rule */
}  /* isr_rule2PyObj() */

/*--------------------------------------------------------------------*/
#if FPSUPP

static PyObject* psp_toPyObj (PATSPEC *psp, double scale, int format)
{                               /* --- report pattern spectrum */
  int        e = 0;             /* error indicator */
  Py_ssize_t i;                 /* list size, result list index */
  ITEM       size;              /* loop variable for sizes */
  RSUPP      min, max;          /* minimum and maximum support */
  PyObject   *res;              /* created Python list object */
  PyObject   *cols[3] = {NULL,NULL,NULL};     /* result columns */
  PyObject   *z, *s, *m, *t;    /* (elements of) result triplet */

  assert(psp);                  /* check the function arguments */
  for (i = 0, size = psp_min(psp); size <= psp_max(psp); size++)
    if (psp_max4sz(psp, size) >= psp_min4sz(psp, size))
      i++;                      /* count the used sizes */
  if       (format == '#')      /* evaluate pattern spectrum format */
    res = PyDict_New();         /* dictionary (size,supp) -> freq */
  else if ((format == '=') || (format == '-'))
    res = PyList_New(i);        /* list of triplets (size,supp,freq) */
  else {                        /* three columns for size/min/max */
    res     = PyList_New(3); if (!res) return NULL;
    cols[0] = PyList_New(i); if (!cols[0]) { clean1(res); return NULL; }
    PyList_SET_ITEM(res, 0, cols[0]);  /* sizes */
    cols[1] = PyList_New(i); if (!cols[1]) { clean1(res); return NULL; }
    PyList_SET_ITEM(res, 1, cols[1]);  /* support values */
    cols[2] = PyList_New(i); if (!cols[2]) { clean1(res); return NULL; }
    PyList_SET_ITEM(res, 2, cols[2]);  /* frequencies */
  }                             /* (create the return object) */
  if (!res) return NULL;        /* check for successful creation */
  for (i = 0, size = psp_min(psp); size <= psp_max(psp); size++) {
    min = psp_min4sz(psp,size); /* traverse the pattern sizes */
    max = psp_max4sz(psp,size); /* get range of support values */
    if (max < min) continue;    /* skip size with empty support range */
    z = PyInt_FromLong((long)size);  /* create objects for storing */
    s = PyFloat_FromDouble((double)min);
    m = PyFloat_FromDouble((double)max);
    if (!z || !s || !m) { clean3(z, s, m); e = -1; break; }
    if      (format == '#') {   /* if pattern spectrum is dictionary */
      t = PyTuple_New(2);       /* create a tuple as value */
      if (!t) { clean3(z, s, m); e = -1; break; }
      PyTuple_SET_ITEM(t, 0, s);/* create value as pair (min,max) */
      PyTuple_SET_ITEM(t, 1, m);/* and map size to support range */
      PyDict_SetItem(res, z, t); }
    else if (format == '=') {   /* if pattern spectrum is list */
      t = PyTuple_New(3);       /* create a result list element */
      if (!t) { clean3(z, s, m); e = -1; break; }
      PyTuple_SET_ITEM(t, 0, z);/* fill element with triplet */
      PyTuple_SET_ITEM(t, 1, s);/* (size, supp_min, supp_max) */
      PyTuple_SET_ITEM(t, 2, m);/* and set it in the result */
      PyList_SET_ITEM(res, i, t); }
    else {                      /* if result is three columns */
      PyList_SET_ITEM(cols[0], i, z);
      PyList_SET_ITEM(cols[1], i, s);
      PyList_SET_ITEM(cols[2], i, m);
    }                           /* set size and support range */
    i += 1;                     /* count the processed pattern size */
  }
  if (e) { clean1(res); return NULL; }
  return res;                   /* return created pattern spectrum */
}  /* psp_toPyObj() */

/*--------------------------------------------------------------------*/
#else                           /* #if FPSUPP ... */

static PyObject* psp_toPyObj (PATSPEC *psp, double scale, int format)
{                               /* --- report pattern spectrum */
  int        e = 0;             /* error indicator */
  Py_ssize_t i;                 /* list size, result list index */
  ITEM       size;              /* loop variable for sizes */
  RSUPP      supp, min, max;    /* loop variable for supports */
  size_t     frq;               /* frequency of a pattern signature */
  PyObject   *res;              /* created Python list object */
  PyObject   *cols[3] = {NULL,NULL,NULL};     /* result columns */
  PyObject   *z, *s, *f, *t;    /* (elements of) result triplet */

  assert(psp);                  /* check the function arguments */
  i = (Py_ssize_t)psp_sigcnt(psp);
  if       (format == '#')      /* evaluate pattern spectrum format */
    res = PyDict_New();         /* dictionary (size,supp) -> freq */
  else if ((format == '=') || (format == '-'))
    res = PyList_New(i);        /* list of triplets (size,supp,freq) */
  else {                        /* three columns for size/supp/freq */
    res     = PyList_New(3); if (!res) return NULL;
    cols[0] = PyList_New(i); if (!cols[0]) { clean1(res); return NULL; }
    PyList_SET_ITEM(res, 0, cols[0]);  /* sizes */
    cols[1] = PyList_New(i); if (!cols[1]) { clean1(res); return NULL; }
    PyList_SET_ITEM(res, 1, cols[1]);  /* support values */
    cols[2] = PyList_New(i); if (!cols[2]) { clean1(res); return NULL; }
    PyList_SET_ITEM(res, 2, cols[2]);  /* frequencies */
  }                             /* (create the return object) */
  if (!res) return NULL;        /* check for successful creation */
  for (i = 0, size = psp_min(psp); size <= psp_max(psp); size++) {
    min = psp_min4sz(psp,size); /* traverse the pattern sizes */
    max = psp_max4sz(psp,size); /* get range of support values */
    if (max < min) continue;    /* skip size with empty support range */
    z = PyInt_FromLong((long)size);
    if (!z) { e = -1; break; }  /* create size object for storing */
    for (supp = min; supp <= max; supp++) {
      if ((frq = psp_getfrq(psp, size, supp)) <= 0)
        continue;               /* traverse the support values */
      s = PyInt_FromLong((long)supp);
      if (!s) {                 e = -1; break; }
      f = PyFloat_FromDouble((double)frq *scale);
      if (!f) {   clean1(s);    e = -1; break; }
      if      (format == '#') { /* if result is dictionary */
        t = PyTuple_New(2);     /* create a tuple as a key */
        if (!t) { clean2(f, s); e = -1; break; }
        Py_INCREF(z);           /* get a new size reference */
        PyTuple_SET_ITEM(t, 0, z);  /* create key (size, support) */
        PyTuple_SET_ITEM(t, 1, s);  /* and map key to frequency */
        PyDict_SetItem(res, t, f); }
      else if (format == '=') { /* if result is list of triplets */
        t = PyTuple_New(3);     /* create a result list element */
        if (!t) { clean2(f, s); e = -1; break; }
        Py_INCREF(z);           /* get a new size reference */
        PyTuple_SET_ITEM(t, 0, z);  /* fill element with triplet */
        PyTuple_SET_ITEM(t, 1, s);  /* (size, support, freq) and */
        PyTuple_SET_ITEM(t, 2, f);  /* set element in the result */
        PyList_SET_ITEM(res, i, t); }
      else {                    /* if result is three columns */
        Py_INCREF(z);           /* get a new size reference */
        PyList_SET_ITEM(cols[0], i, z);
        PyList_SET_ITEM(cols[1], i, s);
        PyList_SET_ITEM(cols[2], i, f);
      }                         /* set size, support, freq. in lists */
      i += 1;                   /* count the processed signature */
    }
    Py_DECREF(z);               /* initial reference no longer needed */
    if (e) break;               /* if an error occurred, abort loop */
  }
  if (e) { clean1(res); return NULL; }
  return res;                   /* return created pattern spectrum */
}  /* psp_toPyObj() */

#endif
/*--------------------------------------------------------------------*/

static int repinit (REPDATA *data, ISREPORT *isrep, CCHAR *report,
                    int target)
{                               /* --- initialize reporting */
  assert(data && isrep && report); /* check the function arguments */
  data->err = 0;                /* initialize the error indicator */
  if ((report[0] == '#')        /* if to get a pattern spectrum */
  ||  (report[0] == '=')        /* '#':     dictionary */
  ||  (report[0] == '-')        /* '=','-': list of triplets */
  ||  (report[0] == '|'))       /* '|':     three columns */
    return isr_addpsp(isrep, NULL);
  data->mode = ((report[0] == '(') || (report[0] == '['))
             ? *report++ : 0;   /* note mode and number of values */
  data->cnt  = (int)strlen(data->rep = report);
  data->res  = PyList_New(0);   /* create an empty output list */
  if (!data->res) return -1;    /* and set the reporting function */
  if (target & ISR_RULES) isr_setrule(isrep, isr_rule2PyObj, data);
  else                    isr_setrepo(isrep, isr_iset2PyObj, data);
  return 0;                     /* return 'ok' */
}  /* repinit() */

/*--------------------------------------------------------------------*/

static int repterm (REPDATA *data, ISREPORT *isrep, CCHAR *report)
{                               /* --- terminate reporting */
  assert(data && isrep && report); /* check the function arguments */
  if ((report[0] == '#')        /* if to get a pattern spectrum */
  ||  (report[0] == '=')        /* '#':     dictionary */
  ||  (report[0] == '-')        /* '=','-': list of triplets */
  ||  (report[0] == '|')) {     /* '|':     three columns */
    data->res = psp_toPyObj(isr_getpsp(isrep), 1.0, report[0]);
    return data->err = (data->res) ? 0 : -1;
  }                             /* make Python pattern spectrum */
  return data->err;             /* return the error status */
}  /* repterm() */

/*--------------------------------------------------------------------*/

static void repfn (long cnt, void *data)
{                               /* --- progress reporting function */
  if ((cnt > *(long*)data) && ((cnt % 20) == 0))
    fprintf(stderr, "%10ld\b\b\b\b\b\b\b\b\b\b", *(long*)data = cnt);
}  /* repfn() */

/*--------------------------------------------------------------------*/
/* fim (tracts, target='s', supp=10, zmin=1, zmax=None,               */
/*      report='a', eval='x', agg='x', thresh=10, border=None)        */
/*--------------------------------------------------------------------*/

static PyObject* py_fim (PyObject *self,
                         PyObject *args, PyObject *kwds)
{                               /* --- frequent item set mining */
  char     *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                        "report", "eval", "agg", "thresh", "border",
                         NULL };
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type as an identifier */
  double   supp    = 10;        /* minimum support of an item set */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    = 'x';       /* evaluation measure as identifier */
  CCHAR    *sagg   = "x";       /* aggregation mode as a string */
  int      agg     =  0;        /* aggregation mode as identifier */
  double   thresh  = 10;        /* threshold for evaluation measure */
  int      algo    = FPG_SIMPLE;   /* algorithm variant */
  int      mode    = FPG_DEFAULT;  /* operation mode/flags */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  FPGROWTH *fpgrowth;           /* fpgrowth miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sdllsssdO", ckwds,
        &tracts, &starg, &supp, &zmin, &zmax, &report,
        &seval, &sagg, &thresh, &border))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascmg");
  if (target < 0) return NULL;  /* translate the target string */
  if (zmin   < 0) return ERR_VALUE("zmin must be >= 0");
  if (zmax   < 0) zmax = LONG_MAX;  /* check the size range */
  if (zmax   < zmin) return ERR_VALUE("zmax must be >= zmin");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  eval = get_evalx(seval);      /* get evaluation measure and */
  if (eval   < 0) return NULL;  /* check whether it is valid */
  agg  = get_agg(sagg);         /* get aggregation mode and */
  if (agg    < 0) return NULL;  /* check whether it is valid */

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromPyObj(tracts, NULL);
  if (!tabag) return ERR_RTN(); /* get the given transactions */
  fpgrowth = fpg_create(target, supp, 100.0, 100.0,
                        (ITEM)zmin, (ITEM)zmax,
                        eval, agg, thresh, algo, mode);
  if (!fpgrowth) { tbg_delete(tabag, 1); return ERR_MEM(); }
  r = fpg_data(fpgrowth, tabag, 0, +2);
  if (r) fpg_delete(fpgrowth,1);/* prepare data for fpgrowth */
  if (r == -1) return ERR_MEM();/* check for error and no items */
  if (r <   0) { sig_remove(); return PyList_New(0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (fpg_report(fpgrowth, isrep) != 0)) {
    fpg_delete(fpgrowth, 1); return ERR_MEM(); }
  if (isr_pyborder(isrep, border) == NULL) {
    fpg_delete(fpgrowth, 1); return ERR_RTN(); }
  if ((repinit(&data, isrep, report, ISR_SETS) != 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    fpg_delete(fpgrowth, 1); return ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = fpg_mine(fpgrowth, ITEM_MIN, 0);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  fpg_delete(fpgrowth, 1);      /* delete the fpgrowth miner */
  if (sig_aborted()) { sig_abort(0);
    clean1(data.res); PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (r < 0) { clean1(data.res); return ERR_MEM(); }
  sig_remove();                 /* remove the signal handler */
  return data.res;              /* return the created result */
}  /* py_fim() */

/*--------------------------------------------------------------------*/
/* arules (tracts, supp=10, conf=80, zmin=1, zmax=None, report='aC',  */
/*         eval='x', thresh=10, mode='', appear=None)                 */
/*--------------------------------------------------------------------*/

static PyObject* py_arules (PyObject *self,
                            PyObject *args, PyObject *kwds)
{                               /* --- association rule mining */
  char     *ckwds[] = { "tracts", "supp", "conf", "zmin", "zmax",
                        "report", "eval", "thresh", "mode",
                        "appear", NULL };
  double   supp    = 10;        /* minimum support    of a rule */
  double   conf    = 80;        /* minimum confidence of a rule */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "aC";      /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    =  0;        /* evaluation measure as identifier */
  double   thresh  = 10;        /* threshold for evaluation measure */
  int      algo    = FPG_SINGLE;/* algorithm variant as identifier */
  CCHAR    *smode  = "";        /* operation mode/flags as a string */
  int      mode    = FPG_DEFAULT;  /* operation mode/flags */
  PyObject *appear = NULL;      /* item appearances indicators */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  FPGROWTH *fpgrowth;           /* fpgrowth miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|ddllssdsO", ckwds,
        &tracts, &supp, &conf, &zmin, &zmax, &report,
        &seval, &thresh, &smode, &appear))
    return NULL;                /* parse the function arguments */
  if ((conf < 0)                /* check the rule confidence */
  ||  (conf > 100)) return ERR_VALUE("invalid confidence");
  if (zmin  < 0)    return ERR_VALUE("zmin must be >= 0");
  if (zmax  < 0)    zmax = LONG_MAX; /* check the size range */
  if (zmax  < zmin) return ERR_VALUE("zmax must be >= zmin");
  if (zmin  > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax  > ITEM_MAX) zmax = ITEM_MAX;
  eval = get_evalx(seval);      /* get evaluation measure and */
  if (eval  < 0) return NULL;   /* check whether it is valid */
  if (strchr(smode, 'o') != NULL)
    mode |= FPG_ORIGSUPP;       /* original rule support definition */

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromPyObj(tracts, appear);
  if (!tabag) return ERR_RTN(); /* get the given transactions */
  fpgrowth = fpg_create(FPG_RULES, supp, 100.0, conf,
                        (ITEM)zmin, (ITEM)zmax,
                        eval, FPG_NONE, thresh, algo, mode);
  if (!fpgrowth) { tbg_delete(tabag, 1); return ERR_MEM(); }
  r = fpg_data(fpgrowth, tabag, 0, +2);
  if (r) fpg_delete(fpgrowth,1);/* prepare data for fpgrowth */
  if (r == -1) return ERR_MEM();/* check for error and no items */
  if (r <   0) { sig_remove(); return PyList_New(0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (fpg_report(fpgrowth, isrep)              != 0)
  ||  (repinit(&data, isrep, report, ISR_RULES) != 0)
  ||  (isr_setup(isrep) < 0)) { /* prepare the item set reporter */
    fpg_delete(fpgrowth, 1); return ERR_MEM(); }

  /* --- association rule mining --- */
  r = fpg_mine(fpgrowth, ITEM_MIN, 0);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  fpg_delete(fpgrowth, 1);      /* delete the fpgrowth miner */
  if (sig_aborted()) { sig_abort(0);
    clean1(data.res); PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (r < 0) { clean1(data.res); return ERR_MEM(); }
  sig_remove();                 /* remove the signal handler */
  return data.res;              /* return the created result */
}  /* py_arules() */

/*--------------------------------------------------------------------*/
/* apriori (tracts, target='s', supp=10, conf=80, zmin=1, zmax=None,  */
/*          report='a', eval='x', agg='x', thresh=10, prune=None,     */
/*          algo='a', mode='', border=None, appear=None)              */
/*--------------------------------------------------------------------*/

static PyObject* py_apriori (PyObject *self,
                             PyObject *args, PyObject *kwds)
{                               /* --- Apriori algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "conf",
                        "zmin", "zmax", "report",
                        "eval", "agg", "thresh", "prune",
                        "algo", "mode", "border", "appear", NULL };
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type as an identifier */
  double   supp    = 10;        /* minimum support of an item set */
  double   conf    = 80;        /* minimum confidence of a rule */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    =  0;        /* evaluation measure as identifier */
  CCHAR    *sagg   = "x";       /* aggregation mode as a string */
  int      agg     =  0;        /* aggregation mode as identifier */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "a";       /* algorithm as a string */
  int      algo    = APR_AUTO;  /* algorithm as an identifier */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = APR_DEFAULT;  /* operation mode/flags */
  long     prune   = LONG_MIN;  /* min. size for evaluation filtering */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *appear = NULL;      /* item appearance indicators */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  APRIORI  *apriori;            /* apriori miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sddllsssdlssOO",ckwds,
        &tracts, &starg, &supp, &conf, &zmin, &zmax, &report, &seval,
        &sagg, &thresh, &prune, &salgo, &smode, &border, &appear))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascmgr");
  if (target < 0) return NULL;  /* translate the target string */
  if (zmin   < 0) return ERR_VALUE("zmin must be >= 0");
  if (zmax   < 0) zmax = LONG_MAX; /* check the size range */
  if (zmax   < zmin) return ERR_VALUE("zmax must be >= zmin");
  if (zmin > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax > ITEM_MAX) zmax = ITEM_MAX;
  eval = get_evalx(seval);      /* get evaluation measure and */
  if (eval   < 0) return NULL;  /* check whether it is valid */
  if (eval   <= RE_NONE) prune = LONG_MIN;
  agg  = get_agg(sagg);         /* get aggregation mode and */
  if (agg    < 0) return NULL;  /* check whether it is valid */
  if (salgo[0] && salgo[1]) {   /* if textual identifier */
    if      (strcmp(salgo, "auto")   == 0) salgo = "a";
    else if (strcmp(salgo, "basic")  == 0) salgo = "b";
    else                                   algo  = -1;
  }                             /* translate the algorithm string */
  if (salgo[0] && !salgo[1]) {  /* if single character, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 'a': algo = APR_AUTO;  break;
      case 'b': algo = APR_BASIC; break;
      default : algo = -1;        break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) return ERR_VALUE("invalid Apriori algorithm variant");
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'o') mode |=  APR_ORIGSUPP;
    else if (*s == 'x') mode &= ~APR_PERFECT;
    else if (*s == 't') mode &= ~APR_TATREE;
    else if (*s == 'T') mode &= ~APR_TATREE;
    else if (*s == 'y') mode |=  APR_POST;
    else if (*s == 'z') eval |=  APR_INVBXS;
  }                             /* adapt the operation mode */

  /* --- get and prepare transactions --- */
  sig_install();                /* install the signal handler */
  if (!(target & ISR_RULES)) appear = NULL;
  tabag = tbg_fromPyObj(tracts, appear);
  if (!tabag) return ERR_RTN(); /* get transactions, create miner */
  apriori = apriori_create(target, supp, 100.0, conf,
                           (ITEM)zmin, (ITEM)zmax,
                           eval, agg, thresh, algo, mode);
  if (!apriori) { tbg_delete(tabag, 1); return ERR_MEM(); }
  r = apriori_data(apriori, tabag, 0, +2);
  if (r) apriori_delete(apriori, 1);   /* prepare data for apriori */
  if (r == -1) return ERR_MEM();/* check for an error and no items */
  if (r <   0) { sig_remove(); return PyList_New(0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  || (apriori_report(apriori, isrep) != 0)) {
    apriori_delete(apriori, 1); return ERR_MEM(); }
  if (isr_pyborder(isrep, border) == NULL) {
    apriori_delete(apriori, 1); return ERR_RTN(); }
  if ((repinit(&data, isrep, report, target) != 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    apriori_delete(apriori, 1); return ERR_MEM(); }

  /* --- frequent item set/association rule mining --- */
  if (prune < ITEM_MIN) prune = ITEM_MIN;
  if (prune > ITEM_MAX) prune = ITEM_MAX;
  r = apriori_mine(apriori, (ITEM)prune, 1.0, 0);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  apriori_delete(apriori, 1);   /* delete the apriori miner */
  if (sig_aborted()) { sig_abort(0);
    clean1(data.res); PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (r < 0) { clean1(data.res); return ERR_MEM(); }
  sig_remove();                 /* remove the signal handler */
  return data.res;              /* return the created result */
}  /* py_apriori() */

/*--------------------------------------------------------------------*/
/* eclat (tracts, target='s', supp=10, conf=80, zmin=1, zmax=None,    */
/*        report='a', eval='x', agg='x', thresh=10, prune=None,       */
/*        algo='a', mode='', border=None, appear=None)                */
/*--------------------------------------------------------------------*/

static PyObject* py_eclat (PyObject *self,
                           PyObject *args, PyObject *kwds)
{                               /* --- Eclat algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "conf",
                        "zmin", "zmax", "report",
                        "eval", "agg", "thresh", "prune",
                        "algo", "mode", "border", "appear", NULL };
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type */
  double   supp    = 10;        /* minimum support of an item set */
  double   conf    = 80;        /* minimum confidence of a rule */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    =  0;        /* evaluation measure as identifier */
  CCHAR    *sagg   = "x";       /* aggregation mode as a string */
  int      agg     =  0;        /* aggregation mode as identifier */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "a";       /* algorithm as a string */
  int      algo    = ECL_OCCDLV;/* algorithm as an identifier */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = ECL_DEFAULT;  /* operation mode/flags */
  long     prune   = LONG_MIN;  /* min. size for evaluation filtering */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *appear = NULL;      /* item appearance indicators */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  ECLAT    *eclat;              /* eclat miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sddllsssdlssOO",ckwds,
        &tracts, &starg, &supp, &conf, &zmin, &zmax, &report, &seval,
        &sagg, &thresh, &prune, &salgo, &smode, &border, &appear))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascmgr");
  if (target < 0) return NULL;  /* translate the target string */
  if ((conf  < 0)               /* check the rule confidence */
  ||  (conf  > 100)) return ERR_VALUE("invalid confidence");
  if (zmin   < 0)    return ERR_VALUE("zmin must be >= 0");
  if (zmax   < 0)    zmax = LONG_MAX; /* check the size range */
  if (zmax   < zmin) return ERR_VALUE("zmax must be >= zmin");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  eval = get_evalx(seval);      /* get evaluation measure and */
  if (eval   < 0) return NULL;  /* check whether it is valid */
  if (eval   <= RE_NONE) prune = LONG_MIN;
  agg  = get_agg(sagg);         /* get the aggregation mode */
  if (agg    < 0) return NULL;
  if (salgo[0] && salgo[1]) {   /* if textual identifier */
    if      (strcmp(salgo, "auto")   == 0) salgo = "a";
    else if (strcmp(salgo, "basic")  == 0) salgo = "e";
    else if (strcmp(salgo, "lists")  == 0) salgo = "i";
    else if (strcmp(salgo, "tids")   == 0) salgo = "i";
    else if (strcmp(salgo, "bits")   == 0) salgo = "b";
    else if (strcmp(salgo, "table")  == 0) salgo = "t";
    else if (strcmp(salgo, "simple") == 0) salgo = "s";
    else if (strcmp(salgo, "ranges") == 0) salgo = "r";
    else if (strcmp(salgo, "occdlv") == 0) salgo = "o";
    else if (strcmp(salgo, "diff")   == 0) salgo = "d";
    else                                   algo  = -1;
  }                             /* translate the algorithm string */
  if (salgo[0] && !salgo[1]) {  /* if single character, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 'a': algo = ECL_AUTO;   break;
      case 'e': algo = ECL_BASIC;  break;
      case 'i': algo = ECL_LISTS;  break;
      case 'b': algo = ECL_BITS;   break;
      case 't': algo = ECL_TABLE;  break;
      case 's': algo = ECL_SIMPLE; break;
      case 'r': algo = ECL_RANGES; break;
      case 'o': algo = ECL_OCCDLV; break;
      case 'd': algo = ECL_DIFFS;  break;
      default : algo = -1;         break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) return ERR_VALUE("invalid Eclat algorithm");
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'o') mode |=  ECL_ORIGSUPP;
    else if (*s == 'l') mode &= ~ECL_FIM16;
    else if (*s == 'x') mode &= ~ECL_PERFECT;
    else if (*s == 'i') mode &= ~ECL_REORDER;
    else if (*s == 'u') mode &= ~ECL_TAIL;
    else if (*s == 'y') mode |=  ECL_HORZ;
    else if (*s == 'Y') mode |=  ECL_VERT;
    else if (*s == 'z') eval |=  ECL_INVBXS;
  }                             /* adapt the operation mode */

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  if (!(target & ISR_RULES)) appear = NULL;
  tabag = tbg_fromPyObj(tracts, appear);
  if (!tabag) return ERR_RTN(); /* get the given transactions */
  eclat = eclat_create(target, supp, 100.0, conf,
                       (ITEM)zmin, (ITEM)zmax,
                       eval, agg, thresh, algo, mode);
  if (!eclat) { tbg_delete(tabag, 1); return ERR_MEM(); }
  r = eclat_data(eclat, tabag, 0, +2);
  if (r) eclat_delete(eclat, 1);/* prepare data for eclat */
  if (r == -1) return ERR_MEM();/* check for an error and no items */
  if (r <   0) { sig_remove(); return PyList_New(0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (eclat_report(eclat, isrep) != 0)) {
    eclat_delete(eclat, 1); return ERR_MEM(); }
  if (isr_pyborder(isrep, border) == NULL) {
    eclat_delete(eclat, 1); return ERR_RTN(); }
  if ((repinit(&data, isrep, report, target) != 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    eclat_delete(eclat, 1); return ERR_MEM(); }

  /* --- frequent item set/association rule mining --- */
  if (prune < ITEM_MIN) prune = ITEM_MIN;
  if (prune > ITEM_MAX) prune = ITEM_MAX;
  r = eclat_mine(eclat, (ITEM)prune, 0);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  eclat_delete(eclat, 1);       /* delete the eclat miner */
  if (sig_aborted()) { sig_abort(0);
    clean1(data.res); PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (r < 0) { clean1(data.res); return ERR_MEM(); }
  sig_remove();                 /* remove the signal handler */
  return data.res;              /* return the created result */
}  /* py_eclat() */

/*--------------------------------------------------------------------*/
/* fpgrowth (tracts, target='s', supp=10, conf=80, zmin=1, zmax=None, */
/*           report='a', eval='x', agg='x', thresh=10, prune=None,    */
/*           algo='s', mode='', border=None, appear=None)             */
/*--------------------------------------------------------------------*/

static PyObject* py_fpgrowth (PyObject *self,
                              PyObject *args, PyObject *kwds)
{                               /* --- FP-growth algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "conf",
                        "zmin", "zmax", "report",
                        "eval", "agg", "thresh", "prune",
                        "algo", "mode", "border", "appear", NULL };
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type as an identifier */
  double   supp    = 10;        /* minimum support of an item set */
  double   conf    = 80;        /* minimum confidence of a rule */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    = 'x';       /* evaluation measure as identifier */
  CCHAR    *sagg   = "x";       /* aggregation mode as a string */
  int      agg     =  0;        /* aggregation mode as identifier */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "s";       /* algorithm as a string */
  int      algo    = FPG_SIMPLE;/* algorithm as an identifier */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = FPG_DEFAULT;  /* operation mode/flags */
  long     prune   = LONG_MIN;  /* min. size for evaluation filtering */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *appear = NULL;      /* item appearance indicators */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  FPGROWTH *fpgrowth;           /* fpgrowth miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sddllsssdlssOO",ckwds,
        &tracts, &starg, &supp, &conf, &zmin, &zmax, &report, &seval,
        &sagg, &thresh, &prune, &salgo, &smode, &border, &appear))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascmgr");
  if (target < 0) return NULL;  /* translate the target string */
  if (zmin   < 0) return ERR_VALUE("zmin must be >= 0");
  if (zmax   < 0) zmax = LONG_MAX; /* check the size range */
  if (zmax   < zmin) return ERR_VALUE("zmax must be >= zmin");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  eval = get_evalx(seval);      /* get the evaluation measure */
  if (eval   < 0) return NULL;
  if (eval   <= RE_NONE) prune = LONG_MIN;
  agg  = get_agg(sagg);         /* get the aggregation mode */
  if (agg    < 0) return NULL;
  if (salgo[0] && salgo[1]) {   /* if textual identifier */
    if      (strcmp(salgo, "simple")  == 0) salgo = "s";
    else if (strcmp(salgo, "complex") == 0) salgo = "c";
    else if (strcmp(salgo, "single")  == 0) salgo = "d";
    else if (strcmp(salgo, "topdown") == 0) salgo = "t";
    else                                    algo  = -1;
  }                             /* translate the algorithm string */
  if (salgo[0] && !salgo[1]) {  /* if single character, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 's': algo = FPG_SIMPLE;  break;
      case 'c': algo = FPG_COMPLEX; break;
      case 'd': algo = FPG_SINGLE;  break;
      case 't': algo = FPG_TOPDOWN; break;
      default : algo = -1;          break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) return ERR_VALUE("invalid FP-growth algorithm");
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'o') mode |=  FPG_ORIGSUPP;
    else if (*s == 'l') mode &= ~FPG_FIM16;
    else if (*s == 'x') mode &= ~FPG_PERFECT;
    else if (*s == 'i') mode &= ~FPG_REORDER;
    else if (*s == 'u') mode &= ~FPG_TAIL;
    else if (*s == 'z') eval |=  FPG_INVBXS;
  }                             /* adapt the operation mode */

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  if (!(target & ISR_RULES)) appear = NULL;
  tabag = tbg_fromPyObj(tracts, appear);
  if (!tabag) return ERR_RTN(); /* get the given transactions */
  fpgrowth = fpg_create(target, supp, 100.0, conf,
                        (ITEM)zmin, (ITEM)zmax,
                        eval, agg, thresh, algo, mode);
  if (!fpgrowth) { tbg_delete(tabag, 1); return ERR_MEM(); }
  r = fpg_data(fpgrowth, tabag, 0, +2);
  if (r) fpg_delete(fpgrowth,1);/* prepare data for fpgrowth */
  if (r == -1) return ERR_MEM();/* check for error and no items */
  if (r <   0) { sig_remove(); return PyList_New(0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (fpg_report(fpgrowth, isrep) != 0)) {
    fpg_delete(fpgrowth, 1); return ERR_MEM(); }
  if (isr_pyborder(isrep, border) == NULL) {
    fpg_delete(fpgrowth, 1); return ERR_RTN(); }
  if ((repinit(&data, isrep, report, target) != 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    fpg_delete(fpgrowth, 1); return ERR_MEM(); }

  /* --- frequent item set mining --- */
  if (prune < ITEM_MIN) prune = ITEM_MIN;
  if (prune > ITEM_MAX) prune = ITEM_MAX;
  r = fpg_mine(fpgrowth, (ITEM)prune, 0);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  fpg_delete(fpgrowth, 1);      /* delete the fpgrowth miner */
  if (sig_aborted()) { sig_abort(0);
    clean1(data.res); PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (r < 0) { clean1(data.res); return ERR_MEM(); }
  sig_remove();                 /* remove the signal handler */
  return data.res;              /* return the created result */
}  /* py_fpgrowth() */

/*--------------------------------------------------------------------*/
/* sam (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',   */
/*      eval='x', thresh=10, algo='b', mode='', border=None)          */
/*--------------------------------------------------------------------*/

static PyObject* py_sam (PyObject *self,
                         PyObject *args, PyObject *kwds)
{                               /* --- SaM algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                        "report", "eval", "thresh", "algo", "mode",
                        "border", NULL };
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type as an identifier */
  double   supp    = 10;        /* minimum support of an item set */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    = 'x';       /* evaluation measure as identifier */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "b";          /* algorithm as a string */
  int      algo    = SAM_BSEARCH;  /* algorithm as an identifier */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = SAM_DEFAULT;  /* operation mode/flags */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  SAM      *sam;                /* split and merge miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sdllssdssO", ckwds,
        &tracts, &starg, &supp, &zmin, &zmax, &report,
        &seval, &thresh, &salgo, &smode, &border))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascm");
  if (target < 0) return NULL;  /* translate the target string */
  if (zmin   < 0) return ERR_VALUE("zmin must be >= 0");
  if (zmax   < 0) zmax = LONG_MAX; /* check the size range */
  if (zmax   < zmin) return ERR_VALUE("zmax must be >= zmin");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  eval = get_eval(seval);       /* get evaluation measure and */
  if (eval   < 0) return NULL;  /* check whether it is valid */
  if (salgo[0] && salgo[1]) {   /* if textual identifier */
    if      (strcmp(salgo, "basic")   == 0) salgo = "s";
    else if (strcmp(salgo, "simple")  == 0) salgo = "s";
    else if (strcmp(salgo, "bsearch") == 0) salgo = "b";
    else if (strcmp(salgo, "double")  == 0) salgo = "d";
    else if (strcmp(salgo, "tree")    == 0) salgo = "t";
    else                                    algo  = -1;
  }                             /* translate the algorithm string */
  if (salgo[0] && !salgo[1]) {  /* if single character, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 's': algo = SAM_BASIC;   break;
      case 'b': algo = SAM_BSEARCH; break;
      case 'd': algo = SAM_DOUBLE;  break;
      case 't': algo = SAM_TREE;    break;
      default : algo = -1;          break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) return ERR_VALUE("invalid SaM algorithm");
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'l') mode &= ~SAM_FIM16;
    else if (*s == 'x') mode &= ~SAM_PERFECT;
  }                             /* adapt the operation mode */

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromPyObj(tracts, NULL);
  if (!tabag) return ERR_RTN(); /* get the given transactions */
  sam = sam_create(target, supp, 0.0, (ITEM)zmin, (ITEM)zmax,
                   0, -1.0, eval, thresh, algo, mode);
  if (!sam) { tbg_delete(tabag, 1); return ERR_MEM(); }
  r = sam_data(sam, tabag, +2); /* create a split and merge miner */
  if (r) sam_delete(sam, 1);    /* prepare data for split and merge */
  if (r == -1) return ERR_MEM();/* check for error and no items */
  if (r <   0) { sig_remove(); return PyList_New(0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (sam_report(sam, isrep) != 0)) {
    sam_delete(sam, 1); return ERR_MEM(); }
  if (isr_pyborder(isrep, border) == NULL) {
    sam_delete(sam, 1); return ERR_RTN(); }
  if ((repinit(&data, isrep, report, target) != 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    sam_delete(sam, 1); return ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = sam_mine(sam, 8192);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  sam_delete(sam, 1);           /* delete the split and merge miner */
  if (sig_aborted()) { sig_abort(0);
    clean1(data.res); PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (r < 0) { clean1(data.res); return ERR_MEM(); }
  sig_remove();                 /* remove the signal handler */
  return data.res;              /* return the created result */
}  /* py_sam() */

/*--------------------------------------------------------------------*/
/* relim (tracts, target='s', supp=10, zmin=1, zmax=None, report='a', */
/*        eval='x', thresh=10, algo='s', mode='', border=None)        */
/*--------------------------------------------------------------------*/

static PyObject* py_relim (PyObject *self,
                           PyObject *args, PyObject *kwds)
{                               /* --- RElim algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                        "report", "eval", "thresh", "algo", "mode",
                        "border", NULL };
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type as an identifier */
  double   supp    = 10;        /* minimum support of an item set */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    = 'x';       /* evaluation measure as identifier */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "s";       /* algorithm as a string */
  int      algo    = REL_BASIC; /* algorithm as an identifier */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = REL_DEFAULT|REL_FIM16; /* operation mode/flags */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  RELIM    *relim;              /* relim miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sdllssdssO", ckwds,
        &tracts, &starg, &supp, &zmin, &zmax, &report,
        &seval, &thresh, &salgo, &smode, &border))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascm");
  if (target < 0) return NULL;  /* translate the target string */
  if (zmin   < 0)    return ERR_VALUE("zmin must be >= 0");
  if (zmax   < 0)    zmax = LONG_MAX;  /* check the size range */
  if (zmax   < zmin) return ERR_VALUE("zmax must be >= zmin");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  eval = get_eval(seval);       /* get evaluation measure and */
  if (eval   < 0) return NULL;  /* check whether it is valid */
  if (salgo[0] && salgo[1]) {   /* if textual identifier */
    if      (strcmp(salgo, "basic")  == 0) salgo = "s";
    else if (strcmp(salgo, "simple") == 0) salgo = "s";
    else                                   algo  = -1;
  }                             /* translate the algorithm string */
  if (salgo[0] && !salgo[1]) {  /* if single character, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 's': algo = REL_BASIC; break;
      default : algo = -1;        break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) return ERR_VALUE("invalid RElim algorithm");
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'l') mode &= ~REL_FIM16;
    else if (*s == 'x') mode &= ~REL_PERFECT;
  }                             /* adapt the operation mode */

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromPyObj(tracts, NULL);
  if (!tabag) return ERR_RTN(); /* get the given transactions */
  relim = relim_create(target, supp, 0.0, (ITEM)zmin, (ITEM)zmax,
                       0, -1.0, eval, thresh, algo, mode);
  if (!relim) { tbg_delete(tabag, 1); return ERR_MEM(); }
  r = relim_data(relim, tabag, +2);
  if (r) relim_delete(relim,1); /* prepare data for relim */
  if (r == -1) return ERR_MEM();/* check for error and no items */
  if (r <   0) { sig_remove(); return PyList_New(0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  || (relim_report(relim, isrep) != 0)) {
    relim_delete(relim, 1); return ERR_MEM(); }
  if (isr_pyborder(isrep, border) == NULL) {
    relim_delete(relim, 1); return ERR_RTN(); }
  if ((repinit(&data, isrep, report, target)!= 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    relim_delete(relim, 1); return ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = relim_mine(relim, 32);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  relim_delete(relim, 1);       /* delete the relim miner */
  if (sig_aborted()) { sig_abort(0);
    clean1(data.res); PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (r < 0) { clean1(data.res); return ERR_MEM(); }
  sig_remove();                 /* remove the signal handler */
  return data.res;              /* return the created result */
}  /* py_relim() */

/*--------------------------------------------------------------------*/
/* carpenter (tracts, target='s', supp=10, zmin=1, zmax=None,         */
/*            report='a', eval='x', thresh=10, algo='a', mode='',     */
/*            border=None)                                            */
/*--------------------------------------------------------------------*/

static PyObject* py_carpenter (PyObject *self,
                               PyObject *args, PyObject *kwds)
{                               /* --- Carpenter algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                        "report", "eval", "thresh", "algo", "mode",
                        "border", NULL };
  CCHAR    *starg  = "c";       /* target type as a string */
  int      target  = ISR_CLOSED;/* target type as an identifier */
  double   supp    = 10;        /* minimum support of an item set */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    = 'x';       /* evaluation measure as identifier */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "a";       /* algorithm as a string */
  int      algo    = CARP_AUTO; /* algorithm as an identifier */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = CARP_DEFAULT; /* operation mode/flags */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  CARP     *carp;               /* carpenter miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sdllssdssO", ckwds,
        &tracts, &starg, &supp, &zmin, &zmax, &report,
        &seval, &thresh, &salgo, &smode, &border))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "cm");
  if (target < 0) return NULL;  /* translate the target string */
  if ((target != ISR_CLOSED) && (target != IST_MAXIMAL)) {
    PyErr_SetString(PyExc_ValueError, "invalid target type");
    return NULL;                /* carpenter only supports mining */
  }                             /* closed and maximal item sets */
  if (zmin < 0)    return ERR_VALUE("zmin must be >= 0");
  if (zmax < 0)    zmax = LONG_MAX;  /* check the size range */
  if (zmax < zmin) return ERR_VALUE("zmax must be >= zmin");
  if (zmin > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax > ITEM_MAX) zmax = ITEM_MAX;
  eval = get_eval(seval);       /* get evaluation measure and */
  if (eval < 0) return NULL;    /* check whether it is valid */
  if (salgo[0] && salgo[1]) {   /* if textual identifier */
    if      (strcmp(salgo, "auto")    == 0) salgo = "a";
    else if (strcmp(salgo, "table")   == 0) salgo = "t";
    else if (strcmp(salgo, "table")   == 0) salgo = "t";
    else if (strcmp(salgo, "tids")    == 0) salgo = "l";
    else if (strcmp(salgo, "tidlist") == 0) salgo = "l";
    else if (strcmp(salgo, "list")    == 0) salgo = "l";
    else                                    algo  = -1;
  }                             /* translate the algorithm string */
  if (salgo[0] && !salgo[1]) {  /* if single character, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 'a': algo = CARP_AUTO;    break;
      case 't': algo = CARP_TABLE;   break;
      case 'l': algo = CARP_TIDLIST; break;
      default : algo = -1;           break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) return ERR_VALUE("invalid Carpenter algorithm");
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'x') mode &= ~CARP_PERFECT;
    else if (*s == 'z') mode |=  CARP_FILTER;
    else if (*s == 'y') mode &= ~CARP_MAXONLY;
    else if (*s == 'p') mode &= ~CARP_COLLATE;
  }                             /* adapt the operation mode */

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromPyObj(tracts, NULL);
  if (!tabag) return ERR_RTN(); /* get the given transactions */
  carp = carp_create(target, supp, 100.0, (ITEM)zmin, (ITEM)zmax,
                     eval, thresh, algo, mode);
  if (!carp) { tbg_delete(tabag, 1); return ERR_MEM(); }
  r = carp_data(carp, tabag, -2);
  if (r) carp_delete(carp, 1);  /* prepare data for carpenter */
  if (r == -1) return ERR_MEM();/* check for error and no items */
  if (r <   0) { sig_remove(); return PyList_New(0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (carp_report(carp, isrep) != 0)) {
    carp_delete(carp, 1); return ERR_MEM(); }
  if (isr_pyborder(isrep, border) == NULL) {
    carp_delete(carp, 1); return ERR_RTN(); }
  if ((repinit(&data, isrep, report, target) != 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    carp_delete(carp, 1); return ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = carp_mine(carp);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  carp_delete(carp, 1);         /* delete the carpenter miner */
  if (sig_aborted()) { sig_abort(0);
    clean1(data.res); PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (r < 0) { clean1(data.res); return ERR_MEM(); }
  sig_remove();                 /* remove the signal handler */
  return data.res;              /* return the created result */
}  /* py_carpenter() */

/*--------------------------------------------------------------------*/
/* ista (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',  */
/*       eval='x', thresh=10, algo='x', mode='', border=None)         */
/*--------------------------------------------------------------------*/

static PyObject* py_ista (PyObject *self,
                          PyObject *args, PyObject *kwds)
{                               /* --- IsTa algorithm */
  char     *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                        "report", "eval", "thresh", "algo", "mode",
                        "border", NULL };
  CCHAR    *starg  = "c";       /* target type as a string */
  int      target  = ISR_CLOSED;/* target type as an identifier */
  double   supp    = 10;        /* minimum support of an item set */
  long     zmin    =  1;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "a";       /* indicators of values to report */
  CCHAR    *seval  = "x";       /* evaluation measure as a string */
  int      eval    = 'x';       /* evaluation measure as identifier */
  double   thresh  = 10;        /* threshold for evaluation measure */
  CCHAR    *salgo  = "x";       /* algorithm as a string */
  int      algo    = ISTA_PREFIX;  /* algorithm as an identifier */
  CCHAR    *smode  = "", *s;    /* operation mode/flags as a string */
  int      mode    = ISTA_DEFAULT; /* operation mode/flags */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  ISTA     *ista;               /* ista miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sdllssdssO", ckwds,
        &tracts, &starg, &supp, &zmin, &zmax, &report,
        &seval, &thresh, &salgo, &smode, &border))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "cm");
  if (target < 0) return NULL;  /* translate the target string */
  if ((target != ISR_CLOSED) && (target != IST_MAXIMAL)) {
    PyErr_SetString(PyExc_ValueError, "invalid target type");
    return NULL;                /* carpenter only supports mining */
  }                             /* closed and maximal item sets */
  if (zmin < 0)    return ERR_VALUE("zmin must be >= 0");
  if (zmax < 0)    zmax = LONG_MAX;  /* check the size range */
  if (zmax < zmin) return ERR_VALUE("zmax must be >= zmin");
  if (zmin > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax > ITEM_MAX) zmax = ITEM_MAX;
  eval = get_eval(seval);       /* get evaluation measure and */
  if (eval < 0) return NULL;    /* check whether it is valid */
  if (salgo[0] && salgo[1]) {   /* if textual identifier */
    if      (strcmp(salgo, "auto")     == 0) salgo = "a";
    else if (strcmp(salgo, "pfx")      == 0) salgo = "x";
    else if (strcmp(salgo, "prefix")   == 0) salgo = "x";
    else if (strcmp(salgo, "pat")      == 0) salgo = "p";
    else if (strcmp(salgo, "patricia") == 0) salgo = "p";
    else                                     algo  = -1;
  }                             /* translate the algorithm string */
  if (salgo[0] && !salgo[1]) {  /* if single character, */
    switch (salgo[0]) {         /* evaluate the algorithm code */
      case 'a': algo = ISTA_AUTO;     break;
      case 'x': algo = ISTA_PREFIX;   break;
      case 'p': algo = ISTA_PATRICIA; break;
      default : algo = -1;            break;
    }                           /* set an error code for all */
  }                             /* other algorithm indicators */
  if (algo < 0) return ERR_VALUE("invalid IsTa algorithm");
  for (s = smode; *s; s++) {    /* traverse the mode characters */
    if      (*s == 'p') mode &= ~ISTA_PRUNE;
    else if (*s == 'z') mode |=  ISTA_FILTER;
  }                             /* adapt the operation mode */

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromPyObj(tracts, NULL);
  if (!tabag) return ERR_RTN(); /* get the given transactions */
  ista = ista_create(target, supp, 100.0, (ITEM)zmin, (ITEM)zmax,
                     eval, thresh, algo, mode);
  if (!ista) { tbg_delete(tabag, 1); return ERR_MEM(); }
  r = ista_data(ista, tabag, -2);
  if (r) ista_delete(ista, 1);  /* prepare data for ista */
  if (r == -1) return ERR_MEM();/* check for error and no items */
  if (r <   0) { sig_remove(); return PyList_New(0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (ista_report(ista, isrep) != 0)) {
    ista_delete(ista, 1); return ERR_MEM(); }
  if (isr_pyborder(isrep, border) == NULL) {
    ista_delete(ista, 1); return ERR_RTN(); }
  if ((repinit(&data, isrep, report, target)!= 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    ista_delete(ista, 1); return ERR_MEM(); }

  /* --- frequent item set mining --- */
  r = ista_mine(ista);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  ista_delete(ista, 1);         /* delete the ista miner */
  if (sig_aborted()) { sig_abort(0);
    clean1(data.res); PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (r < 0) { clean1(data.res); return ERR_MEM(); }
  sig_remove();                 /* remove the signal handler */
  return data.res;              /* return the created result */
}  /* py_ista() */

/*--------------------------------------------------------------------*/
/* apriacc (tracts, supp=-2, zmin=2, zmax=None, report='aP',          */
/*          stat='c', siglvl=1, prune=0, mode='', border=None)        */
/*--------------------------------------------------------------------*/

static PyObject* py_apriacc (PyObject *self,
                             PyObject *args, PyObject *kwds)
{                               /* --- Apriori algorithm */
  char     *ckwds[] = { "tracts", "supp", "zmin", "zmax", "report",
                        "stat", "siglvl", "prune", "mode", "border",
                        NULL };
  double   supp    = -2;        /* minimum support of an item set */
  long     zmin    =  2;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "aP";      /* indicators of values to report */
  CCHAR    *sstat  = "c";       /* test statistic as a string (chi^2) */
  int      stat    =  0;        /* test statistic as an identifier */
  double   siglvl  =  1;        /* minimum evaluation measure value */
  CCHAR    *smode  = "";        /* operation mode/flags as a string */
  int      mode    = APR_DEFAULT;  /* operation mode/flags */
  long     prune   =  0;        /* min. size for evaluation filtering */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  APRIORI  *apriori;            /* apriori miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|dllssdlsO", ckwds,
      &tracts, &supp, &zmin, &zmax, &report,
      &sstat, &siglvl, &prune, &smode, &border))
    return NULL;                /* parse the function arguments */
  if (zmin < 0)    return ERR_VALUE("zmin must be >= 0");
  if (zmax < 0)    zmax = LONG_MAX;  /* check the size range */
  if (zmax < zmin) return ERR_VALUE("zmax must be >= zmin");
  if (zmin > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax > ITEM_MAX) zmax = ITEM_MAX;
  stat = get_stat(sstat);       /* translate the statistic string */
  if (stat < 0)    return NULL; /* and check whether it is valid */
  if (siglvl <= 0) return ERR_VALUE("siglvl must be positive");
  if (strchr(smode, 'z')) stat |= APR_INVBXS;

  /* --- get and prepare transactions --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromPyObj(tracts, NULL);
  if (!tabag) return ERR_RTN(); /* get transactions and create miner */
  apriori = apriori_create(ISR_MAXIMAL, supp, 100.0, 100.0,
                           (ITEM)zmin, (ITEM)zmax,
                           stat, APR_MAX, siglvl, APR_AUTO, mode);
  if (!apriori) { tbg_delete(tabag, 1); return ERR_MEM(); }
  r = apriori_data(apriori, tabag, 0, +2);
  if (r) apriori_delete(apriori, 1);   /* prepare data for apriori */
  if (r == -1) return ERR_MEM();/* check for an error and no items */
  if (r <   0) { sig_remove(); return PyList_New(0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create an item set reporter */
  ||  (apriori_report(apriori, isrep) != 0)) {
    apriori_delete(apriori, 1); return ERR_MEM(); }
  if (isr_pyborder(isrep, border) == NULL) {
    apriori_delete(apriori, 1); return ERR_RTN(); }
  if ((repinit(&data, isrep, report, ISR_SETS) != 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    apriori_delete(apriori, 1); return ERR_MEM(); }

  /* --- frequent item set mining --- */
  if (prune < ITEM_MIN) prune = ITEM_MIN;
  if (prune > ITEM_MAX) prune = ITEM_MAX;
  r = apriori_mine(apriori, (ITEM)prune, 1.0, 0);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  apriori_delete(apriori, 1);   /* delete the apriori miner */
  if (sig_aborted()) { sig_abort(0);
    clean1(data.res); PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (r < 0) { clean1(data.res); return ERR_MEM(); }
  sig_remove();                 /* remove the signal handler */
  return data.res;              /* return the created result */
}  /* py_apriacc() */

/*--------------------------------------------------------------------*/
/* accretion (tracts, supp=-2, zmin=2, zmax=None, report='aP',        */
/*            stat='c', siglvl=1, maxext=2, mode='', border=None)     */
/*--------------------------------------------------------------------*/

static PyObject* py_accretion (PyObject *self,
                               PyObject *args, PyObject *kwds)
{                               /* --- Accretion algorithm */
  char     *ckwds[] = { "tracts", "supp", "zmin", "zmax", "report",
                        "stat", "siglvl", "maxext", "mode", "border",
                        NULL };
  double   supp    =  1;        /* minimum support of an item set */
  long     zmin    =  2;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "aP";      /* indicators of values to report */
  CCHAR    *sstat  = "c";       /* test statistic as a string (chi^2) */
  int      stat    =  0;        /* test statistic as an identifier */
  double   siglvl  =  1;        /* significance level (max. p-value) */
  CCHAR    *smode  = "";        /* operation mode/flags as a string */
  int      mode    = ACC_DEFAULT;  /* operation mode/flags */
  long     maxext  =  2;        /* maximum number of extension items */
  PyObject *border = NULL;      /* support border for filtering */
  PyObject *tracts;             /* transaction database */
  TABAG    *tabag;              /* transaction bag */
  ISREPORT *isrep;              /* item set reporter */
  ACCRET   *accret;             /* accretion miner */
  REPDATA  data;                /* data for item set reporting */
  int      r;                   /* result of function call */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|dllssdlsO", ckwds,
        &tracts, &supp, &zmin, &zmax, &report,
        &sstat, &siglvl, &maxext, &smode, &border))
    return NULL;                /* parse the function arguments */
  if (zmin < 0)    return ERR_VALUE("zmin must be >= 0");
  if (zmax < 0)    zmax = LONG_MAX;  /* check the size range */
  if (zmax < zmin) return ERR_VALUE("zmax must be >= zmin");
  if (zmin > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax > ITEM_MAX) zmax = ITEM_MAX;
  stat = get_stat(sstat);       /* translate the statistic string */
  if (stat < 0)    return NULL; /* and check whether it is valid */
  if (strchr(smode, 'z')) stat |= ACC_INVBXS;
  if (siglvl <= 0) return ERR_VALUE("siglvl must be positive");
  if (maxext < 0)               /* a negative value means that */
    maxext = LONG_MAX;          /* there is no limit on extensions */

  /* --- create transaction bag --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromPyObj(tracts, NULL);
  if (!tabag) return ERR_RTN(); /* get the given transactions */
  accret = accret_create(ISR_MAXIMAL, supp, 100.0,
                         (ITEM)zmin, (ITEM)zmax, stat, siglvl, mode);
  if (!accret) { tbg_delete(tabag, 1); ERR_MEM(); }
  r = accret_data(accret, tabag, +2);
  if (r) accret_delete(accret, 1);  /* prepare data for accretion */
  if (r == -1) return ERR_MEM();/* check for error and no items */
  if (r <   0) { sig_remove(); return PyList_New(0); }

  /* --- create item set reporter --- */
  isrep = isr_create(tbg_base(tabag));
  if (!isrep                    /* create item set reporter */
  ||  (accret_report(accret, isrep) != 0)) {
    accret_delete(accret, 1); return ERR_MEM(); }
  if (isr_pyborder(isrep, border) == NULL) {
    accret_delete(accret, 1); return ERR_RTN(); }
  if ((repinit(&data, isrep, report, ISR_SETS) != 0)
  ||  (isr_setup(isrep) < 0)) { /* set up the item set reporter */
    accret_delete(accret, 1); return ERR_MEM(); }

  /* --- frequent item set mining --- */
  if (maxext > ITEM_MAX) maxext = ITEM_MAX;
  r = accret_mine(accret, (ITEM)maxext);
  if (r >= 0) r = repterm(&data, isrep, report);

  /* --- clean up --- */
  accret_delete(accret, 1);     /* delete the accretion miner */
  if (sig_aborted()) { sig_abort(0);
    clean1(data.res); PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (r < 0) { clean1(data.res); return ERR_MEM(); }
  sig_remove();                 /* remove the signal handler */
  return data.res;              /* return the created result */
} /* py_accretion() */

/*--------------------------------------------------------------------*/
/* genpsp (tracts, target='c', supp=2, zmin=2, zmax=None,             */
/*         report='#', cnt=1000, surr='p', seed=0, cpus=0)            */
/*--------------------------------------------------------------------*/

static PyObject* py_genpsp (PyObject *self,
                            PyObject *args, PyObject *kwds)
{                               /* --- generate a pattern spectrum */
  char     *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                        "report", "cnt", "surr", "seed", "cpus", NULL};
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_CLOSED;/* target type as an identifier */
  double   supp    = -2;        /* minimum support of an item set */
  long     zmin    =  2;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  CCHAR    *report = "#";       /* indicators of reporting format */
  long     cnt     = 1000;      /* number of data sets to generate */
  CCHAR    *ssurr  = "p";       /* surrogate method as a string */
  int      surr    =  2;        /* surrogate method identifier */
  long     seed    =  0;        /* seed for random number generator */
  int      cpus    =  0;        /* number of cpus */
  long     done    =  0;        /* number of completed surrogates */
  PATSPEC  *psp    = NULL;      /* created pattern spectrum (C)*/
  PyObject *pypsp  = NULL;      /* created pattern spectrum (Python) */
  PyObject *tracts;             /* transaction bag (Python) */
  TABAG    *tabag;              /* transaction bag (C) */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sdllslsli", ckwds,
         &tracts, &starg, &supp, &zmin, &zmax, &report,
         &cnt, &ssurr, &seed, &cpus))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "ascm");
  if (target <  0) return NULL; /* translate the target string */
  if (zmin   <  1) return ERR_VALUE("zmin must be positive");
  if (zmax   <  1) zmax = LONG_MAX;  /* check the size range */
  if (zmax   <  zmin) return ERR_VALUE("zmax must be >= zmin");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  if (cnt    <= 0) cnt = 1;     /* check the number of data sets */
  surr   = get_surr(ssurr);     /* translate the surrogate string */
  if (surr   <  0) return NULL; /* and check for a valid code */
  if (surr   == 0) cnt = 1;     /* only one surrogate for identity */
  if (seed   == 0) seed = (long)time(NULL);

  /* --- generate pattern spectrum --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromPyObj(tracts, NULL);
  if (!tabag) return ERR_RTN(); /* get the given transactions */
  if ((surr == FPG_SHUFFLE) && !tbg_istab(tabag)) {
    tbg_delete(tabag, 1);       /* if shuffle surrogates requested */
    return ERR_VALUE("for shuffle surrogates transactions "
                     "must form a table");
  }                             /* check for table-derived data */
  psp  = fpg_genpsp(tabag, target, supp, (ITEM)zmin, (ITEM)zmax,
                    FPG_SIMPLE, FPG_DEFAULT, (size_t)cnt, surr, seed,
                    cpus, repfn, &done);
  if (psp) { pypsp = psp_toPyObj(psp, 1.0/(double)cnt, report[0]);
             psp_delete(psp); } /* generate a pattern spectrum */
  tbg_delete(tabag, 1);         /* delete the created train set */
  if (sig_aborted()) { sig_abort(0);
    clean1(pypsp); PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (!pypsp) return ERR_MEM(); /* check for an error */
  sig_remove();                 /* remove the signal handler */
  return pypsp;                 /* return created pattern spectrum */
}  /* py_genpsp() */

/*--------------------------------------------------------------------*/
/* estpsp (tracts, target='s', supp=2, zmin=2, zmax=None, report='#', */
/*         equiv=10000, alpha=0.5, smpls=1000, seed=0)                */
/*--------------------------------------------------------------------*/

static PyObject* py_estpsp (PyObject *self,
                            PyObject *args, PyObject *kwds)
{                               /* --- estimate a pattern spectrum */
  char    *ckwds[] = { "tracts", "target", "supp", "zmin", "zmax",
                       "report", "equiv", "alpha", "smpls", "seed",
                       NULL };
  long     equiv   = 10000;     /* equivalent number of surrogates */
  CCHAR    *starg  = "s";       /* target type as a string */
  int      target  = ISR_SETS;  /* target type identifier */
  double   supp    =  2;        /* minimum support of an item set */
  double   alpha   =  0.5;      /* probability dispersion factor */
  long     smpls   = 1000;      /* number of samples per set size */
  long     zmin    =  2;        /* minimum size of an item set */
  long     zmax    = -1;        /* maximum size of an item set */
  long     seed    =  0;        /* seed for random number generator */
  CCHAR    *report = "#";       /* indicators of reporting format */
  PyObject *pypsp  = NULL;      /* created Python pattern spectrum */
  PyObject *tracts;             /* transaction database (Python) */
  TABAG    *tabag;              /* transaction database (C) */
  PATSPEC  *psp;                /* created pattern spectrum */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sdllsldll", ckwds,
        &tracts, &starg, &supp, &zmin, &zmax, &report,
        &equiv, &alpha, &smpls, &seed))
    return NULL;                /* parse the function arguments */
  target = get_target(starg, "as");
  if (target <  0) return NULL; /* translate the target string */
  if (zmin   <  1) return ERR_VALUE("zmin must be positive");
  if (zmax   <  1) zmax = LONG_MAX;  /* check the size range */
  if (zmax   < zmin) return ERR_VALUE("zmax must be >= zmin");
  if (zmin   > ITEM_MAX) zmin = ITEM_MAX;
  if (zmax   > ITEM_MAX) zmax = ITEM_MAX;
  if (equiv  <= 0) equiv = 1;   /* check the number of data sets */
  if (smpls  <= 0) return ERR_VALUE("smpls must be positive");
  if (seed   == 0) seed = (long)time(NULL);

  /* --- estimate pattern spectrum --- */
  sig_install();                /* install the signal handler */
  tabag = tbg_fromPyObj(tracts, NULL);
  if (!tabag) return ERR_RTN(); /* get the given transactions */
#if 0
  if (tbg_recode(tabag, smin, -1, -1, -2) < 0) {
    tbg_delete(tabag, 1); return ERR_MEM(); }
  tbg_filter(tabag, (ITEM)zmin, NULL, 0);
#endif
  psp = fpg_estpsp(tabag, target, supp, (ITEM)zmin, (ITEM)zmax,
                   (size_t)equiv, alpha, (size_t)smpls, seed);
  if (psp) { pypsp = psp_toPyObj(psp, 1/(double)equiv, report[0]);
             psp_delete(psp); } /* estimate a pattern spectrum */
  tbg_delete(tabag, 1);         /* and the train set */
  if (sig_aborted()) { sig_abort(0);
    clean1(pypsp); PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (!pypsp) return ERR_MEM(); /* check for an error and return */
  sig_remove();                 /* remove the signal handler */
  return pypsp;                 /* the created pattern spectrum */
}  /* py_estpsp() */

/*--------------------------------------------------------------------*/
/* psp2bdr (psp)                                                      */
/*--------------------------------------------------------------------*/

static PyObject* py_psp2bdr (PyObject *self,
                             PyObject *args, PyObject *kwds)
{                               /* --- extract a decision border */
  char     *ckwds[] = { "psp", NULL };
  PyObject *psp;                /* pattern spectrum */
  PyObject *ei;                 /* pattern spectrum element iterator */
  PyObject *elem;               /* to traverse the elements */
  PyObject *size;               /* to traverse the sizes */
  PyObject *supp;               /* to traverse the support values */
  ITEM     z, zmax;             /* (maximum) pattern size */
  RSUPP    s, *bdr;             /* pattern support, decision border */
  PyObject *pybdr;              /* Python decision border */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", ckwds, &psp))
    return NULL;                /* parse the function arguments */

  /* --- find maximum size --- */
  zmax = 2;                     /* initialize the maximum size */
  ei   = PyObject_GetIter(psp); /* get an iterator for the trains */
  if (!ei) return ERR_TYPE("pattern spectrum must be iterable");
  while ((elem = PyIter_Next(ei)) != NULL) {
    if (!PySequence_Check(elem) /* elements must be sequences */
    ||  (PySequence_Length(elem) < 2)) { clean1(elem);
      return ERR_TYPE("pattern spectrum elements "
                      "must have length >= 2"); }
    size = PySequence_GetItem(elem, 0);
    if (!size) { clean2(elem, ei); return NULL; }
    if      (PyLong_Check(size))/* if size is a long integer */
      z = (ITEM)PyLong_AsLong(size);
    else if (PyInt_Check(size)) /* if size is an integer */
      z = (ITEM)PyInt_AsLong(size);
    else if (PyFloat_Check(size))  /* if size is a float */
      z = (ITEM)PyFloat_AsDouble(size);
    else z = 0;                 /* interpret all else as 0 */
    Py_DECREF(size);            /* drop the size reference */
    Py_DECREF(elem);            /* drop the element reference */
    if (z >= zmax) zmax = z;    /* update the maximum size */
  }
  Py_DECREF(ei);                /* drop the element iterator */

  /* --- find decision border --- */
  bdr = (RSUPP*)calloc((size_t)(zmax+1), sizeof(RSUPP));
  if (!bdr) return ERR_MEM();   /* create a decision border */
  ei = PyObject_GetIter(psp);   /* get an iterator for the trains */
  if (!ei) return ERR_TYPE("pattern spectrum must be iterable");
  while ((elem = PyIter_Next(ei)) != NULL) {
    if (!PySequence_Check(elem) /* elements must be sequences */
    ||  (PySequence_Length(elem) < 2)) { clean1(elem);
      return ERR_TYPE("pattern spectrum elements "
                      "must have length >= 2"); }
    size = PySequence_GetItem(elem, 0);
    if (!size) { clean2(elem, ei); return NULL; }
    if      (PyLong_Check(size))/* if size is a long integer */
      z = (ITEM)PyLong_AsLong(size);
    else if (PyInt_Check(size)) /* if size is an integer */
      z = (ITEM)PyInt_AsLong(size);
    else if (PyFloat_Check(size))  /* if size is a float */
      z = (ITEM)PyFloat_AsDouble(size);
    else z = 0;                 /* interpret all else as 0 */
    Py_DECREF(size);            /* drop the size reference */
    supp = PySequence_GetItem(elem, 1);
    if (!supp) { clean3(size, elem, ei); return NULL; }
    if      (PyLong_Check(size))/* if support is a long integer */
      s = (RSUPP)PyLong_AsLong(supp);
    else if (PyInt_Check(supp)) /* if support is an integer */
      s = (RSUPP)PyInt_AsLong(supp);
    else if (PyFloat_Check(supp))  /* if support is a float */
      s = (RSUPP)PyFloat_AsDouble(supp);
    else s = 0;                 /* interpret all else as 0 */
    Py_DECREF(supp);            /* drop the support reference */
    Py_DECREF(elem);            /* drop the element reference */
    if (s > bdr[z]) bdr[z] = s; /* update the decision border */
  }                             /* (maximum support per size) */

  /* --- build decision border --- */
  for (z = zmax; --z >= 0; )    /* ensure monotone border */
    if (bdr[z+1] > bdr[z]) bdr[z] = bdr[z+1];
  pybdr = PyList_New((Py_ssize_t)(zmax+1));
  if (!pybdr) return ERR_MEM(); /* create a border list */
  supp  = PyFloat_FromDouble(INFINITY);
  PyList_SET_ITEM(pybdr, 0, supp);
  PyList_SET_ITEM(pybdr, 1, supp);
  for (z = 2; z <= zmax; z++) { /* traverse the pattern sizes */
    supp = PyInt_FromLong((long)bdr[z]+1);
    if (!supp) { clean1(pybdr); return ERR_MEM(); }
    PyList_SET_ITEM(pybdr, z, supp);
  }                             /* set minimum support value */
  return pybdr;                 /* return created decision border */
}  /* py_psp2bdr() */

/*--------------------------------------------------------------------*/
/* patred (pats, method='S', border=None, addis=False)                */
/*--------------------------------------------------------------------*/

static PyObject* py_patred (PyObject *self,
                            PyObject *args, PyObject *kwds)
{                               /* --- extract a decision border */
  char     *ckwds[] = { "pats", "method", "border", "addis", NULL };
  PyObject *pats   = NULL;      /* pattern set to reduce */
  CCHAR    *smeth  = "S";       /* pattern set reduction method */
  int      method  = PSR_COVER1;/* pattern set reduction method */
  PyObject *border = NULL;      /* support border for filtering */
  char     addis   = 0;         /* whether to add intersections */
  IDMAP    *map;                /* item to identifier map */
  PATSET   *patset;             /* internal pattern set */
  size_t   i, n;                /* loop variable, number of patterns */
  size_t   k, z;                /* (maximal) size of a pattern */
  size_t   x;                   /* pattern extent (item instances) */
  PyObject *pi;                 /* pattern iterator */
  PyObject *pat;                /* to traverse the patterns */
  PyObject *iset;               /* to traverse the item sets */
  PyObject *ii;                 /* item iterator */
  PyObject *item;               /* to traverse the items */
  PyObject *supp;               /* to traverse the support values */
  PyObject *red;                /* reduced pattern set */
  int      dict;                /* whether pattern set is dictionary */

  /* --- evaluate the function arguments --- */
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|sOb", ckwds,
        &pats, &smeth, &border, &addis))
    return NULL;                /* parse the function arguments */
  method = get_red(smeth);      /* translate the method string */
  if (method < 0) return NULL;  /* into an identifier/offset */

  /* --- find maximum size and item count --- */
  sig_install();                /* install the signal handler */
  n = x = 0; k = 1;             /* init. counters and max. size */
  dict = PyDict_Check(pats);    /* check whether pattern set is dict. */
  pi = PyObject_GetIter(pats);  /* get a pattern iterator */
  if (!pi) return ERR_TYPE("pattern set must be iterable");
  while ((pat = PyIter_Next(pi)) != NULL) {
    if (dict) iset = pat;       /* if dictionary, key is item set */
    else {                      /* if some other iterable, */
      if (!PySequence_Check(pat)/* patterns must be sequences */
      ||  (PySequence_Length(pat) < 2)) { clean2(pat, pi);
        return ERR_TYPE("patterns in non-dictionary must be pairs "
                        "(sequences with length >= 2)"); }
      iset = PySequence_GetItem(pat, 0);
      if (!iset) { clean2(pat, pi); return NULL; }
      Py_DECREF(pat);           /* get the item set and */
    }                           /* drop the pattern reference */
    if      (PySequence_Check(iset))     /* if sequence, */
      z = (size_t)PySequence_Size(iset); /* get its length */
    else if (PyAnySet_Check(iset))       /* if set or frozenset, */
      z = (size_t)PySet_Size(iset);      /* get its size */
    else {                      /* if any other iterable object */
      ii = PyObject_GetIter(iset);  /* check for an iterable */
      if (!ii) { clean2(iset, pi);  /* (list, tuple, set etc.) */
        return ERR_TYPE("item set of a pattern must be iterable"); }
      for (z = 0; (item = PyIter_Next(ii)) != NULL; z++)
        Py_DECREF(item);        /* count the number of items and */
      Py_DECREF(ii);            /* drop the iterartor reference */
      if (PyErr_Occurred()) { clean2(iset, pi); return NULL; }
    }
    Py_DECREF(iset);            /* drop the item set reference */
    if (z > k) k = z;           /* update maximum size and extent */
    x += z; n += 1;             /* and count the pattern */
    if (sig_aborted()) break;   /* check for user abort */
  }
  Py_DECREF(pi);                /* drop the pattern iterator */
  if (PyErr_Occurred()) return NULL;
  if (sig_aborted()) { sig_abort(0);
    PyErr_SetInterrupt(); return ERR_ABORT(); }
  /* n: number of patterns (item set/support pairs) */
  /* k: maximum size of a pattern (number of items) */
  /* x: total number of item istances (extent) */

  /* --- create pattern set --- */
  z   = (n > 255) ? n : 255;    /* compute identifier map size */
  map = idm_create(z, 0, hashitem, cmpitems, NULL, NULL);
  if (!map) return ERR_MEM();   /* create item map and pattern set */
  patset = psr_create(n, (ITEM)k, x, map);
  if (!patset) { idm_delete(map); return ERR_MEM(); }
  if (border) {                 /* if a decision border is given */
    if (!PySequence_Check(border)) { psr_delete(patset, 1);
      return ERR_TYPE("border must be a sequence"); }
    i = (size_t)PySequence_Length(border);
    if (i > k+1) i = k+1;       /* limit border to pattern size */
    for (z = 2; z < i; z++) {   /* traverse the border */
      supp = PySequence_GetItem(border, (Py_ssize_t)z);
      if (!supp) { psr_delete(patset, 1); return NULL; }
      if      (PyLong_Check(supp))  /* if support is a long integer */
        psr_setbdr(patset, (ITEM)z, (RSUPP)PyLong_AsLong(supp));
      else if (PyInt_Check(supp))   /* if support is an integer */
        psr_setbdr(patset, (ITEM)z, (RSUPP)PyInt_AsLong(supp));
      else if (PyFloat_Check(supp)) /* if support is a float */
        psr_setbdr(patset, (ITEM)z, (RSUPP)PyFloat_AsDouble(supp));
      Py_DECREF(supp);          /* drop the support reference */
    }                           /* (fill the decision border */
  }                             /*  with minimum support values) */
  if (sig_aborted()) { sig_abort(0);
    psr_delete(patset, 1); PyErr_SetInterrupt(); return ERR_ABORT(); }

  /* --- collect patterns --- */
  pi = PyObject_GetIter(pats);  /* get a pattern iterator */
  if (!pi) { psr_delete(patset, 1);
    return ERR_TYPE("pattern set must be iterable"); }
  for (i = 0; (pat = PyIter_Next(pi)) != NULL; i++) {
    psr_addorig(patset, pat);   /* note the input pattern */
    if (dict) {                 /* if dictionary, key is item set */
      iset = pat;               /* that is mapped to the support */
      supp = PyDict_GetItem(pats, iset); Py_INCREF(supp); }
    else {                      /* if some other iterable, */
      if (!PySequence_Check(pat)/* patterns must be sequences */
      ||  (PySequence_Length(pat) < 2)) {
        clean2(pat, pi); psr_delete(patset, 1);
        return ERR_TYPE("patterns in non-dictionary must be pairs "
                        "(sequences with length >= 2)"); }
      iset = PySequence_GetItem(pat, 0);
      if (!iset) { clean2(pat, pi); psr_delete(patset,1); return NULL; }
      supp = PySequence_GetItem(pat, 1);
      Py_DECREF(pat);           /* get item set and its support */
      if (!supp) { clean2(iset,pi); psr_delete(patset,1); return NULL; }
    }                           /* drop the pattern reference */
    ii = PyObject_GetIter(iset);/* check for an iterable */
    if (!ii) { clean3(supp, iset, pi); psr_delete(patset, 1);
      return ERR_TYPE("item set of a pattern must be an iterable"); }
    Py_DECREF(iset);            /* drop the item set reference */
    for (z = 0; (item = PyIter_Next(ii)) != NULL; z++) {
      if (psr_additem(patset, &item) != 0) { clean4(item, supp, ii, pi);
        psr_delete(patset, 1); return ERR_MEM(); }
      Py_DECREF(item);          /* add item to current pattern */
    }                           /* and drop the item reference */
    Py_DECREF(ii);              /* drop the item iterator */
    if      (PyInt_Check(supp)) /* get the pattern support */
      psr_addsupp(patset, (RSUPP)PyInt_AsLong(supp));
    else if (PyLong_Check(supp))
      psr_addsupp(patset, (RSUPP)PyLong_AsLong(supp));
    else if (PyFloat_Check(supp))
      psr_addsupp(patset, (RSUPP)PyFloat_AsDouble(supp));
    else { clean2(supp, pi); psr_delete(patset, 1);
      return ERR_TYPE("pattern support must be a number"); }
    Py_DECREF(supp);            /* drop the support reference */
    if (sig_aborted()) break;   /* check for user abort */
  }
  Py_DECREF(pi);                /* drop the pattern iterator */
  if (PyErr_Occurred()) { psr_delete(patset, 1); return NULL; }
  if (sig_aborted()) { sig_abort(0); psr_delete(patset, 1);
    PyErr_SetInterrupt(); return ERR_ABORT(); }

  /* --- pattern set reduction --- */
  k = psr_reduce(patset, method, addis);
  if (sig_aborted()) { sig_abort(0); psr_delete(patset, 1);
    PyErr_SetInterrupt(); return ERR_ABORT(); }
  if (dict) {                   /* if patterns in dictionary */
    red = PyDict_New();         /* create reduced pattern set */
    if (!red) { psr_delete(patset, 1); return ERR_MEM(); }
    for (i = 0; i < n; i++) {   /* traverse the patterns */
      pat = psr_getorig(patset, i);
      if (!pat) continue;       /* skip removed patterns */
      supp = PyDict_GetItem(pats, pat);
      if (!supp) { clean1(red); psr_delete(patset, 1); return NULL; }
      if (PyDict_SetItem(red, pat, supp) != 0) {
        clean1(red); psr_delete(patset, 1); return ERR_MEM(); }
      if (sig_aborted()) break; /* check for user abort */
    } }                         /* (return pattern set in same form) */
  else {                        /* if patterns in other iterable */
    red = PyList_New((Py_ssize_t)k);
    if (!red) { psr_delete(patset, 1); return ERR_MEM(); }
    for (i = k = 0; i < n; i++){/* traverse the patterns */
      pat = psr_getorig(patset, i);
      if (!pat) continue;       /* skip removed patterns */
      Py_INCREF(pat);           /* add pair as taken from pattern */
      PyList_SET_ITEM(red, (Py_ssize_t)k, pat);
      k += 1;                   /* count the output pattern */
      if (sig_aborted()) break; /* check for user abort */
    }
  }
  psr_delete(patset, 1);        /* clean up the pattern set */
  if (sig_aborted()) { sig_abort(0);
    clean1(pat); PyErr_SetInterrupt(); return ERR_ABORT(); }
  sig_remove();                 /* remove the signal handler */
  return red;                   /* return the reduced pattern set */
}  /* py_patred() */

/*--------------------------------------------------------------------*/
/* Python Function List                                               */
/*--------------------------------------------------------------------*/

static PyMethodDef fim_methods[] = {
  { "fim", (PyCFunction)py_fim, METH_VARARGS|METH_KEYWORDS,
    "fim (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "     eval='x', agg='x', thresh=10, border=None)\n"
    "Find frequent item sets (simplified interface).\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    #if FPSUPP
    "        the keys, the values their (numeric) multiplicities.\n"
    #else
    "        the keys, the values their (integer) multiplicities.\n"
    #endif
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a   sets/all   all     frequent item sets\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "        g     gens       generators\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        Q     support of the empty set (total number of transactions)\n"
    "        (     combine values in a tuple (must be first character)\n"
    "        [     combine values in a list  (must be first character)\n"
    "        #     pattern spectrum as a dictionary  (no patterns)\n"
    "        =     pattern spectrum as a list        (no patterns)\n"
    "        |     pattern spectrum as three columns (no patterns)\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient       (+)\n"
    "        c     conf       rule confidence                            (+)\n"
    "        d     confdiff   absolute confidence difference to prior    (+)\n"
    "        l     lift       lift value (confidence divided by prior)   (+)\n"
    "        a     liftdiff   absolute difference of lift value to 1     (+)\n"
    "        q     liftquot   difference of lift quotient to 1           (+)\n"
    "        v     cvct       conviction (inverse lift for negated head) (+)\n"
    "        e     cvctdiff   absolute difference of conviction to 1     (+)\n"
    "        r     cvctquot   difference of conviction quotient to 1     (+)\n"
    "        k     cprob      conditional probability ratio              (+)\n"
    "        j     import     importance (binary log. of prob. ratio)    (+)\n"
    "        z     cert       certainty factor (relative conf. change)   (+)\n"
    "        n     chi2       normalized chi^2 measure                   (+)\n"
    "        p     chi2pval   p-value from (unnormalized) chi^2 measure  (-)\n"
    "        y     yates      normalized chi^2 with Yates' correction    (+)\n"
    "        t     yatespval  p-value from Yates-corrected chi^2 measure (-)\n"
    "        i     info       information difference to prior            (+)\n"
    "        g     infopval   p-value from G statistic/info. difference  (-)\n"
    "        f     fetprob    Fisher's exact test (table probability)    (-)\n"
    "        h     fetchi2    Fisher's exact test (chi^2 measure)        (-)\n"
    "        m     fetinfo    Fisher's exact test (mutual information)   (-)\n"
    "        s     fetsupp    Fisher's exact test (support)              (-)\n"
    "        Measures marked with (+) must meet or exceed the threshold,\n"
    "        measures marked with (-) must not exceed the threshold\n"
    "        in order for the item set to be reported.\n"
    "agg     evaluation measure aggregation mode    (default: x)\n"
    "        x     none       no aggregation (use first value)\n"
    "        m     min        minimum of individual measure values\n"
    "        n     max        maximum of individual measure values\n"
    "        a     avg        average of individual measure values\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns a list of patterns (i.e. tuples with one or more elements),\n"
    "        each consisting of a tuple with a found frequent item set\n"
    "        and the values selected by 'report', which may be combined\n"
    "        into a tuple or list if report[0] is '(' or '[', respectively,\n"
    "        *or* a pattern spectrum as a dictionary mapping pairs\n"
    "        (size, supp) to the corresponding occurrence frequency,\n"
    "        as a list of triplets (size, supp, freq) or as three\n"
    "        columns for sizes, support values and frequencies\n"
  },
  { "arules", (PyCFunction)py_arules, METH_VARARGS|METH_KEYWORDS,
    "arules (tracts, supp=10, conf=80, zmin=1, zmax=None, report='aC',\n"
    "        eval='x', thresh=10, mode='', appear=None)\n"
    "Find association rules (simplified interface).\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    #if FPSUPP
    "        the keys, the values their (numeric) multiplicities.\n"
    #else
    "        the keys, the values their (integer) multiplicities.\n"
    #endif
    "supp    minimum support    of an assoc. rule   (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "conf    minimum confidence of an assoc. rule   (default: 80%)\n"
    "zmin    minimum number of items per rule       (default: 1)\n"
    "zmax    maximum number of items per rule       (default: no limit)\n"
    "report  values to report with a assoc. rule    (default: aC)\n"
    "        a     absolute item set  support (number of transactions)\n"
    "        s     relative item set  support as a fraction\n"
    "        S     relative item set  support as a percentage\n"
    "        b     absolute body set  support (number of transactions)\n"
    "        x     relative body set  support as a fraction\n"
    "        X     relative body set  support as a percentage\n"
    "        h     absolute head item support (number of transactions)\n"
    "        y     relative head item support as a fraction\n"
    "        Y     relative head item support as a percentage\n"
    "        c     rule confidence as a fraction\n"
    "        C     rule confidence as a percentage\n"
    "        l     lift value of a rule (confidence/prior)\n"
    "        L     lift value of a rule as a percentage\n"
    "        e     value of rule evaluation measure\n"
    "        E     value of rule evaluation measure as a percentage\n"
    "        Q     support of the empty set (total number of transactions)\n"
    "        (     combine values in a tuple (must be first character)\n"
    "        [     combine values in a list  (must be first character)\n"
    "        #     pattern spectrum as a dictionary  (no patterns)\n"
    "        =     pattern spectrum as a list        (no patterns)\n"
    "        |     pattern spectrum as three columns (no patterns)\n"
    "eval    measure for rule evaluation            (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient       (+)\n"
    "        c     conf       rule confidence                            (+)\n"
    "        d     confdiff   absolute confidence difference to prior    (+)\n"
    "        l     lift       lift value (confidence divided by prior)   (+)\n"
    "        a     liftdiff   absolute difference of lift value to 1     (+)\n"
    "        q     liftquot   difference of lift quotient to 1           (+)\n"
    "        v     cvct       conviction (inverse lift for negated head) (+)\n"
    "        e     cvctdiff   absolute difference of conviction to 1     (+)\n"
    "        r     cvctquot   difference of conviction quotient to 1     (+)\n"
    "        k     cprob      conditional probability ratio              (+)\n"
    "        j     import     importance (binary log. of prob. ratio)    (+)\n"
    "        z     cert       certainty factor (relative conf. change)   (+)\n"
    "        n     chi2       normalized chi^2 measure                   (+)\n"
    "        p     chi2pval   p-value from (unnormalized) chi^2 measure  (-)\n"
    "        y     yates      normalized chi^2 with Yates' correction    (+)\n"
    "        t     yatespval  p-value from Yates-corrected chi^2 measure (-)\n"
    "        i     info       information difference to prior            (+)\n"
    "        g     infopval   p-value from G statistic/info. difference  (-)\n"
    "        f     fetprob    Fisher's exact test (table probability)    (-)\n"
    "        h     fetchi2    Fisher's exact test (chi^2 measure)        (-)\n"
    "        m     fetinfo    Fisher's exact test (mutual information)   (-)\n"
    "        s     fetsupp    Fisher's exact test (support)              (-)\n"
    "        Measures marked with (+) must meet or exceed the threshold,\n"
    "        measures marked with (-) must not exceed the threshold\n"
    "        in order for the item set to be reported.\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        o     use original rule support definition (body & head)\n"
    "appear  dictionary mapping items to item appearance indicators,\n"
    "        with the key None referring to the default item appearance.\n"
    "        (If None does not occur as a key or no dictionary is given,\n"
    "        the default item appearance indicator is 'both'.)\n"
    "        * item may not appear anywhere in a rule:\n"
    "          '-', 'n', 'none', 'neither', 'ignore'\n"
    "        * item may appear only in rule body/antecedent:\n"
    "          'i', 'in', 'inp', 'input', 'b', 'body',\n"
    "          'a', 'ante', 'antecedent'\n"
    "        * item may appear only in rule head/consequent:\n"
    "          'o', 'out',      'output', 'h', 'head',\n"
    "          'c', 'cons', 'consequent'\n"
    "        * item may appear anywhere in a rule:\n"
    "          'io', 'i&o', 'inout', 'in&out', 'bh', 'b&h', 'both'\n"
    "returns a list of rules (i.e. tuples with two or more elements),\n"
    "        each consisting of a head/consequent item, a tuple with\n"
    "        a body/antecedent item set, and the values selected by\n"
    "        the parameter 'report', which may be combined into a tuple\n"
    "        or a list if report[0] is '(' or '[', respectively."
  },
  { "apriori", (PyCFunction)py_apriori, METH_VARARGS|METH_KEYWORDS,
    "apriori (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "         eval='x', agg='x', thresh=10, prune=None, algo='b', mode='',\n"
    "         border=None)\n"
    "Find frequent item sets with the Apriori algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    #if FPSUPP
    "        the keys, the values their (numeric) multiplicities.\n"
    #else
    "        the keys, the values their (integer) multiplicities.\n"
    #endif
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a   sets/all   all     frequent item sets\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "        g     gens       generators\n"
    "        r     rules      association rules\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "conf    minimum confidence of an assoc. rule   (default: 80%)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        (     combine values in a tuple (must be first character)\n"
    "        [     combine values in a list  (must be first character)\n"
    "        #     pattern spectrum as a dictionary  (no patterns)\n"
    "        =     pattern spectrum as a list        (no patterns)\n"
    "        |     pattern spectrum as three columns (no patterns)\n"
    "        for target 'r' (association rules) also available:\n"
    "        b     absolute body set  support (number of transactions)\n"
    "        x     relative body set  support as a fraction\n"
    "        X     relative body set  support as a percentage\n"
    "        h     absolute head item support (number of transactions)\n"
    "        y     relative head item support as a fraction\n"
    "        Y     relative head item support as a percentage\n"
    "        c     rule confidence as a fraction\n"
    "        C     rule confidence as a percentage\n"
    "        l     lift value of a rule (confidence/prior)\n"
    "        L     lift value of a rule as a percentage\n"
    "        Q     support of the empty set (total number of transactions)\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient       (+)\n"
    "        c     conf       rule confidence                            (+)\n"
    "        d     confdiff   absolute confidence difference to prior    (+)\n"
    "        l     lift       lift value (confidence divided by prior)   (+)\n"
    "        a     liftdiff   absolute difference of lift value to 1     (+)\n"
    "        q     liftquot   difference of lift quotient to 1           (+)\n"
    "        v     cvct       conviction (inverse lift for negated head) (+)\n"
    "        e     cvctdiff   absolute difference of conviction to 1     (+)\n"
    "        r     cvctquot   difference of conviction quotient to 1     (+)\n"
    "        k     cprob      conditional probability ratio              (+)\n"
    "        j     import     importance (binary log. of prob. ratio)    (+)\n"
    "        z     cert       certainty factor (relative conf. change)   (+)\n"
    "        n     chi2       normalized chi^2 measure                   (+)\n"
    "        p     chi2pval   p-value from (unnormalized) chi^2 measure  (-)\n"
    "        y     yates      normalized chi^2 with Yates' correction    (+)\n"
    "        t     yatespval  p-value from Yates-corrected chi^2 measure (-)\n"
    "        i     info       information difference to prior            (+)\n"
    "        g     infopval   p-value from G statistic/info. difference  (-)\n"
    "        f     fetprob    Fisher's exact test (table probability)    (-)\n"
    "        h     fetchi2    Fisher's exact test (chi^2 measure)        (-)\n"
    "        m     fetinfo    Fisher's exact test (mutual information)   (-)\n"
    "        s     fetsupp    Fisher's exact test (support)              (-)\n"
    "        Measures marked with (+) must meet or exceed the threshold,\n"
    "        measures marked with (-) must not exceed the threshold\n"
    "        in order for the item set to be reported.\n"
    "agg     evaluation measure aggregation mode    (default: x)\n"
    "        x     none       no aggregation (use first value)\n"
    "        m     min        minimum of individual measure values\n"
    "        n     max        maximum of individual measure values\n"
    "        a     avg        average of individual measure values\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "prune   min. size for evaluation filtering     (default: no pruning)\n"
    "        = 0   backward filtering       (no subset check)\n"
    "        < 0   weak   forward filtering (one subset  must qualify)\n"
    "        > 0   strong forward filtering (all subsets must qualify)\n"
    "algo    algorithm variant to use               (default: a)\n"
    "        b     basic      standard algorithm (only choice)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        x     do not use perfect extension pruning\n"
    "        t/T   do not organize transactions as a prefix tree\n"
    "        y     a-posteriori pruning of infrequent item sets\n"
    "        z     invalidate evaluation below expected support\n"
    "        o     use original rule support definition (body & head)\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "appear  dictionary mapping items to item appearance indicators,\n"
    "        with the key None referring to the default item appearance.\n"
    "        (If None does not occur as a key or no dictionary is given,\n"
    "        the default item appearance indicator is 'both'.)\n"
    "        This parameter is only used if the target type is rules.\n"
    "        * item may not appear anywhere in a rule:\n"
    "          '-', 'n', 'none', 'neither', 'ignore'\n"
    "        * item may appear only in rule body/antecedent:\n"
    "          'i', 'in', 'inp', 'input', 'b', 'body',\n"
    "          'a', 'ante', 'antecedent'\n"
    "        * item may appear only in rule head/consequent:\n"
    "          'o', 'out',      'output', 'h', 'head',\n"
    "          'c', 'cons', 'consequent'\n"
    "        * item may appear anywhere in a rule:\n"
    "          'io', 'i&o', 'inout', 'in&out', 'bh', 'b&h', 'both'\n"
    "returns if report is not in ['#','=','|']:\n"
    "          if the target is association rules:\n"
    "            a list of rules (i.e. tuples with two or more elements),\n"
    "            each consisting of a head/consequent item, a tuple with\n"
    "            a body/antecedent item set, and the values selected by\n"
    "            the parameter 'report', which may be combined into a\n"
    "            tuple or a list if report[0] is '(' or '[', respectively.\n"
    "          if the target is a type of item sets:\n"
    "            a list of patterns (i.e. tuples with one or more elements),\n"
    "            each consisting of a tuple with a found frequent item set\n"
    "            and the values selected by the parameter 'report', which\n"
    "            may be combined into a tuple or list if report[0] is '('\n"
    "            or '[', respectively\n"
    "        if report in ['#','=','|']:\n"
    "          a pattern spectrum as a dictionary mapping pattern sizes\n"
    "          to the corresponding occurrence support ranges, as a list\n"
    "          of triplets (size, min. support, max. support) or as three\n"
    "          columns for sizes and minimum and maximum support values\n"
  },
  { "eclat", (PyCFunction)py_eclat, METH_VARARGS|METH_KEYWORDS,
    "eclat (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "       eval='x', agg='x', thresh=10, prune=None, algo='a', mode='',\n"
    "       border=None)\n"
    "Find frequent item sets with the Eclat algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    #if FPSUPP
    "        the keys, the values their (numeric) multiplicities.\n"
    #else
    "        the keys, the values their (integer) multiplicities.\n"
    #endif
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a   sets/all   all     frequent item sets\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "        g     gens       generators\n"
    "        r     rules      association rules\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "conf    minimum confidence of an assoc. rule   (default: 80%)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        (     combine values in a tuple (must be first character)\n"
    "        [     combine values in a list  (must be first character)\n"
    "        #     pattern spectrum as a dictionary  (no patterns)\n"
    "        =     pattern spectrum as a list        (no patterns)\n"
    "        |     pattern spectrum as three columns (no patterns)\n"
    "        for target 'r' (association rules) also available:\n"
    "        b     absolute body set  support (number of transactions)\n"
    "        x     relative body set  support as a fraction\n"
    "        X     relative body set  support as a percentage\n"
    "        h     absolute head item support (number of transactions)\n"
    "        y     relative head item support as a fraction\n"
    "        Y     relative head item support as a percentage\n"
    "        c     rule confidence as a fraction\n"
    "        C     rule confidence as a percentage\n"
    "        l     lift value of a rule (confidence/prior)\n"
    "        L     lift value of a rule as a percentage\n"
    "        Q     support of the empty set (total number of transactions)\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient       (+)\n"
    "        c     conf       rule confidence                            (+)\n"
    "        d     confdiff   absolute confidence difference to prior    (+)\n"
    "        l     lift       lift value (confidence divided by prior)   (+)\n"
    "        a     liftdiff   absolute difference of lift value to 1     (+)\n"
    "        q     liftquot   difference of lift quotient to 1           (+)\n"
    "        v     cvct       conviction (inverse lift for negated head) (+)\n"
    "        e     cvctdiff   absolute difference of conviction to 1     (+)\n"
    "        r     cvctquot   difference of conviction quotient to 1     (+)\n"
    "        k     cprob      conditional probability ratio              (+)\n"
    "        j     import     importance (binary log. of prob. ratio)    (+)\n"
    "        z     cert       certainty factor (relative conf. change)   (+)\n"
    "        n     chi2       normalized chi^2 measure                   (+)\n"
    "        p     chi2pval   p-value from (unnormalized) chi^2 measure  (-)\n"
    "        y     yates      normalized chi^2 with Yates' correction    (+)\n"
    "        t     yatespval  p-value from Yates-corrected chi^2 measure (-)\n"
    "        i     info       information difference to prior            (+)\n"
    "        g     infopval   p-value from G statistic/info. difference  (-)\n"
    "        f     fetprob    Fisher's exact test (table probability)    (-)\n"
    "        h     fetchi2    Fisher's exact test (chi^2 measure)        (-)\n"
    "        m     fetinfo    Fisher's exact test (mutual information)   (-)\n"
    "        s     fetsupp    Fisher's exact test (support)              (-)\n"
    "        Measures marked with (+) must meet or exceed the threshold,\n"
    "        measures marked with (-) must not exceed the threshold\n"
    "        in order for the item set to be reported.\n"
    "agg     evaluation measure aggregation mode    (default: x)\n"
    "        x     none       no aggregation (use first value)\n"
    "        m     min        minimum of individual measure values\n"
    "        n     max        maximum of individual measure values\n"
    "        a     avg        average of individual measure values\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "prune   min. size for evaluation filtering     (default: no pruning)\n"
    "        = 0   backward filtering       (no subset check)\n"
    "        < 0   weak   forward filtering (one subset  must qualify)\n"
    "        > 0   strong forward filtering (all subsets must qualify)\n"
    "algo    algorithm variant to use               (default: a)\n"
    "        a     auto       automatic choice based on data properties\n"
    "        e     basic      transaction id lists intersection (basic)\n"
    "        i     tids       transaction id lists intersection (improved)\n"
    "        b     bits       transaction id lists as bit vectors\n"
    "        t     table      item occurrence table (standard)\n"
    "        s     simple     item occurrence table (simplified)\n"
    "        r     ranges     transaction id range lists intersection\n"
    "        o     occdlv     occurrence deliver from transaction lists\n"
    "        d     diff       transaction id difference sets (diffsets)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        l     do not use a 16-items machine\n"
    "        x     do not use perfect extension pruning\n"
    "        i     do not sort items w.r.t. conditional support\n"
    "        u     do not head union tail (hut) pruning (maximal)\n"
    "        y     check extensions for closed/maximal item sets\n"
    "        z     invalidate evaluation below expected support\n"
    "        o     use original rule support definition (body & head)\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "appear  dictionary mapping items to item appearance indicators,\n"
    "        with the key None referring to the default item appearance.\n"
    "        (If None does not occur as a key or no dictionary is given,\n"
    "        the default item appearance indicator is 'both'.)\n"
    "        This parameter is only used if the target type is rules.\n"
    "        * item may not appear anywhere in a rule:\n"
    "          '-', 'n', 'none', 'neither', 'ignore'\n"
    "        * item may appear only in rule body/antecedent:\n"
    "          'i', 'in', 'inp', 'input', 'b', 'body',\n"
    "          'a', 'ante', 'antecedent'\n"
    "        * item may appear only in rule head/consequent:\n"
    "          'o', 'out',      'output', 'h', 'head',\n"
    "          'c', 'cons', 'consequent'\n"
    "        * item may appear anywhere in a rule:\n"
    "          'io', 'i&o', 'inout', 'in&out', 'bh', 'b&h', 'both'\n"
    "returns if report is not in ['#','=','|']:\n"
    "          if the target is association rules:\n"
    "            a list of rules (i.e. tuples with two or more elements),\n"
    "            each consisting of a head/consequent item, a tuple with\n"
    "            a body/antecedent item set, and the values selected by\n"
    "            the parameter 'report', which may be combined into a\n"
    "            tuple or a list if report[0] is '(' or '[', respectively.\n"
    "          if the target is a type of item sets:\n"
    "            a list of patterns (i.e. tuples with one or more elements),\n"
    "            each consisting of a tuple with a found frequent item set\n"
    "            and the values selected by the parameter 'report', which\n"
    "            may be combined into a tuple or list if report[0] is '('\n"
    "            or '[', respectively\n"
    "        if report in ['#','=','|']:\n"
    "          a pattern spectrum as a dictionary mapping pattern sizes\n"
    "          to the corresponding occurrence support ranges, as a list\n"
    "          of triplets (size, min. support, max. support) or as three\n"
    "          columns for sizes and minimum and maximum support values\n"
  },
  { "fpgrowth", (PyCFunction)py_fpgrowth, METH_VARARGS|METH_KEYWORDS,
    "fpgrowth (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "          eval='x', agg='x', thresh=10, prune=Nobe, algo='s', mode='',\n"
    "          border=None)\n"
    "Find frequent item sets with the FP-growth algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    #if FPSUPP
    "        the keys, the values their (numeric) multiplicities.\n"
    #else
    "        the keys, the values their (integer) multiplicities.\n"
    #endif
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a   sets/all   all     frequent item sets\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "        g     gens       generators\n"
    "        r     rules      association rules\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "conf    minimum confidence of an assoc. rule   (default: 80%)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        (     combine values in a tuple (must be first character)\n"
    "        [     combine values in a list  (must be first character)\n"
    "        #     pattern spectrum as a dictionary  (no patterns)\n"
    "        =     pattern spectrum as a list        (no patterns)\n"
    "        |     pattern spectrum as three columns (no patterns)\n"
    "        for target 'r' (association rules) also available:\n"
    "        b     absolute body set  support (number of transactions)\n"
    "        x     relative body set  support as a fraction\n"
    "        X     relative body set  support as a percentage\n"
    "        h     absolute head item support (number of transactions)\n"
    "        y     relative head item support as a fraction\n"
    "        Y     relative head item support as a percentage\n"
    "        c     rule confidence as a fraction\n"
    "        C     rule confidence as a percentage\n"
    "        l     lift value of a rule (confidence/prior)\n"
    "        L     lift value of a rule as a percentage\n"
    "        Q     support of the empty set (total number of transactions)\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient       (+)\n"
    "        c     conf       rule confidence                            (+)\n"
    "        d     confdiff   absolute confidence difference to prior    (+)\n"
    "        l     lift       lift value (confidence divided by prior)   (+)\n"
    "        a     liftdiff   absolute difference of lift value to 1     (+)\n"
    "        q     liftquot   difference of lift quotient to 1           (+)\n"
    "        v     cvct       conviction (inverse lift for negated head) (+)\n"
    "        e     cvctdiff   absolute difference of conviction to 1     (+)\n"
    "        r     cvctquot   difference of conviction quotient to 1     (+)\n"
    "        k     cprob      conditional probability ratio              (+)\n"
    "        j     import     importance (binary log. of prob. ratio)    (+)\n"
    "        z     cert       certainty factor (relative conf. change)   (+)\n"
    "        n     chi2       normalized chi^2 measure                   (+)\n"
    "        p     chi2pval   p-value from (unnormalized) chi^2 measure  (-)\n"
    "        y     yates      normalized chi^2 with Yates' correction    (+)\n"
    "        t     yatespval  p-value from Yates-corrected chi^2 measure (-)\n"
    "        i     info       information difference to prior            (+)\n"
    "        g     infopval   p-value from G statistic/info. difference  (-)\n"
    "        f     fetprob    Fisher's exact test (table probability)    (-)\n"
    "        h     fetchi2    Fisher's exact test (chi^2 measure)        (-)\n"
    "        m     fetinfo    Fisher's exact test (mutual information)   (-)\n"
    "        s     fetsupp    Fisher's exact test (support)              (-)\n"
    "        Measures marked with (+) must meet or exceed the threshold,\n"
    "        measures marked with (-) must not exceed the threshold\n"
    "        in order for the item set to be reported.\n"
    "agg     evaluation measure aggregation mode    (default: x)\n"
    "        x     none       no aggregation (use first value)\n"
    "        m     min        minimum of individual measure values\n"
    "        n     max        maximum of individual measure values\n"
    "        a     avg        average of individual measure values\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "prune   min. size for evaluation filtering     (default: no pruning)\n"
    "        = 0   backward filtering       (no subset check)\n"
    "        < 0   weak   forward filtering (one subset  must qualify)\n"
    "        > 0   strong forward filtering (all subsets must qualify)\n"
    "algo    algorithm variant to use               (default: s)\n"
    "        s     simple     simple  tree nodes (only link and parent)\n"
    "        c     complex    complex tree nodes (children and siblings)\n"
    "        d     single     top-down processing on a single prefix tree\n"
    "        t     topdown    top-down processing of the prefix trees\n"
    "        Variant d does not support closed/maximal item set mining.\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        l     do not use a 16-items machine\n"
    "        x     do not use perfect extension pruning\n"
    "        i     do not sort items w.r.t. conditional support\n"
    "        u     do not head union tail (hut) pruning (maximal)\n"
    "        z     invalidate evaluation below expected support\n"
    "        o     use original rule support definition (body & head)\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "appear  dictionary mapping items to item appearance indicators,\n"
    "        with the key None referring to the default item appearance.\n"
    "        (If None does not occur as a key or no dictionary is given,\n"
    "        the default item appearance indicator is 'both'.)\n"
    "        This parameter is only used if the target type is rules.\n"
    "        * item may not appear anywhere in a rule:\n"
    "          '-', 'n', 'none', 'neither', 'ignore'\n"
    "        * item may appear only in rule body/antecedent:\n"
    "          'i', 'in', 'inp', 'input', 'b', 'body',\n"
    "          'a', 'ante', 'antecedent'\n"
    "        * item may appear only in rule head/consequent:\n"
    "          'o', 'out',      'output', 'h', 'head',\n"
    "          'c', 'cons', 'consequent'\n"
    "        * item may appear anywhere in a rule:\n"
    "          'io', 'i&o', 'inout', 'in&out', 'bh', 'b&h', 'both'\n"
    "returns if report is not in ['#','=','|']:\n"
    "          if the target is association rules:\n"
    "            a list of rules (i.e. tuples with two or more elements),\n"
    "            each consisting of a head/consequent item, a tuple with\n"
    "            a body/antecedent item set, and the values selected by\n"
    "            the parameter 'report', which may be combined into a\n"
    "            tuple or a list if report[0] is '(' or '[', respectively.\n"
    "          if the target is a type of item sets:\n"
    "            a list of patterns (i.e. tuples with one or more elements),\n"
    "            each consisting of a tuple with a found frequent item set\n"
    "            and the values selected by the parameter 'report', which\n"
    "            may be combined into a tuple or list if report[0] is '('\n"
    "            or '[', respectively\n"
    "        if report in ['#','=','|']:\n"
    "          a pattern spectrum as a dictionary mapping pattern sizes\n"
    "          to the corresponding occurrence support ranges, as a list\n"
    "          of triplets (size, min. support, max. support) or as three\n"
    "          columns for sizes and minimum and maximum support values\n"
  },
  { "sam", (PyCFunction)py_sam, METH_VARARGS|METH_KEYWORDS,
    "sam (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "     eval='x', thresh=10, algo='b', mode='', border=None)\n"
    "Find frequent item sets with the SaM algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    #if FPSUPP
    "        the keys, the values their (numeric) multiplicities.\n"
    #else
    "        the keys, the values their (integer) multiplicities.\n"
    #endif
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a   sets/all   all     frequent item sets\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        Q     support of the empty set (total number of transactions)\n"
    "        (     combine values in a tuple (must be first character)\n"
    "        [     combine values in a list  (must be first character)\n"
    "        #     pattern spectrum as a dictionary  (no patterns)\n"
    "        =     pattern spectrum as a list        (no patterns)\n"
    "        |     pattern spectrum as three columns (no patterns)\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "algo    algorithm variant to use               (default: o)\n"
    "        s     simple     basic split and merge algorithm\n"
    "        b     bsearch    split and merge with binary search\n"
    "        d     double     SaM with double source buffering\n"
    "        t     tree       SaM with transaction prefix tree\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        l     do not use a 16-items machine\n"
    "        x     do not use perfect extension pruning\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns if report is not in ['#','=','|']:\n"
    "          a list of patterns (i.e. tuples with one or more elements),\n"
    "          each consisting of a tuple with a found frequent item set\n"
    "          and the values selected by the parameter 'report', which\n"
    "          may be combined into a tuple or list if report[0] is '('\n"
    "          or '[', respectively\n"
    "        if report in ['#','=','|']:\n"
    "          a pattern spectrum as a dictionary mapping pattern sizes\n"
    "          to the corresponding occurrence support ranges, as a list\n"
    "          of triplets (size, min. support, max. support) or as three\n"
    "          columns for sizes and minimum and maximum support values\n"
  },
  { "relim", (PyCFunction)py_relim, METH_VARARGS|METH_KEYWORDS,
    "relim (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "       eval='x', thresh=10, algo='s', mode='', border=None)\n"
    "Find frequent item sets with the RElim algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    #if FPSUPP
    "        the keys, the values their (numeric) multiplicities.\n"
    #else
    "        the keys, the values their (integer) multiplicities.\n"
    #endif
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a   sets/all   all     frequent item sets\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        Q     support of the empty set (total number of transactions)\n"
    "        (     combine values in a tuple (must be first character)\n"
    "        [     combine values in a list  (must be first character)\n"
    "        #     pattern spectrum as a dictionary  (no patterns)\n"
    "        =     pattern spectrum as a list        (no patterns)\n"
    "        |     pattern spectrum as three columns (no patterns)\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "algo    algorithm variant to use               (default: o)\n"
    "        s     simple     basic recursive elimination algorithm\n"
    "        (this parameter is essentially a placeholder for extensions)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        l     do not use a 16-items machine\n"
    "        x     do not use perfect extension pruning\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns if report is not in ['#','=','|']:\n"
    "          a list of patterns (i.e. tuples with one or more elements),\n"
    "          each consisting of a tuple with a found frequent item set\n"
    "          and the values selected by the parameter 'report', which\n"
    "          may be combined into a tuple or list if report[0] is '('\n"
    "          or '[', respectively\n"
    "        if report in ['#','=','|']:\n"
    "          a pattern spectrum as a dictionary mapping pattern sizes\n"
    "          to the corresponding occurrence support ranges, as a list\n"
    "          of triplets (size, min. support, max. support) or as three\n"
    "          columns for sizes and minimum and maximum support values\n"
  },
  { "carpenter", (PyCFunction)py_carpenter, METH_VARARGS|METH_KEYWORDS,
    "carpenter (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "           eval='x', thresh=10, algo='a', mode='', border=None)\n"
    "Find frequent item sets with the Carpenter algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    #if FPSUPP
    "        the keys, the values their (numeric) multiplicities.\n"
    #else
    "        the keys, the values their (integer) multiplicities.\n"
    #endif
    "target  type of frequent item sets to find     (default: s)\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        Q     support of the empty set (total number of transactions)\n"
    "        (     combine values in a tuple (must be first character)\n"
    "        [     combine values in a list  (must be first character)\n"
    "        #     pattern spectrum as a dictionary  (no patterns)\n"
    "        =     pattern spectrum as a list        (no patterns)\n"
    "        |     pattern spectrum as three columns (no patterns)\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "algo    algorithm variant to use               (default: s)\n"
    "        a     auto       automatic choice based on table size\n"
    "        t     table      item occurrence counter table\n"
    "        l     tidlist    transaction identifier lists\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        x     do not use perfect extension pruning\n"
    "        z     filter maximal item sets with repository\n"
    "        y     add only maximal item sets to repository\n"
    "        p     do not collate equal transactions\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns if report is not in ['#','=','|']:\n"
    "          a list of patterns (i.e. tuples with one or more elements),\n"
    "          each consisting of a tuple with a found frequent item set\n"
    "          and the values selected by the parameter 'report', which\n"
    "          may be combined into a tuple or list if report[0] is '('\n"
    "          or '[', respectively\n"
    "        if report in ['#','=','|']:\n"
    "          a pattern spectrum as a dictionary mapping pattern sizes\n"
    "          to the corresponding occurrence support ranges, as a list\n"
    "          of triplets (size, min. support, max. support) or as three\n"
    "          columns for sizes and minimum and maximum support values\n"
  },
  { "ista", (PyCFunction)py_ista, METH_VARARGS|METH_KEYWORDS,
    "ista (tracts, target='s', supp=10, zmin=1, zmax=None, report='a',\n"
    "      eval='x', thresh=10, algo='x', mode='', border=None)\n"
    "Find frequent item sets with the IsTa algorithm.\n"
    "tracts  transaction database to mine (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    #if FPSUPP
    "        the keys, the values their (numeric) multiplicities.\n"
    #else
    "        the keys, the values their (integer) multiplicities.\n"
    #endif
    "target  type of frequent item sets to find     (default: s)\n"
    "        c     closed     closed  frequent item sets\n"
    "        m     maximal    maximal frequent item sets\n"
    "supp    minimum support of an item set         (default: 10)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 1)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: a)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        e     value of item set evaluation measure\n"
    "        E     value of item set evaluation measure as a percentage\n"
    "        Q     support of the empty set (total number of transactions)\n"
    "        (     combine values in a tuple (must be first character)\n"
    "        [     combine values in a list  (must be first character)\n"
    "        #     pattern spectrum as a dictionary  (no patterns)\n"
    "        =     pattern spectrum as a list        (no patterns)\n"
    "        |     pattern spectrum as three columns (no patterns)\n"
    "eval    measure for item set evaluation        (default: x)\n"
    "        x     none       no measure / zero (default)\n"
    "        b     ldratio    binary logarithm of support quotient\n"
    "thresh  threshold for evaluation measure       (default: 10%)\n"
    "algo    algorithm variant to use               (default: x)\n"
    "        x     prefix     use a standard prefix tree\n"
    "        p     patricia   use a patricia tree\n"
    "        (a patricia tree may be faster for very many items\n"
    "        and very few transactions)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        p     do not prune the prefix/patricia tree\n"
    "        z     filter maximal item sets with repository\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns if report is not in ['#','=','|']:\n"
    "          a list of patterns (i.e. tuples with one or more elements),\n"
    "          each consisting of a tuple with a found frequent item set\n"
    "          and the values selected by the parameter 'report', which\n"
    "          may be combined into a tuple or list if report[0] is '('\n"
    "          or '[', respectively\n"
    "        if report in ['#','=','|']:\n"
    "          a pattern spectrum as a dictionary mapping pattern sizes\n"
    "          to the corresponding occurrence support ranges, as a list\n"
    "          of triplets (size, min. support, max. support) or as three\n"
    "          columns for sizes and minimum and maximum support values\n"
  },
  { "apriacc", (PyCFunction)py_apriacc, METH_VARARGS|METH_KEYWORDS,
    "apriacc (tracts, supp=-2, zmin=2, zmax=None, report='aP',\n"
    "         stat='c', siglvl=1, prune=0, mode='', border=None)\n"
    "Find frequent item sets with an accretion-style apriori algorithm\n"
    "(that is, with an interface analog to accretion()).\n"
    "tracts  transaction database to mine           (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    #if FPSUPP
    "        the keys, the values their (numeric) multiplicities.\n"
    #else
    "        the keys, the values their (integer) multiplicities.\n"
    #endif
    "supp    minimum support of an item set         (default: -2)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 2)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: aP)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        p     p-value of item set test as a fraction\n"
    "        P     p-value of item set test as a percentage\n"
    "        Q     support of the empty set (total number of transactions)\n"
    "        (     combine values in a tuple (must be first character)\n"
    "        [     combine values in a list  (must be first character)\n"
    "        #     pattern spectrum as a dictionary  (no patterns)\n"
    "        =     pattern spectrum as a list        (no patterns)\n"
    "        |     pattern spectrum as three columns (no patterns)\n"
    "stat    test statistic for item set evaluation (default: c)\n"
    "        x     none     no statistic / zero\n"
    "        c/n/p chi2     chi^2 measure (default)\n"
    "        y/t   yates    chi^2 measure with Yates' correction\n"
    "        i/g   info     mutual information / G statistic\n"
    "        f     fetprob  Fisher's exact test (table probability)\n"
    "        h     fetchi2  Fisher's exact test (chi^2 measure)\n"
    "        m     fetinfo  Fisher's exact test (mutual information)\n"
    "        s     fetsupp  Fisher's exact test (support)\n"
    "siglvl  significance level (maximum p-value)   (default: 1%)\n"
    "prune   min. size for evaluation filtering     (default: 0)\n"
    "        = 0   backward filtering       (no subset checks)\n"
    "        < 0   weak   forward filtering (one subset  must qualify)\n"
    "        > 0   strong forward filtering (all subsets must qualify)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        z     invalidate evaluation below expected support\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns if report is not in ['#','=','|']:\n"
    "          a list of patterns (i.e. tuples with one or more elements),\n"
    "          each consisting of a tuple with a found frequent item set\n"
    "          and the values selected by the parameter 'report', which\n"
    "          may be combined into a tuple or list if report[0] is '('\n"
    "          or '[', respectively\n"
    "        if report in ['#','=','|']:\n"
    "          a pattern spectrum as a dictionary mapping pattern sizes\n"
    "          to the corresponding occurrence support ranges, as a list\n"
    "          of triplets (size, min. support, max. support) or as three\n"
    "          columns for sizes and minimum and maximum support values\n"
  },
  { "accretion", (PyCFunction)py_accretion, METH_VARARGS|METH_KEYWORDS,
    "accretion (tracts, supp=-2, zmin=2, zmax=None, report='aP',\n"
    "           stat='c', siglvl=1, maxext=2, mode='', border=None)\n"
    "Find frequent item sets with the accretion algorithm.\n"
    "tracts  transaction database to mine           (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    #if FPSUPP
    "        the keys, the values their (numeric) multiplicities.\n"
    #else
    "        the keys, the values their (integer) multiplicities.\n"
    #endif
    "supp    minimum support of an item set         (default: -2)\n"
    "        (positive: percentage, negative: absolute number)\n"
    "zmin    minimum number of items per item set   (default: 2)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  values to report with an item set      (default: aP)\n"
    "        a     absolute item set support (number of transactions)\n"
    "        s     relative item set support as a fraction\n"
    "        S     relative item set support as a percentage\n"
    "        p     p-value of item set test as a fraction\n"
    "        P     p-value of item set test as a percentage\n"
    "        Q     support of the empty set (total number of transactions)\n"
    "        (     combine values in a tuple (must be first character)\n"
    "        [     combine values in a list  (must be first character)\n"
    "        #     pattern spectrum as a dictionary  (no patterns)\n"
    "        =     pattern spectrum as a list        (no patterns)\n"
    "        |     pattern spectrum as three columns (no patterns)\n"
    "stat    test statistic for item set evaluation (default: c)\n"
    "        x     none     no statistic / zero\n"
    "        c/p/n chi2     chi^2 measure (default)\n"
    "        y/t   yates    chi^2 measure with Yates' correction\n"
    "        i/g   info     mutual information / G statistic\n"
    "        f     fetprob  Fisher's exact test (table probability)\n"
    "        h     fetchi2  Fisher's exact test (chi^2 measure)\n"
    "        m     fetinfo  Fisher's exact test (mutual information)\n"
    "        s     fetsupp  Fisher's exact test (support)\n"
    "siglvl  significance level (maximum p-value)   (default: 1%)\n"
    "maxext  maximum number of extension items      (default: 2)\n"
    "mode    operation mode indicators/flags        (default: None)\n"
    "        z     invalidate evaluation below expected support\n"
    "border  support border for filtering item sets (default: None)\n"
    "        Must be a list or tuple of (absolute) minimum support values\n"
    "        per item set size (by which the list/tuple is indexed).\n"
    "returns if report is not in ['#','=','|']:\n"
    "          a list of patterns (i.e. tuples with one or more elements),\n"
    "          each consisting of a tuple with a found frequent item set\n"
    "          and the values selected by the parameter 'report', which\n"
    "          may be combined into a tuple or list if report[0] is '('\n"
    "          or '[', respectively\n"
    "        if report in ['#','=','|']:\n"
    "          a pattern spectrum as a dictionary mapping pattern sizes\n"
    "          to the corresponding occurrence support ranges, as a list\n"
    "          of triplets (size, min. support, max. support) or as three\n"
    "          columns for sizes and minimum and maximum support values\n"
  },
  { "genpsp", (PyCFunction)py_genpsp, METH_VARARGS|METH_KEYWORDS,
    "genpsp (tracts, target='s', supp=2, zmin=2, zmax=None,\n"
    "        report='#', cnt=1000, surr='p', seed=0, cpus=0)\n"
    "Generate a pattern spectrum from surrogate data sets.\n"
    "tracts  transaction database to mine           (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    #if FPSUPP
    "        the keys, the values their (numeric) multiplicities.\n"
    #else
    "        the keys, the values their (integer) multiplicities.\n"
    #endif
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a  sets/all   all     frequent item sets\n"
    "        c    closed     closed  frequent item sets\n"
    "        m    maximal    maximal frequent item sets\n"
    "supp    minimum support of an item set         (default: 2)\n"
    "zmin    minimum number of items per item set   (default: 2)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  pattern spectrum reporting format      (default: #)\n"
    "        #    pattern spectrum as a dictionary (size,supp) -> freq\n"
    "        =    pattern spectrum as a list of triplets\n"
    "        |    pattern spectrum as three columns\n"
    "cnt     number of surrogate data sets          (default: 1000)\n"
    "surr    surrogate data generation method       (default: p)\n"
    "        i    ident      identity (keep original data)\n"
    "        r    random     random transaction generation\n"
    "        p    swap       permutation by pair swaps\n"
    "        s    shuffle    shuffle table-derived data (columns)\n"
    "seed    seed for random number generator       (default: 0)\n"
    "        (seed = 0: use system time as a seed)\n"
    "cpus    number of cpus to use                  (default: 0)\n"
    "        A value <= 0 means all cpus reported as available.\n"
    #ifdef FPSUPP
    "returns a pattern spectrum as a dictionary mapping pattern sizes\n"
    "        to the corresponding occurrence support ranges, as a list\n"
    "        of triplets (size, min. support, max. support) or as three\n"
    "        columns for sizes and minimum and maximum support values\n"
    #else
    "returns a pattern spectrum as a dictionary mapping pairs\n"
    "        (size, supp) to the corresponding occurrence frequency,\n"
    "        as a list of triplets (size, supp, freq) or as three\n"
    "        columns for sizes, support values and frequencies\n"
    #endif
  },
  { "estpsp", (PyCFunction)py_estpsp, METH_VARARGS|METH_KEYWORDS,
    "estpsp (tracts, target='s', supp=2, zmin=2, zmax=None,\n"
    "        report='#', equiv=10000, alpha=0.5, smpls=1000, seed=0)\n"
    "Estimate a pattern spectrum from data characteristics.\n"
    "tracts  transaction database to mine           (mandatory)\n"
    "        The database must be an iterable of transactions;\n"
    "        each transaction must be an iterable of items;\n"
    "        each item must be a hashable object.\n"
    "        If the database is a dictionary, the transactions are\n"
    #if FPSUPP
    "        the keys, the values their (numeric) multiplicities.\n"
    #else
    "        the keys, the values their (integer) multiplicities.\n"
    #endif
    "target  type of frequent item sets to find     (default: s)\n"
    "        s/a  sets/all   all     frequent item sets\n"
    "supp    minimum support of an item set         (default: 2)\n"
    "zmin    minimum number of items per item set   (default: 2)\n"
    "zmax    maximum number of items per item set   (default: no limit)\n"
    "report  pattern spectrum reporting format      (default: #)\n"
    "        #    pattern spectrum as a dictionary\n"
    "        =    pattern spectrum as a list of triplets\n"
    "        |    pattern spectrum as three columns\n"
    "equiv   equivalent number of surrogates        (default: 10000)\n"
    "alpha   probability dispersion factor          (default: 0.5)\n"
    "smpls   number of samples per item set size    (default: 1000)\n"
    "seed    seed for random number generator       (default: 0)\n"
    "        (seed = 0: use system time as a seed)\n"
    #ifdef FPSUPP
    "returns a pattern spectrum as a dictionary mapping pattern sizes\n"
    "        to the corresponding occurrence support ranges, as a list\n"
    "        of triplets (size, min. support, max. support) or as three\n"
    "        columns for sizes and minimum and maximum support values\n"
    #else
    "returns a pattern spectrum as a dictionary mapping pairs\n"
    "        (size, supp) to the corresponding occurrence frequency,\n"
    "        as a list of triplets (size, supp, freq) or as three\n"
    "        columns for sizes, support values and frequencies\n"
    #endif
  },
  { "psp2bdr", (PyCFunction)py_psp2bdr, METH_VARARGS|METH_KEYWORDS,
    "psp2bdr (psp)\n"
    "Extract a decision border from a pattern spectrum.\n"
    "psp     a pattern spectrum as returned by one of the functions\n"
    "        genpsp() or estpsp() in the form of dictionary or a list\n"
    "        of triplets (the three column format is not supported;\n"
    "        for a 3-column pattern spectrum, use psp2bdr(zip(*psp)))\n"
    "returns a decision border as a list of minimum support values\n"
    "        that is to be index by the pattern size (starting at 0)\n"
  },
  { "patred", (PyCFunction)py_patred, METH_VARARGS|METH_KEYWORDS,
    "patred (pats, method='S', border=None, addis=False)\n"
    "Perform pattern set reduction based on a preference relation.\n"
    "pats    set of patterns to reduced             (mandatory)\n"
    "        The pattern set must be an iterable (list, tuple etc.).\n"
    "        If it is a dictionary, the keys must be item sets\n"
    "        that are mapped to their support (number).\n"
    "        Otherwise each iterable element must be a pair, consisting\n"
    "        of a item set as an iterable of item and a support (number).\n"
    "        Each item identifier must be a hashable object.\n"
    "method  index or identifier of pattern comparison method or\n"
    "        a pattern comparison function patcmp(zA, cA, zB, cB, border)\n"
    "border  detection border as a list of minimum support values per size\n"
    "addis   whether to add intersections of patterns\n"
    "returns a reduced set of patterns, either as a dictionary that maps\n"
    "        item sets to support values (if the input was a dictionary)\n"
    "        or as a list of patterns in the form in which they were\n"
    "        passed into this function (all parts are maintained)"
  },
  { NULL }                      /* sentinel */
};

/*----------------------------------------------------------------------
  Initialization Function
----------------------------------------------------------------------*/
#define FIM_DESC \
  "Frequent Item Set Mining and Association Rule Induction for Python\n" \
  "version 6.28 (2017.03.24)     (c) 2011-2017   Christian Borgelt"

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef fimdef = {
  PyModuleDef_HEAD_INIT, "fim", FIM_DESC,
  -1, fim_methods, NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC PyInit_fim (void)
{ return PyModule_Create(&fimdef); }

#else

PyMODINIT_FUNC initfim (void)
{ Py_InitModule3("fim", fim_methods, FIM_DESC); }

#endif
