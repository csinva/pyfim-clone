/*----------------------------------------------------------------------
  File    : patred.c
  Contents: pattern set reduction
  Author  : Christian Borgelt
  History : 2015.08.19 file created
            2015.08.19 bugs in function psr_reduce() fixed
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "arrays.h"
#include "patred.h"
#ifdef PSR_ABORT
#include "sigint.h"
#endif

/*--------------------------------------------------------------------*/
/* Type Definitions                                                   */
/*--------------------------------------------------------------------*/
typedef int PATCMPFN (const FRQPAT *A,
                      const FRQPAT *B, RSUPP *border);

/*--------------------------------------------------------------------*/
/* Auxiliary Functions                                                */
/*--------------------------------------------------------------------*/

static int patcmp (const void *A, const void *B, void *data)
{                               /* --- compare two patterns */
  ITEM i, n;                    /* loop variables */
  ITEM *a, *b;                  /* items of the two patterns */

  n = ((const FRQPAT*)A)->size; /* compare the pattern sizes */
  if (n > ((const FRQPAT*)B)->size) return +1;
  if (n < ((const FRQPAT*)B)->size) return -1;
  a = ((const FRQPAT*)A)->items;
  b = ((const FRQPAT*)B)->items;
  for (i = 0; i < n; i++)       /* traverse and compare the items */
    if (a[i] != b[i]) return (a[i] > b[i]) ? +1 : -1;
  return 0;                     /* if patterns are equal, return 0 */
}  /* patcmp() */

/*--------------------------------------------------------------------*/

static int isect (const FRQPAT *A, const FRQPAT *B, FRQPAT *buf)
{                               /* --- intersect two patterns */
  ITEM iA, iB;                  /* loop variables */

  assert(A && B);               /* check the function arguments */
  for (iA = iB = buf->size = 0; (iA < A->size) && (iB < B->size); ) {
    if      (A->items[iA] < B->items[iB]) iA++;
    else if (A->items[iA] > B->items[iB]) iB++;
    else { buf->items[buf->size++] = A->items[iA++]; iB++; }
  }                             /* collect items that occur in both */
  buf->supp = (A->supp > B->supp) ? A->supp : B->supp;
  return buf->size;             /* set pattern support and size */
}  /* isect() */

/*--------------------------------------------------------------------*/

static int subset (const FRQPAT *A, const FRQPAT *B)
{                               /* --- test for A subset B */
  ITEM iA, iB;                  /* loop variables */

  assert(A && B);               /* check the function arguments */
  if (A->size >= B->size) return 0;  /* check size relationship */
  for (iA = iB = 0; (iA < A->size) && (iB < B->size); ) {
    if (A->items[iA] < B->items[iB]) return 0;
    if (A->items[iA] > B->items[iB]) iB++;
    else                     { iA++; iB++; }
  }                             /* check whether all items in A */
  return (iA < A->size) ? 0 : -1;            /* also occur in B */
}  /* subset() */

/*--------------------------------------------------------------------*/
/* Pattern Preference Functions                                       */
/*--------------------------------------------------------------------*/

static int psr_none (const FRQPAT *A,
                     const FRQPAT *B, RSUPP *border)
{ return 0; }                   /* --- no preference */

/*--------------------------------------------------------------------*/

static int psr_coins0 (const FRQPAT *A,
                       const FRQPAT *B, RSUPP *border)
{                               /* --- excess coincidences 0 */
  if (A->supp >= B->supp) return +1;
  return (B->supp -A->supp < border[B->size]) ? +1 : -1;
}  /* psr_coin0() */

/*--------------------------------------------------------------------*/

static int psr_coins1 (const FRQPAT *A,
                       const FRQPAT *B, RSUPP *border)
{                               /* --- excess coincidences 1 */
  if (A->supp >= B->supp) return +1;
  return (B->supp -A->supp +1 < border[B->size]) ? +1 : -1;
}  /* psr_coin1() */

/*--------------------------------------------------------------------*/

static int psr_items2 (const FRQPAT *A,
                       const FRQPAT *B, RSUPP *border)
{                               /* --- excess items */
  if (A->supp >= B->supp) return +1;
  return (A->supp < border[A->size -B->size +2]) ? -1 : +1;
}  /* psr_items2() */

/*--------------------------------------------------------------------*/

static int psr_cover0 (const FRQPAT *A,
                       const FRQPAT *B, RSUPP *border)
{                               /* --- covered events/spikes 0 */
  if (A->supp >= B->supp) return +1;
  return (A->size *A->supp >= B->size *B->supp) ? +1 : -1;
}  /* psr_cover0() */

/*--------------------------------------------------------------------*/

static int psr_cover1 (const FRQPAT *A,
                       const FRQPAT *B, RSUPP *border)
{                               /* --- covered events/spikes 1 */
  if (A->supp >= B->supp) return +1;
  return ((A->size-1) *A->supp >= (B->size-1) *B->supp) ? +1 : -1;
}  /* psr_cover1() */

/*--------------------------------------------------------------------*/

static int psr_leni0 (const FRQPAT *A,
                      const FRQPAT *B, RSUPP *border)
{                               /* --- lenient combination 0 */
  int xA, xB;                   /* excess items/coins. explainable */

  if (A->supp >= B->supp) return +1;
  xA = (A->supp < border[A->size -B->size +2]);
  xB = (B->supp -A->supp +1 < border[B->size]);
  if ( xA && !xB) return -1;
  if (!xA &&  xB) return +1;
  if (!xA && !xB) return  0;
  return (A->size *A->supp >= B->size *B->supp) ? +1 : -1;
}  /* psr_leni0() */

/*--------------------------------------------------------------------*/

static int psr_leni1 (const FRQPAT *A,
                      const FRQPAT *B, RSUPP *border)
{                               /* --- lenient combination 1 */
  int xA, xB;                   /* excess items/coins. explainable */

  if (A->supp >= B->supp) return +1;
  xA = (A->supp < border[A->size -B->size +2]);
  xB = (B->supp -A->supp +1 < border[B->size]);
  if ( xA && !xB) return -1;
  if (!xA &&  xB) return +1;
  if (!xA && !xB) return  0;
  return ((A->size-1) *A->supp >= (B->size-1) *B->supp) ? +1 : -1;
}  /* psr_leni1() */

/*--------------------------------------------------------------------*/

static int psr_strict0 (const FRQPAT *A,
                        const FRQPAT *B, RSUPP *border)
{                               /* --- strict combination 0 */
  int xA, xB;                   /* excess items/coins. explainable */

  if (A->supp >= B->supp) return +1;
  xA = (A->supp < border[A->size -B->size +2]);
  xB = (B->supp -A->supp +1 < border[B->size]);
  if ( xA && !xB) return -1;
  if (!xA &&  xB) return +1;
  return (A->size *A->supp >= B->size *B->supp) ? +1 : -1;
}  /* psr_strict0() */

/*--------------------------------------------------------------------*/

static int psr_strict1 (const FRQPAT *A,
                        const FRQPAT *B, RSUPP *border)
{                               /* --- strict combination 1 */
  int xA, xB;                   /* excess items/coins. explainable */

  if (A->supp >= B->supp) return +1;
  xA = (A->supp < border[A->size -B->size +2]);
  xB = (B->supp -A->supp +1 < border[B->size]);
  if ( xA && !xB) return -1;
  if (!xA &&  xB) return +1;
  return ((A->size-1) *A->supp >= (B->size-1) *B->supp) ? +1 : -1;
}  /* psr_strict1() */

/*--------------------------------------------------------------------*/

static PATCMPFN *psr_tab[] = {  /* pattern preference functions */
  /* PSR_NONE      0 */  psr_none,
  /* PSR_COINS0    1 */  psr_coins0,
  /* PSR_COINS1    2 */  psr_coins1,
  /* PSR_ITEMS     3 */  psr_items2,
  /* PSR_COVER0    4 */  psr_cover0,
  /* PSR_COVER1    5 */  psr_cover1,
  /* PSR_LENIENT0  6 */  psr_leni0,
  /* PSR_LENIENT1  7 */  psr_leni1,
  /* PSR_STRICT0   8 */  psr_strict0,
  /* PSR_STRICT1   9 */  psr_strict1,
};

/*----------------------------------------------------------------------
  Pattern Set Reduction Functions
----------------------------------------------------------------------*/

PATSET* psr_create (size_t patcnt, ITEM patmax,
                    size_t extent, IDMAP *map)
{                               /* --- create pattern set red. object */
  PATSET *psr;                  /* created pattern set red. object */
  FRQPAT *b;                    /* to access the pattern buffer */

  psr = (PATSET*)malloc(sizeof(PATSET) +(patcnt-1) *sizeof(FRQPAT));
  if (!psr) return NULL;        /* create a pattern set */
  psr->map = map;               /* and note the identifier map */
  if (patmax < 2) patmax = 2;   /* create a decision border */
  psr->border = (RSUPP*)calloc(((size_t)patmax+1), sizeof(RSUPP));
  if (!psr->border) { free(psr); return NULL; }
  psr->border[0] = psr->border[1] = RSUPP_MAX;
  psr->max = patmax;            /* store the size variables */
  psr->cnt = patcnt; psr->cur = 0;
  psr->rem = extent;
  b = &psr->buf;                /* initialize the pattern buffer */
  b->size  = 0; b->supp  = 0; b->orig  = NULL;
  b->items = (ITEM*)malloc((extent+(size_t)patmax)*sizeof(ITEM));
  psr->next = (extent > 0) ? b->items +patmax : NULL;
  return psr;                   /* return the created pattern set */
}  /* psr_create() */

/*--------------------------------------------------------------------*/

void psr_delete (PATSET *psr, int delmap)
{                               /* --- delete pattern set red. object */
  free(psr->buf.items);         /* delete the item buffer, the border */
  free(psr->border);            /* and the identifier map */
  if (psr->map && delmap) idm_delete(psr->map);
  free(psr);                    /* delete the base structure */
}  /* psr_delete() */

/*--------------------------------------------------------------------*/

void psr_setbdr (PATSET *psr, ITEM size, RSUPP supp)
{                               /* --- set the decision border */
  assert(psr                    /* check the function arguments */
  &&    (size >= 0) && (supp >= 0));
  if (size <= psr->max) psr->border[size] = supp;
}  /* psr_setbdr() */

/*--------------------------------------------------------------------*/

void psr_addpat (PATSET *psr,
                 ITEM *items, ITEM size, RSUPP supp, void *orig)
{                               /* --- add a pattern to the set */
  FRQPAT *p;                    /* current pattern */

  assert(psr                    /* check the function arguments */
  &&    (psr->cur < psr->cnt)  && !psr->next
  &&    (items || (size <= 0)) && (supp >= 0) && orig);
  p = psr->pats +psr->cur;      /* get the current pattern */
  p->items  = items; p->size = size; p->supp = supp;
  p->orig   = orig;             /* store the pattern elements */
  int_qsort(p->items, (size_t)p->size, +1);
  psr->cur += 1;                /* switch to the next pattern */
}  /* psr_addpat() */

/*--------------------------------------------------------------------*/

void psr_addorig (PATSET *psr, void *orig)
{                               /* --- set the original pattern */
  FRQPAT *p;                    /* current pattern */

  assert(psr && psr->map        /* check the function arguments */
  &&    (psr->cur < psr->cnt) && psr->next && orig);
  p = psr->pats +psr->cur;      /* get the current pattern */
  p->items = psr->next;         /* set the current item location, */
  p->size  = 0;                 /* clear the pattern size, */
  p->orig  = orig;              /* and store the original pattern */
}  /* psr_addorig() */

/*--------------------------------------------------------------------*/

int psr_additem (PATSET *psr, const void *item)
{                               /* --- add an item to current pattern */
  FRQPAT *p;                    /* current pattern */
  IDENT   *id;                  /* to access the item identifiers */

  assert(psr && psr->map        /* check the function arguments */
  &&    (psr->cur < psr->cnt) && psr->next && (psr->rem > 0) && item);
  id = (ITEM*)idm_bykey(psr->map, item);
  if (!id) {                    /* get the item identifier */
    id = idm_add(psr->map, item, sizeof(void*), sizeof(IDENT));
    if (!id) return -1;         /* if the item does not exist, */
  }                             /* add it to the identifier map */
  *psr->next++ = (ITEM)*id;     /* store the item */
  p = psr->pats +psr->cur;      /* get the current pattern */
  p->size  += 1;                /* count the added item and */
  psr->rem -= 1;                /* reduce the remaining instances */
  return 0;                     /* return 'ok' */
}  /* psr_setbdr() */

/*--------------------------------------------------------------------*/

void psr_addsupp (PATSET *psr, RSUPP supp)
{                               /* --- finalize a pattern */
  FRQPAT *p;                    /* current pattern */

  assert(psr && psr->map        /* check the function arguments */
  &&    (psr->cur < psr->cnt) && psr->next && (supp >= 0));
  p = psr->pats +psr->cur;      /* get the current pattern */
  p->supp = supp;               /* store the pattern support */
  int_qsort(p->items, (size_t)p->size, +1);
  psr->cur += 1;                /* switch to the next pattern */
}  /* psr_addsupp() */

/*--------------------------------------------------------------------*/

size_t psr_reduce (PATSET *psr, int method, int addis)
{                               /* --- reduce the pattern set */
  size_t   i, k, n;             /* loop variables */
  FRQPAT   *p, *b;              /* to access the patterns */
  PATCMPFN *cmpfn;              /* pattern comparison function */
  int      r;                   /* result of pattern comparison */

  assert(psr                    /* check the function arguments */
  && (method >= PSR_NONE) && (method <= PSR_STRICT1));
  if (method <= PSR_NONE) return psr->cnt;
  cmpfn = psr_tab[method];      /* get the preference function */
  b     = &psr->buf;            /* and the pattern buffer */
  p     = psr->pats;            /* sort the patterns by size */
  obj_qsort(p, psr->cnt, sizeof(FRQPAT), +1, patcmp, NULL);
  for (i = 1; i < psr->cnt; i++) {
    for (k = 0; k < i; k++) {   /* traverse the pattern pairs */
      if (!p[i].orig && !p[k].orig)
        continue;               /* if both to be discarded, skip */
      if (isect(p+k, p+i, b) <= 0)
        continue;               /* if patterns are disjoint, skip */
      if (b->size < p[k].size){ /* if proper intersection exists */
        if (!addis              /* if to ignore intersections or */
        || (b->supp < psr->border[b->size]))
          continue;             /* if filtered out by border, skip */
        n = obj_bisect(b, p, psr->cnt, sizeof(FRQPAT), patcmp, NULL);
        if (patcmp(b, p+n, NULL) == 0)
          continue;             /* check whether intersection exists */
        for ( ; n < psr->cnt; n++) {
          if (subset(b, p+n) && (cmpfn(p+n, b, psr->border) < 0))
            p[n].orig = NULL;   /* filter superset-patterns */
        } }                     /* with pattern intersection */
      else {                    /* if pattern k is a subset of i */
        r = cmpfn(p+i, p+k, psr->border);
        if      (r > 0) p[k].orig = NULL;
        else if (r < 0) p[i].orig = NULL;
      }                         /* compare the two patterns and */
    }                           /* unmark the disfavored pattern */
    #ifdef PSR_ABORT
    if (sig_aborted()) break;   /* check for user abort */
    #endif
  }
  for (i = n = 0; i < psr->cnt; i++)
    if (p[i].orig) n++;         /* count the marked patterns */
  return n;                     /* return the number of patterns */
}  /* psr_reduce() */

