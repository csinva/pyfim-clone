/*----------------------------------------------------------------------
  File    : patred.h
  Contents: pattern set reduction
  Author  : Christian Borgelt
  History : 2015.08.19 file created
----------------------------------------------------------------------*/
#ifndef __PATRED__
#define __PATRED__
#include "tract.h"
#include "report.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- pattern set reduction methods --- */
#define PSR_NONE        0       /* none (keep all patterns) */
#define PSR_COINS0      1       /* excess coincidences (zb,cb-ca) */
#define PSR_COINS1      2       /* excess coincidences (zb,cb-ca+1) */
#define PSR_ITEMS2      3       /* excess items        (za-zb+2,ca) */
#define PSR_COVER0      4       /* covered events       z   *c */
#define PSR_COVER1      5       /* covered events      (z-1)*c */
#define PSR_LENIENT0    6       /* combined lenient    (C+i+s) */
#define PSR_LENIENT1    7       /* combined lenient    (C+i+S) */
#define PSR_STRICT0     8       /* combined strict     (C+i+s) */
#define PSR_STRICT1     9       /* combined strict     (C+i+S) */

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- frequent pattern --- */
  ITEM    size;                 /* number of items (pattern size) */
  RSUPP   supp;                 /* support (number of occurrences) */
  ITEM    *items;               /* array of items (size elements) */
  void    *orig;                /* original item set object/marker */
} FRQPAT;                       /* (frequent pattern) */

typedef struct {                /* --- frequent pattern set --- */
  IDMAP   *map;                 /* item identifier map */
  RSUPP   *border;              /* decision border */
  ITEM    max;                  /* maximum pattern size */
  size_t  cnt;                  /* number of patterns */
  size_t  cur;                  /* index of current pattern */
  size_t  rem;                  /* remaining item instances */
  ITEM    *next;                /* next item position */
  FRQPAT  buf;                  /* buffer for a pattern */
  FRQPAT  pats[1];              /* pattern array */
} PATSET;                       /* (frequent pattern set) */

/*----------------------------------------------------------------------
  Pattern Set Reduction
----------------------------------------------------------------------*/
extern PATSET* psr_create  (size_t patcnt, ITEM patmax,
                            size_t extent, IDMAP *map);
extern void    psr_delete  (PATSET *psr, int delmap);
extern size_t  psr_patcnt  (PATSET *psr);
extern size_t  psr_patmax  (PATSET *psr);
extern void    psr_setbdr  (PATSET *psr, ITEM size, RSUPP supp);
extern RSUPP   psr_getbdr  (PATSET *psr, ITEM size);

extern void    psr_addpat  (PATSET *psr, ITEM *items, ITEM size,
                            RSUPP supp, void *orig);

extern void    psr_addorig (PATSET *psr, void *orig);
extern int     psr_additem (PATSET *psr, const void *item);
extern void    psr_addsupp (PATSET *psr, RSUPP supp);

extern size_t  psr_curcnt  (PATSET *psr);

extern size_t  psr_reduce  (PATSET *psr, int method, int addis);

extern void*   psr_getorig (PATSET *psr, size_t i);

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define psr_patcnt(s)     ((s)->cnt)
#define psr_patmax(s)     ((s)->max)
#define psr_getbdr(s)     ((s)->border)
#define psr_curcnt(s)     ((s)->cur)
#define psr_getorig(s,i)  ((s)->pats[i].orig)

#endif
