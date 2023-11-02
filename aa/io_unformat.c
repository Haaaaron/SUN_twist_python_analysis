/*
 *  UNFORMATTED IO SYSTEM
 *  Kari Rummukainen 1990 - 1999
 */
 

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
/* #include <float.h> */
#include <memory.h>

#include "stuff.h"


static int inv_bytes,long_mode,block,dblock,lblock,fblock,cblock,iblock;
static double *dtmparr;
static float *ftmparr;
static long *ltmparr;
static int *itmparr;
static char *ctmparr;
static int msg = 0;

#define l_h (sizeof(e_header)/sizeof(long))

typedef union {
  e_header h;
  long l[l_h];
} h_union;

typedef union {
  i_header h;
  int l[l_h];
} i_union;

#define ll_h (sizeof(ll_header)/sizeof(int))

typedef union {
  ll_header h;
  int l[ll_h];
} ll_union;



/**************************************************
 * invert the byte ordering
 */

long 
swaplong(long a)
{
  union { 
    long l;
    char c[sizeof(long)];
  } t1,t2;
  int i;

  t1.l = a;
  for (i=0; i<sizeof(long); i++) t2.c[i] = t1.c[sizeof(long)-1-i];
  return(t2.l);
}


int
swapint(int a)
{
  union { 
    int l;
    char c[sizeof(int)];
  } t1,t2;
  int i;

  t1.l = a;
  for (i=0; i<sizeof(int); i++) t2.c[i] = t1.c[sizeof(int)-1-i];
  return(t2.l);
}


double 
swapdouble(double a)
{
  union { 
    double d;
    char c[sizeof(double)];
  } t1,t2;
  int i;

  t1.d = a;
  for (i=0; i<sizeof(double); i++) t2.c[i] = t1.c[sizeof(double)-1-i];
  return(t2.d);
}


float
swapfloat(float a)
{
  union { 
    float d;
    char c[sizeof(double)];
  } t1,t2;
  int i;

  t1.d = a;
  for (i=0; i<sizeof(float); i++) t2.c[i] = t1.c[sizeof(float)-1-i];
  return(t2.d);
}


/************************************************************
 *  Main I/O stuff here
 */

int
readheader(FILE *ff,e_header *h)
{
  i_union i;
  h_union l;
  ll_union ll;
  int j,k;

  /* first, debug the type of the header by reading 1st 4 ints */

  long_mode = 1;
  if (sizeof(long) == 2*sizeof(int)) { 
    long b[2];
    /* now check if we can do longs here */
    fread(&b,sizeof(long),2,ff); rewind(ff);

    if (b[0] == E_HEADER_ID && b[1] == sizeof(e_header)) {
      long_mode = 0; inv_bytes = 0;
    } else if (swaplong(b[0]) == E_HEADER_ID && 
	       swaplong(b[1]) == sizeof(e_header)) {
      long_mode = 0; inv_bytes = 1;
    }
  }
  
  if (long_mode == 1) {
    int err = 0;
    int a[4];
    
    fread(&a,sizeof(int),4,ff); rewind(ff);

    if (a[0] == E_HEADER_ID) {
      inv_bytes = 0;   /* normal byte ordering */
      if (a[1] == sizeof(i_header)) long_mode = 1;
      else if (a[2] == sizeof(ll_header)) long_mode = 2;
      else { fprintf(stderr,"header size error 1\n"); err = 1; }
    } else if (swapint(a[0]) == E_HEADER_ID) {
      inv_bytes = 1;
      if (swapint(a[1]) == sizeof(i_header)) long_mode = 1;
      else if (swapint(a[2]) == sizeof(ll_header)) long_mode = 2;
      else { fprintf(stderr,"header size error 2\n"); err = 1; }
    } else if (a[1] == E_HEADER_ID) {
      inv_bytes = 0;   /* normal byte ordering */
      if (a[3] == sizeof(ll_header)) long_mode = 3;
      else { fprintf(stderr,"header size error 3\n"); err = 1; }
    } else if (swapint(a[1]) == E_HEADER_ID) {
      inv_bytes = 1;
      if (swapint(a[3]) == sizeof(ll_header)) long_mode = 3;
      else { fprintf(stderr,"header size error 3\n"); err = 1; }
    } else {
      // Assume now plain ascii file!
      return(-1);
    }   
    
    if (err) { 
      if (err) fprintf(stderr,"Header error\n");
      fprintf(stderr,"4 first:       %d  %d  %d  %d\n",a[0],a[1],a[2],a[3]);
      fprintf(stderr,"swapped bytes  %d  %d  %d  %d\n",
	      swapint(a[0]),swapint(a[1]),swapint(a[2]),swapint(a[3]));
      if (err) exit(0);
    }
  }

  if (!msg) {
    if (inv_bytes)      fprintf(stderr,"* inverted byte ordering ");
    if (long_mode == 0) fprintf(stderr,"* (long) ");
    if (long_mode == 1) fprintf(stderr,"* (int) ");
    if (long_mode == 2) fprintf(stderr,"* (int)(fill) ");
    if (long_mode == 3) fprintf(stderr,"* (fill)(int) ");
    fprintf(stderr,"\n");
    msg = 1;
  }

  if (long_mode == 0) {           /* normal longs */
    fread(&l.h,sizeof(e_header),1,ff);
    if (inv_bytes) for (j=0; j<l_h; j++) l.l[j] = swaplong(l.l[j]);
  } else if (long_mode == 1) {    /* ints */
    fread(&i.h,sizeof(i_header),1,ff);

    if (!inv_bytes) 
      for (j=0; j<l_h; j++) l.l[j] = i.l[j];
    else 
      for (j=0; j<l_h; j++) l.l[j] = swapint(i.l[j]);
  } else if (long_mode >= 2) {    /* now long is longer than our long */
    fread(&ll.h,sizeof(ll_header),1,ff);
    
    for (k=0, j=( (long_mode == 2) ? 0 : 1 ); j<ll_h; j+=2, k++) {
      l.l[k] = ll.l[j];
      if (inv_bytes) l.l[k] = swaplong(l.l[k]);
    }
  }

  *h = l.h;
  if (h->lz < 1) h->lz = 1; 
  if (h->lt < 1) h->lt = 1;

  lblock = iblock = 0;
  block  = dblock = l.h.n_double;
  if (dblock) dtmparr = dblarr(l.h.n_double);
  block += fblock = l.h.n_float;
  if (fblock) ftmparr = (float *)calloc(l.h.n_float,sizeof(float));
  block += cblock = l.h.n_char;
  if (cblock) ctmparr = (char *)calloc(l.h.n_char,sizeof(char));

  if (l.h.n_long) {
    if (long_mode == 0) {
      block += lblock = l.h.n_long;
      ltmparr = (long *)calloc(l.h.n_long,sizeof(long));
    } else if (long_mode == 1) {
      block += iblock = l.h.n_long;
      itmparr = (int *)calloc(l.h.n_long,sizeof(int));
    } else if (long_mode >= 2) {
      block += iblock = l.h.n_long;
      itmparr = (int *)calloc(2*l.h.n_long,sizeof(int));
    }
  }

  return (block);

}

int
skipheader(FILE *ff)
{
  e_header e;
  i_header i;
  
  if (long_mode == 0) fread(&e,sizeof(e_header),1,ff);
  else if (long_mode == 1) fread(&i,sizeof(i_header),1,ff);
  else if (long_mode >= 2) {
    fread(&i,sizeof(i_header),1,ff);
    fread(&i,sizeof(i_header),1,ff);
  }
  return(1);
}


#define itmp_index(i) ((long_mode < 2) ? i : 2*i + long_mode-2)

long
readdata(FILE *ff,double *arr)
{
  int ik,k,i;
  long lk;
  int end;
  
  if (dblock) end = (fread(dtmparr,sizeof(double),dblock,ff) != dblock);
  if (fblock) end = (fread(ftmparr,sizeof(float),fblock,ff) != fblock);
  if (lblock) end = (fread(ltmparr,sizeof(long),lblock,ff) != lblock);
  if (iblock && long_mode < 2)  
    end = (fread(itmparr,sizeof(int),iblock,ff) != iblock);
  if (iblock && long_mode >= 2) 
    end = (fread(itmparr,sizeof(int),2*iblock,ff) != 2*iblock);
  if (cblock) 
    end = (fread(ctmparr,sizeof(char),cblock,ff) != cblock);

  if (!end) {
    k = 0;
    if (long_mode == 1) {
      end = (fread(&ik,sizeof(int),1,ff) != 1); 
      if (inv_bytes) lk = swapint(ik); else lk = ik;
    } else if (long_mode == 0) {
      end = (fread(&lk,sizeof(long),1,ff) != 1); 
      if (inv_bytes) lk = swaplong(lk);
    } else {
      int a[2];
      end = (fread(&a,sizeof(int),2,ff) != 2);
      if (inv_bytes) lk = swapint(a[long_mode-2]);
      else lk = a[long_mode-2];
    }
  }

  if (end) return(0);

  if (inv_bytes) {
    for (i=0; i<dblock; i++) arr[k++] = swapdouble(dtmparr[i]);
    for (i=0; i<fblock; i++) arr[k++] = swapfloat(ftmparr[i]);
    for (i=0; i<lblock; i++) arr[k++] = swaplong(ltmparr[i]);
    for (i=0; i<iblock; i++) arr[k++] = swapint(itmparr[itmp_index(i)]);
    for (i=0; i<cblock; i++) arr[k++] = ctmparr[i];
  } else {
    for (i=0; i<dblock; i++) arr[k++] = dtmparr[i];
    for (i=0; i<fblock; i++) arr[k++] = ftmparr[i];
    for (i=0; i<lblock; i++) arr[k++] = ltmparr[i];
    for (i=0; i<iblock; i++) arr[k++] = itmparr[itmp_index(i)];
    for (i=0; i<cblock; i++) arr[k++] = ctmparr[i];
  }
    
  return(lk);
}


