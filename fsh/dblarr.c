#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
/* #include <float.h> */
#include <memory.h>
 
/* #include <io.h> */

#include "stuff.h"

double * dblarr(int size)
{
  double * p;

  p = (double *)calloc(size,sizeof(double));
  if (p == NULL) {
    fprintf(stderr," --- could not allocate double array of size %d\n",size);
    exit(0);
  }
  return (p);
}

float * fltarr(int size)
{
  float * p;

  p = (float *)calloc(size,sizeof(float));
  if (p == NULL) {
    fprintf(stderr," --- could not allocate float array of size %d\n",size);
    exit(0);
  }
  return (p);
}


int * intarr(int size)
{
  int * ip;
  
  ip = (int *)calloc(size,sizeof(int));
  if (ip == NULL) {
    fprintf(stderr," --- could not allocate int array of size %d\n",size);
    exit(0);
  }
  return(ip);
}
