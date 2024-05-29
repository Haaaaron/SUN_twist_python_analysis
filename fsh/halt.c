#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
/* #include <float.h> */
#include <memory.h>

int
halt(char *s,void *p)
{
  fprintf(stderr,s,p);
  fprintf(stderr,"\n");
  exit(0);
  return(0);  /* just to get rid of an warning.. */
}
