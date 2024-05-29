/************************************************************

  subroutine polyfit - does a polynomial fit to a vector

************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
/* #include <float.h> */
#include <memory.h>
 
/* #include <io.h> */
 
#include "stuff.h"  


#define tol 1e-17

double polyfit(int ndata,double *x,double *y,double *sig,
	       int deg,double *par,double *ep)
{
  int iserr,i,j,ma;
  double tmp,wmax,*u,*v,*w,*b;
  double thresh,chisq,sum;

  u = dblarr(ndata*(deg+1)); v = dblarr((deg+1)*(deg+1));
  w = dblarr(deg+1); b = dblarr(ndata);

  iserr = (sig != NULL);

  if (deg < 0) { fprintf(stderr,"polynomial degree < 0?\n"); exit(-1); }

  ma = deg+1;

  if (ma > ndata) {
    fprintf(stderr,"Underdetermined problem, filling with zeros\n");
    chisq = polyfit(ndata,x,y,sig,ndata-1,par,ep);
    for (i=ndata; i<ma; i++) par[i] = ep[i] = 0;
    return(chisq);
  }
    
  for (i=0; i<ndata; i++) {
    if (iserr) tmp = 1./sig[i]; else tmp = 1;
    b[i] = y[i]*tmp;
    for (j=0; j<ma; j++) {
      u[i*(deg+1)+j] = tmp;
      tmp *= x[i];
    }
  }

  svdcmp(u,ndata,ma,w,v);
  for (wmax=j=0; j<ma; j++) if (w[j] > wmax) wmax = w[j];
  thresh = tol*wmax;
  for (j=i=0; j<ma; j++) if (w[j] < thresh) { w[j] = 0; i++; }
  if (i) fprintf(stderr,"-- %d nearly singular values\n",i);
  svbksb(u,w,v,ndata,ma,b,par);
  chisq = 0;
  for (i=0; i<ndata; i++) {
    sum = 0;
    tmp = 1;
    for (j=0; j<ma; j++) {
      sum += par[j]*tmp;
      tmp *= x[i];
    }
    if (iserr) chisq += sqr((y[i]-sum)/sig[i]); else chisq += sqr(y[i]-sum);
  }

  for (j=0; j<ma; j++) {
    ep[j] = 0;
    for (i=0; i<ma; i++) if (w[i] > 0.) ep[j] += sqr(v[j*(deg+1)+i]/w[i]);
    ep[j] = sqrt(ep[j]);
  }

  free(b); free(u); free(w); free(v);

  return(chisq);
}
