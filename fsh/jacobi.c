/************************************************
  
  computes the eigenvalues and eigenvectors of
  a real symmetric matrix a, size n^2.
  On return, d[n] contains the eigenvalues, and
  v[n*n] the eigenvectors.
  Returns the number of jacobi rotations

  + orders vectors and eigenvalues with descending order

************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double fabs(double);

#define verysmall(a,b) ((fabs(a)+b) - fabs(a) == 0)

void eigsrt(int n,double d[],double *v);

int jacobi(int n, double *ap,double d[],double *vp,int is_ordered)
{
  double *a,*b,*z,*v;
  int ip,iq,i,j,nrot;
  double sm,tresh,g,h,t,theta,tau,c,s;

  a = (double *)calloc(n*n,sizeof(double));
  v = (double *)calloc(n*n,sizeof(double));
  b = (double *)calloc(n,sizeof(double));
  z = (double *)calloc(n,sizeof(double));

  for (i=0; i<n; i++) for (j=0; j<n; j++) a[i*n+j] = ap[i*n+j];

  for (ip=0; ip<n; ip++) {
    for (iq=0; iq<n; iq++) v[ip*n+iq] = 0;
    v[ip*n+ip] = 1.0;
  }

  /* init b and d to the diagonal of a */
  
  for (ip=0; ip<n; ip++) {
    d[ip] = b[ip] = a[ip*n+ip];
    z[ip] = 0;
  }

  nrot = 0;
  for (i=1; i<=50; i++) {
    sm = 0.0;
    for (ip=0; ip<n-1; ip++) for (iq=ip+1; iq<n; iq++) sm += fabs(a[ip*n+iq]);
    if (sm == 0.0) {
      /* normal return, order the matrix */
      if (is_ordered) eigsrt(n,d,(double *)v);
      /* note here put the elements in natural order */
      for (i=0; i<n; i++) for (j=0; j<n; j++) vp[j*n+i] = v[i*n+j];
      free(z); free(b); free(v); free(a);
      return(nrot);
    }

    if (i < 4) tresh=0.2*sm/(n*n); else tresh = 0;
    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++) {
	g = 100.*fabs(a[ip*n+iq]);
	/* after 4 sweeps, skip the rotation if the off-diag. element
	   is small */
	if (i > 4 && verysmall(d[ip],g) && verysmall(d[iq],g)) a[ip*n+iq] = 0;
	else if (fabs(a[ip*n+iq]) > tresh) {
	  h = d[iq]-d[ip];
	  if (verysmall(h,g)) t = a[ip*n+iq]/h;
	  else {
	    theta = 0.5*h/a[ip*n+iq];
	    t = 1./(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0) t = -t;
	  }
	  c = 1./sqrt(1+t*t);
	  s = t*c;
	  tau = s/(1.0+c);
	  h = t*a[ip*n+iq];
	  z[ip] -= h; z[iq] += h;
	  d[ip] -= h; d[iq] += h;
	  a[ip*n+iq] = 0;
	  for (j=0; j<=ip-1; j++) {
	    g = a[j*n+ip]; h = a[j*n+iq];
	    a[j*n+ip] = g - s*(h + g*tau);
	    a[j*n+iq] = h + s*(g - h*tau);
	  }
	  for (j=ip+1; j<=iq-1; j++) {
	    g = a[ip*n+j]; h = a[j*n+iq];
	    a[ip*n+j] = g - s*(h + g*tau);
	    a[j*n+iq] = h + s*(g - h*tau);
	  }
	  for (j=iq+1; j<n; j++) {
	    g = a[ip*n+j]; h = a[iq*n+j];
	    a[ip*n+j] = g - s*(h + g*tau);
	    a[iq*n+j] = h + s*(g - h*tau);
	  }
	  for (j=0; j<n; j++) {
	    g = v[j*n+ip]; h = v[j*n+iq];
	    v[j*n+ip] = g - s*(h + g*tau);
	    v[j*n+iq] = h + s*(g - h*tau);
	  }
	  nrot++;
	} /* if fabs(a[ip*n+iq]) > tresh) */
      } /* iq =ip+1..n-1 */
    } /* ip = 0..n-2 */
    
    for (ip=0; ip<n; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0;
    }

  } /* i=1..50 */

  fprintf(stderr,"** Jacobi: 50 iterations - should not happen\n");
  free(z); free(b); free(v); free(a);
  return(nrot);
}



void eigsrt(int n,double d[],double *v)
{
  int i,k,j;
  double p;

  for (i=0; i<n-1; i++) {
    k = i;
    p = d[i];
    for (j=i+1; j<n; j++) if (d[j] >= p) p = d[k = j];
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j=0; j<n; j++) {
	p = v[j*n+i];
	v[j*n+i] = v[j*n+k];
	v[j*n+k] = p;
      }
    }
  }
}
