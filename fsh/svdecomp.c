/***********************************************************

int svdecomp(a,n,b,m)
solves a*x=b;
a: input n*n matrix, in an array of n*n
b: rhs, n*m vector

returns solution in b.

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

#define zerolimit 1e-15

void svdecomp(a,n,b,m)
int n,m;
double *a,*b;
{
  int i,j,k;
  double *u, *w, *v, *x;
  double wmin,wmax,condnum;

  u = dblarr(n*n); w = dblarr(n); v = dblarr(n*n); x = dblarr(n);

  /* not transpose the guy.. */
  for (j=0; j<n; j++) for (i=0; i<n; i++) u[i*n+j] = a[i*n+j];

  /* decompose with numerical receipes..*/
  svdcmp(u,n,n,w,v);

  if (n == 0) {
    fprintf(stderr,"* svd: no convergence\n");
    exit(-1);
  }

  wmax = 0;
  j = 0;
  condnum = 1;
  for (i=0; i<n; i++) if (w[i] > wmax) wmax = w[i];
  wmin = zerolimit*wmax;
  for (i=0; i<n; i++) {
    condnum = (condnum < w[i]/wmax) ? condnum : w[i]/wmax;
    if (w[i] < wmin) {
      w[i] = 0.;
      j++;
    }
  }

  if (j) fprintf(stderr,"* svd: Zeroed %d elements\n",j);
  if (condnum < 10*zerolimit)
    fprintf(stderr,"* svd: condition %g\n",1./condnum);

  /* singular value back substitution */

  for (i=0; i<m; i++) {
    svbksb(u,w,v,n,n,b+(n*i),x);
    for (j=0; j<n; j++) b[n*i+j] = x[j];
  }

  /* invert the matrix to a */

  for (i=0; i<n; i++) for (j=0; j<n; j++) {
    a[i*n+j] = 0;
    for (k=0; k<n; k++) if (w[k] != 0.0) a[i*n+j] += v[i*n+k]*u[j*n+k]/w[k];
  }
  
  free(x); free(v); free(w); free(u);

}



#define sign(a,b) ((b)*(fabs((a)/(b))))

void svdcmp(a,m,n,w,v)
int m,n;
double *a,*w,*v;
{
  double *rv1;
  double c,g,s,scale,f,h,x,y,z,anorm;
  int i,l,j,jj,k,nm,its,breakval;

  rv1 = dblarr(n);

  if (m < n) {
    fprintf(stderr,"Underdetermined system - fill matrix with zeros\n");
    exit(-1);
  }

  g=scale=anorm=0.0;
  for (i=0; i<n; i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m-1) {
      for (k=i; k<m; k++) scale=scale+fabs(a[k*n+i]);
      if (fabs(scale) > 1e-20) {
	for (k=i; k<m; k++) {
	  a[k*n+i]=a[k*n+i]/scale;
	  s=s+a[k*n+i]*a[k*n+i];
	}
	f=a[i*n+i];
	g=-sign(sqrt(s),f);
	h=f*g-s;
	a[i*n+i]=f-g;	
	if (i != n-1) {
	  for (j=l; j<n; j++) {
	    s=0.0;
	    for (k=i; k<m; k++) s=s+a[k*n+i]*a[k*n+j];
	    f=s/h;
	    for (k=i; k<m; k++) a[k*n+j]=a[k*n+j]+f*a[k*n+i];
	  }
	}
	for (k=i; k<m; k++) a[k*n+i]=scale*a[k*n+i];
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if ((i <= m-1) && (i != n-1)) {
      for (k=l; k<n; k++) scale=scale+fabs(a[i*n+k]);
      if (scale > 1e-15) {
	for (k=l; k<n; k++) {
	  a[i*n+k]=a[i*n+k]/scale;
	  s=s+a[i*n+k]*a[i*n+k];
	}
	f=a[i*n+l];
	g=-sign(sqrt(s),f);
	h=f*g-s;
	a[i*n+l]=f-g;
	for (k=l; k<n; k++) rv1[k]=a[i*n+k]/h;
	if (i != m-1) {
	  for (j=l; j<m; j++) {
	    s=0.0;
	    for (k=l; k<n; k++) s=s+a[j*n+k]*a[i*n+k];
	    for (k=l; k<n; k++) a[j*n+k]=a[j*n+k]+s*rv1[k];
	  }
	}
	for (k=l; k<n; k++) a[i*n+k]=scale*a[i*n+k];
      }
    }

    anorm=greater(anorm,(fabs(w[i])+fabs(rv1[i])));
  }

  for (i=n-1; i>=0; i--) {
    if (i < n-1) {
      if (fabs(g) > 1e-15) {
	for (j=l; j<n; j++) v[j*n+i]=(a[i*n+j]/a[i*n+l])/g;
	for (j=l; j<n; j++) {
	  s=0.0;
	  for (k=l; k<n; k++) s=s+a[i*n+k]*v[k*n+j];
	  for (k=l; k<n; k++) v[k*n+j]=v[k*n+j]+s*v[k*n+i];
	}
      }
      for (j=l; j<n; j++) v[i*n+j]=v[j*n+i]=0.;
    }
    v[i*n+i]=1.0;
    g=rv1[i];
    l=i;
  }

  for (i=n-1; i>=0; i--) {
    l=i+1;
    g=w[i];
    if (i < n-1) for (j=l; j<n; j++) a[i*n+j]=0.0;
    if (fabs(g) > 1e-15) {
      g=1.0/g;
      if (i != n-1) {
	for (j=l; j<n; j++) {
	  s=0.0;
	  for (k=l; k<m; k++) s=s+a[k*n+i]*a[k*n+j];
	  f=(s/a[i*n+i])*g;
	  for (k=i; k<m; k++) a[k*n+j]=a[k*n+j]+f*a[k*n+i];
	}
      }
      for (j=i; j<m; j++) a[j*n+i]=a[j*n+i]*g;
    } else {
      for (j=i; j<m; j++) a[j*n+i]=0.0;
    }
    a[i*n+i]=a[i*n+i]+1.0;
  }
  for (k=n-1; k>=0; k--) {
    for (its=0; its<30; its++) {
      breakval = 0;
      for (l=k; l>=0; l--) {
	nm = l-1;
/*	
	  if (fabs(rv1[l]/anorm) < 1e-13) goto label2;
	  if (fabs(w[nm]/anorm)  < 1e-13) goto label1;
*/	  
	if (fabs(rv1[l])+anorm == anorm) goto label2;
	if (fabs(w[nm])+anorm == anorm)  goto label1;

      }
      fprintf(stderr,"Svd: should not be here!! anorm %.10lg \n",anorm);
      exit(0);
      

    label1:
      if (breakval != 2) {
	c=0.0;
	s=1.0;
	for (i=l; i<=k; i++) {
	  f=s*rv1[i];
	  if (fabs(f)+anorm != anorm) {
	    g=w[i];
	    h=sqrt(f*f+g*g);
	    w[i]=h;
	    h=1.0/h;
	    c= (g*h);
	    s=-(f*h);
	    for (j=0; j<m; j++) {
	      y=a[j*n+nm];
	      z=a[j*n+i];
	      a[j*n+nm]=(y*c)+(z*s);
	      a[j*n+i]=-(y*s)+(z*c);
	    }
	  }
	}
      }
      /* entry for breakval 2 */
    label2:      
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k]=-z;
	  for (j=0; j<n; j++) v[j*n+k]=-v[j*n+k];
	}
	goto label3;
      } 
      
      if (its == 200-1) {
	fprintf(stderr,"svd: no solution in 200 iterations\n");
	exit(-1);
      }
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=sqrt(f*f+1.0);
      f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
      c=1.0;
      s=1.0;
      for (j=l; j<=nm; j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=sqrt(f*f+h*h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f= (x*c)+(g*s);
	g=-(x*s)+(g*c);
	h=y*s;
	y=y*c;
	for (jj=0; jj<n; jj++) {
	  x=v[jj*n+j];
	  z=v[jj*n+i];
	  v[jj*n+j]= (x*c)+(z*s);
	  v[jj*n+i]=-(x*s)+(z*c);
	}
	z=sqrt(f*f+h*h);
	w[j]=z;
	if (fabs(z) > 1e-15) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f= (c*g)+(s*y);
	x=-(s*g)+(c*y);
	for (jj=0; jj<m; jj++) {
	  y=a[jj*n+j];
	  z=a[jj*n+i];
	  a[jj*n+j]= (y*c)+(z*s);
	  a[jj*n+i]=-(y*s)+(z*c);
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  label3: continue;
    /* entry here from the break 3 */
  }

  free(rv1);
}




void svbksb(u,w,v,m,n,b,x)
int m,n;
double *u,*w,*v,*b,*x;
{
  double *tmp,s;
  int i,j,jj;

  tmp = dblarr(n);
  for (j=0; j<n; j++) {
    s=0.;
    if (w[j] != 0.) {
      for (i=0; i<m; i++) s=s+u[i*n+j]*b[i];
      s=s/w[j];
    }
    tmp[j]=s;
  }
  for (j=0; j<n; j++) {
    s=0.;
    for (jj=0; jj<n; jj++) s=s+v[j*n+jj]*tmp[jj];
    x[j]=s;
  }
  free(tmp);
}

