/*******************************************************

stuff.h - contains prototypes for some math functions

*******************************************************/


#define sqr(a) ((a)*(a))
#define swap(a,b) { double swap_buf = a; a = b; b = swap_buf; }



#define greater(x,y) (((x) > (y)) ? (x) : (y))
#define bigger(x,y) (((x) > (y)) ? (x) : (y))
#define smaller(x,y) (((x) < (y)) ? (x) : (y))

#define getoptnum(par1,par2,val){		\
  char *sp;					\
  sp = ss;					\
  if (!*ss && (argc-1)) ss = (argv+1)[0];	\
  if (sscanf(ss,par1,par2) != 1) {		\
    *par2 = val; ss = sp;			\
  } else {					\
    if (!*sp && --argc) argv++;			\
    ss = strchr(ss,0);				\
  }						\
}

#define getnum(par1,par2){                \
  if (!*ss && --argc) ss = (++argv)[0];   \
  if (sscanf(ss,par1,par2) != 1) {        \
     fprintf(stderr,"%s",usage);               \
     exit(-1);                            \
  }                                       \
  ss = strchr(ss,0);}

#define get2num(str,p1,p2){               \
  if (!(*ss) && --argc) ss = (++argv)[0]; \
  if (sscanf(ss,str,p1,p2) != 2) {        \
    fprintf(stderr,"%s",usage);                \
    exit(-1);                             \
  }                                       \
  ss = strchr(ss,0);}

#define get3num(str,p1,p2,p3){            \
  if (!(*ss) && --argc) ss = (++argv)[0]; \
  if (sscanf(ss,str,p1,p2,p3) != 3) {     \
    fprintf(stderr,"%s",usage);                \
    exit(-1);                             \
  }                                       \
  ss = strchr(ss,0);}

#define get4num(str,p1,p2,p3,p4){		\
  if (!(*ss) && --argc) ss = (++argv)[0];	\
  if (sscanf(ss,str,p1,p2,p3,p4) != 4) {	\
    fprintf(stderr,"%s",usage);			\
    exit(-1);					\
  }						\
  ss = strchr(ss,0);}

#define get5num(str,p1,p2,p3,p4,p5){		\
  if (!(*ss) && --argc) ss = (++argv)[0];	\
  if (sscanf(ss,str,p1,p2,p3,p4,p5) != 5) {	\
    fprintf(stderr,"%s",usage);			\
    exit(-1);					\
  }						\
  ss = strchr(ss,0);}


#define getlist(v,i){					\
  if (!(*ss) && --argc) ss = (++argv)[0];		\
  i = 0;						\
  while (*ss) {						\
    v[i] = ss; 						\
    i++;						\
    if (strchr(ss,';') == NULL) break;			\
    ss = strchr(ss,';');				\
    *ss = 0;						\
    ss++;						\
  }							\
  if (i <= 0) {						\
    fprintf(stderr,"%s",usage);				\
    exit(-1);						\
  }							\
  ss = strchr(ss,0);}

#define getstring(v){					\
  if (!(*ss) && --argc) ss = (++argv)[0];		\
  if (!(*ss)) {						\
    fprintf(stderr,"%s",usage);				\
    exit(-1);						\
  }							\
  v = ss;						\
  ss = strchr(ss,0);}


#define getnumlist(v,i,format){			\
  if (!(*ss) && --argc) ss = (++argv)[0];	\
  i = 0;					\
  while (*ss) {					\
    if (sscanf(ss,format,&v[i++]) != 1) {	\
      fprintf(stderr,"%s",usage); exit(-1);		\
    }						\
    if (strchr(ss,',') == NULL) break;		\
    ss = strchr(ss,',');			\
    ss++;					\
  }						\
  if (i <= 0) {					\
    fprintf(stderr,"%s",usage);			\
    exit(-1);					\
  }						\
  ss = strchr(ss,0);}


#define get1or3num(str1,str3,p1,p2,p3,n){	\
  if (!(*ss) && --argc) ss = (++argv)[0];	\
  if (sscanf(ss,str3,&p1,&p2,&p3) != 3) {	\
    if (sscanf(ss,str1,&p1) != 1) {		\
      fprintf(stderr,"%s",usage);			\
      exit(-1);					\
    }						\
    n = 1; 					\
    p3 = p1; p2 = 1;				\
  } else n = 3;					\
  ss = strchr(ss,0);}



/********* headers for unformatted io **********/

typedef struct {
  int headerid,f1,headersize,f2;
  int n_double,f3,n_long,f4,n_int,f5,n_char,f6;
  int lx,f7,ly,f8,lz,f9,lt,f10;
  int d1,f11,d2,f12,d3,f13,d4,f14,d5,f15,d6,f16,d7,f17,d8,f18;
} ll_header;

typedef struct {
  long headerid,headersize;
  long n_double,n_long,n_float,n_char;
  long lx,ly,lz,lt;
  long d1,d2,d3,d4,d5,d6,d7,d8;
} e_header;

typedef struct {
  int headerid,headersize;
  int n_double,n_long,n_float,n_char;
  int lx,ly,lz,lt;
  int d1,d2,d3,d4,d5,d6,d7,d8;
} i_header;

#define E_HEADER_ID 91919191

/* and a couple of protos */
int readheader(FILE *ff,e_header *h);
int skipheader(FILE *ff);
long readdata(FILE *ff,double *tmparr);

/* read_ascii */
int read_ascii_info( FILE *f,       // input file
		     char *label,   // label string - if NULL, plain line-ascii
		     int *rows,      // number of rows
                     char **legend,  // legend string
                     int iskip
		     );
int read_ascii( FILE *f,        // input file
		char *label,    // label string - if NULL, plain line-ascii
		double *dat,     // data array
                int iskip
		);

/* other prototypes */

double calclist(double d[],int dl,char *cmd);
double * dblarr(int size);float * fltarr(int size);
int * intarr(int size);
double confidence(double chisq,int dof);
int gaussj(double* a,int n, int np, double* b,int m);
void fitfun(double *x,double *y,double *sig,int ndata,double *a,int ma,
	    double *covar,int *lista,int mfit,double *chisq,int print,
	    double funcs());
void covarfit(double *x,double *y,double *cmat,int ndata,double *a,int ma,
	      double *covar,int *lista,int mfit,double *chisq,int print,
	      double funcs());
void jackfit(double *x,double *y,double * sig,int ndata,int n1, int n2,
	     int jack,double *a,int ma,
	     double *covar,int *lista,int mfit,double *chisq,int print,
	     int fullcov,int simplex,double funcs());
void jack_fit(double *x,double *y,double * sig,int ndata,int n1, int n2,
	      int jack,double *a,int ma,
	      double *covar,int *lista,int mfit,double *chisq,int print,
	      int fullcov,int simplex,double funcs(),double *av);
double fitfun_s(double *x,double *y,double *sig,int ndata,double *a,int ma,
		double funcs());
double brent(double ax,double bx,double cx,double f(),double tol,double *xmin);
double polyfit(int ndata,double x[],double y[],double sig[],
	       int deg,double par[],double ep[]);
double nelder(int ndim,double p[],double ftol,double funk(),int *i);
void simplexfit(double *x,double *y,double *sig,int ndata,double *a,int ma,
		int *lista,int mfit,double *chisq,int print,double funcs());
int jacobi(int n, double *ap,double d[],double *v,int is_ordered);
double simpson(double x[], double y[], double res[], int r);

void svdecomp(double *a,int n,double *b, int m);
void svdcmp(double *a, int m, int n, double *w, double *v);
void svbksb(double *u, double *w, double *v, int m, int n, double *b, double *x);


/* these defined in diag_corr.c */
void diag_corr( double ***matrix, // input matrix
		int msize,        // size of matrix
		int d1,           // min. distance
		int d2,           // max. distance
		int block,        // number of jackknife blocks
		int dist_invert,  // distance at which we form inversion
		                  // < 0 : not do it at all
		int usesqrt,      // Use symmetric/left/right sqrt -operation
		int diag_dist,    // distance at which we take eigenvectors
	                          // if 0, eigenvectors at all distances
		double **resv,    // output eigenvalues
		double **err,     // output eigenvalue errors
		double ***evalue, // output eigenvalue blocks, if desired
		double **rvec,    // output vectors, if desired
		double **evec     // and errors of those
		);

void localmass( double ***evalue, int msize,  // vector value and size
		int block,                    // jack block 
		int d1, int d2,               // min max distance to do avg
		double *amass, double *err,   // avg. mass + error
		double *chisq,                // chisqr
		double **jmass,               // jack block mass
		double **dmass, double **derr,// dist. mass + error
		int latlen);                  // periodicity distance

double ***alloc3darray(int d2,int s2, int s3);
void free3darray(double ***r,int d2,int s2);
double **alloc2darray(int d1,int s );
void free2darray(double **r,int d1);

