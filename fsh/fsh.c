/*************************************************************************
 *
 *  fsh.c
 *
 *  Kari Rummukainen 1991-1999
 *
 *  This c-program calculates histograms written with special
 *  c-unformatted output
 *
 *************************************************************************/

/* TO COMPILE:
 *
 * cc -O -o fsh fsh.c io_unformat.c calclist.c dblarr.c svdecomp.c polyfit.c jacobi.c halt.c -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <memory.h>

#include "stuff.h"

/* #include <io.h> */

#define BETAMAX 2000
#define MAXRUN  300  /* Max number of runs */
#define MAX_DAT 50   /* Max number of data analyzed at one go */

#define max(x,y) (((x) > (y)) ? (x) : (y))

typedef struct {
  double r,i;
} complex;

#define pi  3.141592653589793
#define pi2 (2.0*pi)

#define scan_whitespace(r)     { while (*r == ' ' || *r == '\t') r++; }
#define scan_not_whitespace(r) { while (*r && *r != ' ' && *r != '\t' && *r != '\n') r++; }


/* prototypes... */

int FSinit(int num,int *nvec,double *beta,double *beta2,
	   double **act,double **act2,
	   double epsilon,int tc,int maxl,int init,char *ss,
	   int jack,int irange,double *w[]);

void FSCalc(double **dv[],int nd,int isfree,double mom,int binder,int cross);
void FSCpart(double *imbeta,int ni);
void FShgvec(double **dv[],
	     int bins,double cut_value,int eq_weight,int eq_height);
void FShgram(double **d[],int nd,double beta, double beta2, 
	     int bins, double cut, int eqw, int eqh);
void FShgram2(double **d1,double **d2,double beta, double beta2,
	      int bins1, int bins2);
void FSeig(double **dv[],int nd,int eprint,char **lists);
void FSm2mat(double **dv[],int nd,double *val);
void FSJack(double **dv[],int nd,int jack,int isfree,int jprint,char *js,
	    char *ss,double mom,int binder,int cross);
void FShgvecJ(double **dv[],int jack,int jprint,char *js,
	      char *ss,
	      int bins,double cut_value,int eq_weight,int eq_height);
void FShgJack(double **dat[],int nd,double beta, double beta2, int bins, double cut,
	      int jack,int jprint,char *js,char *ss, int eqw, int eqh);
void FSiniJack(int iter, int jack, char *pre);
void FSeigJack(double **dv[],int nd,int jack,int jprint,char *js,char *ss);
void errorcalc(int nconf,double svect[],double *avep,double *sigp,
	       double *fsp, double *fser, double *w);

int inifil(char *fn,double *inv_vol);
void getdbldat(double d[]);
double getdouble();
int closefil();

void map_beta(double *b1,double *b2);
double beta2_Higgs(double mh,double betah);
double beta_u1su2(double p_y,double p_x);
double beta_su3h(double p_y,double *beta4);

int errors(int n,int nv[], double beta[], double *a[],int tc,int jack,double *w[]);

void getextrema(double hist[],int bins,double min1,double max1, double cut,
		double *maxv1,double *minv,double *maxv2,
		double *maxl1,double *minl,double *maxl2,
		double *cu1,double *cu2);

char expnum[] = " expecting a number after the option\n";
char usage[] = " usage: fsp [-opt] datafile > results\n\
  v              : verbose\n\
  X              : raw data without FS analysis\n\
  b value OR min:step:max : reweight to the given value OR vector of values\n\
  B val OR min:step:max : reweight also wrt. second reweight component \n\
  s str          : use str for aux. filenames [datafile]\n\
  S(+)           : load settings from .values (append to it) [save to it]\n\
  x              : fit data 2 at a time\n\
  f/q            : F/S\n\
  z b_im/min:step:max : partition function, with imaginary beta\n\
  d num          : data number\n\
  D '#1+2*#2;#3' : get data and do arithmetics [+-*/^() ,]\n\
  2              : second moment of data\n\
  M power        : moment power of data\n\
  Z              : calculate cross-correlation <(1-<1>)(2-<2>)> from D\n\
  4              : calculate Binder's 4th cumulant\n\
  a              : for moments: use jackknife block averages (normally global)\n\
  V              : divide value by volume\n\
  E              : 2nd moment eigenmatrix analysis\n\
    EC           : print in calc-format\n\
  p/P value      : low/highpass filter in data\n\
  C cut_value    : discontinuity value in data\n\
  G cut1:cut2    : a gap in the data to define cut\n\
  h              : fs-histogram\n\
    A bins[,b2]  : number of bins\n\
    W            : equal weight beta search, from cut [W+ search beta2 too]\n\
    H            : equal height beta search, from cut [H+] \n\
    L l1,l2,l3   : number of points to find the extrema [5,10,5]\n\
    F            : assume constant (flat) central part\n\
    R            : print the extrema data\n\
    o/O value    : low/highpass x-value in histograms\n\
  g              : double histogram with datas d1,d2\n\
  Q              : print delta F, using cut_value(s) to separate phases\n\
  c              : do not calculate time correlations\n\
  j num          : jackknife with num blocks\n\
  J name/=       : print jackknifed data/to stdout\n\
  l number       : maximum loop limit, 0-no limit [100]\n\
  r double       : relative accuracy used [1e-8]\n\
  t number       : thermalize n iterations\n\
  I n1:n2        : use only lines n1 to n2\n";


#define NO_CUT_VALUE -1313
static int noauto = 0;
static int verbose = 0;
static int is_cut = 0, delta = 0, deltaF = 0;
static int jstdout = 0, hp_low = 0, hp_high = 0;
static double hprint_low = NO_CUT_VALUE, hprint_high = NO_CUT_VALUE;
static double cut_value = NO_CUT_VALUE;
static double cut_v1,cut_v2;
static double inv_vol = 1;
static int is_ascii = 0;
static int max_r1 = 5, max_r2 = 5, min_r = 10, print_ext = 0;
static e_header h;
static int flat = 0;
static int is_susy = 0, is_susy23 = 0, is_u1su2 = 0, Higgs = 0, is_su3h = 0;
static int isact2 = 0,b2search = 0;
static double mass_U = 0.0,betag = 0.0, susy23_amul, susy23_mh, p_x = 0.0, p_x_global = 0.0;
static double betav[BETAMAX],beta2v[BETAMAX],imbetav[BETAMAX];
static int betanumber = 0;
static int is_global_ave = 1;
static int read_2_beta = 0;

#define p_beta(x) (is_susy ? 1.0/sqrt(x) : x)

#define RW_STRLEN 600

int
main(int argc,char *argv[])
{
  int nmeas,thermo;
  double dum,tmp[500],bmin,bstep,bmax;
  FILE *dataf;
  double *act[MAXRUN], *w[MAXRUN], *act2[MAXRUN], *ap, *a2p, *dp;
  double **dv[MAX_DAT];
  double epsilon;
  double beta[MAXRUN],beta2[MAXRUN];
  int nv[MAXRUN];
  char  name[500],npus[500],hname[500],js[200];
  char *ss, *chp,*lists[MAXRUN];
  int i,j,n,raw,nd,datdim;
  int ntop,n1,n2,reweight_data,reweight_cmd,act2_data;
  char reweight_string[RW_STRLEN],reweight_string2[RW_STRLEN];
  double sh,mh,reweight_mult,mom, *pt;
  int tc,maxl,init,irange;
  int jack,jprint;
  int hg2,bins,bins2,dnum,eq_weight,eq_height;
  int sok,fok,binder,data,voldiv,eig,hgvec,ishg,ispart,cross;

  voldiv = flat = reweight_data = hgvec = ishg = reweight_cmd = 0;
  sok = fok = mom = binder = eq_weight = eq_height
    = irange = data = eig = hg2 = ispart = 0;
  dnum = nd = cross = 0;
  
  n1 = 0;
  n2 = 100000;
  betanumber = 0;
  epsilon = 1e-8;
  thermo = 0;
  maxl = 100;
  name[0] = npus[0] = 0;
  init = 0;
  tc = 1;
  raw = jack = jprint = 0;
  bins = bins2 = 0;

  if (argc <= 1) halt(usage,NULL);

  while (--argc > 0 && (*++argv)[0] == '-') {
    ss = argv[0] + 1;

    while (*ss) {
      switch(*ss++) {

      case 'v': verbose = 1; break;
	
      case 'x': irange  = 1; break;
      case 'q': sok     = 1; break;
      case 'f': fok     = 1; break;

      case '2': mom     = 2.0; break;
      case 'M': getnum("%lg",&mom); break;
      case 'Z': cross   = 1; break;
      case 'a': is_global_ave = 0; break;
      case '4': binder  = 1; 
	break;
      case 'X': raw     = 1; break;
      case 'E': 
	eig = 1; 
	if (*ss == 'C') { ss++; eig = -1; } 
	break;

      case 'c': tc      = 0; noauto = 1; break;

      case 'd': getnum("%d",&dnum); nd = 1; break;
      case 'D': getlist(lists,nd); 
	if (nd > MAX_DAT) {
	  printf("Too many data elements! max \n", MAX_DAT); 
	  exit(0);
	}
	break;

      case 'p': is_cut = -1; 
	getnum("%lg",&cut_value); 
	cut_v1 = cut_v2 = cut_value;
	break;
      case 'P': is_cut =  1; 
	getnum("%lg",&cut_value);
	cut_v1 = cut_v2 = cut_value;
	break;
      case 'C': is_cut = delta = 1; 
	getnum("%lg",&cut_value);
	cut_v1 = cut_v2 = cut_value;
	break;
      case 'G': is_cut = delta = 1; 
	if (!(*ss) && --argc) ss = (++argv)[0];
	if (sscanf(ss,"%lg:%lg",&cut_v1,&cut_v2) != 2) halt(usage,NULL);
	ss = strchr(ss,0);
	break;
	
      case 'Q': deltaF = 1; break;

      case 'S':
	init = 1;
	if (*ss == '+') { init = -1; ss++; }
	break;

      case 'b':
	if (!(*ss) && --argc) ss = (++argv)[0];
	if (sscanf(ss,"%lg:%lg:%lg",&bmin,&bstep,&bmax) != 3) {
	  if (sscanf(ss,"%lg",&bmin) != 1) halt(usage,NULL);
	  betav[0] = bmin;
	  betanumber = 1; 
	  bmax = bmin;
	  bstep = 1;
	} else {
	  betanumber = 0;
	  while (bmin <= bmax) {
	    beta2v[betanumber] = 0;
	    betav[betanumber++] = bmin;
	    bmin += bstep;
	    if (betanumber > BETAMAX) halt(" + Too long beta-vector (max %d)\n",BETAMAX);
	  }
	}
	ss = strchr(ss,0);
	break;

      case 'z': 
      case 'B': 
	if (*(ss-1) == 'z') { ispart = 1; pt = imbetav; }
	else               { isact2 = 1; pt = beta2v;  }
	if (betanumber <= 0) halt("Give first primary beta (-b)\n",NULL);
	if (!(*ss) && --argc) ss = (++argv)[0];
	if (sscanf(ss,"%lg:%lg:%lg",&bmin,&bstep,&bmax) != 3) {
	  if (sscanf(ss,"%lg",&bmin) != 1) halt(usage,NULL);
	  for (i=0; i<betanumber; i++) pt[i] = bmin;
	} else if (betanumber > 1) {
	  /* Here's the case where we have 2 vectors of betas */
	  for (i=0; i<betanumber; i++,bmin+=bstep) pt[i] = bmin;
	  if (bmin < bmax || bmin > bmax-bstep) 
	    fprintf(stderr," + NOTE: beta2-vector adjusted to end at %g\n",
		    pt[betanumber-1]);
	} else {
	  betanumber = 0;
	  while (bmin <= bmax) {
	    pt[betanumber] = bmin;
	    betav[betanumber++] = betav[0];
	    bmin += bstep;
	    if (betanumber > BETAMAX) halt(" + Too long beta-vector (max %d)\n",BETAMAX);
	  }
	}
	ss = strchr(ss,0);
	break;

      case 's':
	if (!(*ss) && --argc) ss = (++argv)[0];
	if (!*ss) halt(usage,NULL);
	strcpy(npus,ss);
	ss = strchr(ss,0);
	break;

      case 'J':
	jprint = 1;
	if (!(*ss) && --argc) ss = (++argv)[0];
	if (!*ss) halt(usage,NULL);
	if (*ss == '=') jstdout = 1;
	strcpy(js,ss);
	ss = strchr(ss,0);
	break;

      case 'h': ishg = 1; break;
      case 'g': ishg = hg2 = 1; break;
      case 'V': voldiv = 1; break;

      case 'A':
	if (!(*ss) && --argc) ss = (++argv)[0];
	i = sscanf(ss,"%d,%d",&bins,&bins2);
	if (i == 1) bins2 = 0;
	else if (i != 2) halt(usage,NULL);
	ss = strchr(ss,0);
	break;

      case 'o': hp_low = 1; getnum("%lg",&hprint_low); break;
      case 'O': hp_high =1; getnum("%lg",&hprint_high); break;
      case 'W': eq_weight = 1;
	if (*ss == '+') {ss++; b2search=1;} 
	else if (*ss == 'b') {ss++; b2search = -1;}
	else if (*ss == 'B') {ss++; b2search = -2;}
	break;
      case 'H': eq_height = print_ext = 1; 
	if (*ss == '+') {ss++; b2search = 1;}
	else if (*ss == 'b') {ss++; b2search = -1;}
	else if (*ss == 'B') {ss++; b2search = -2;}
	break;
      case 'R': print_ext = 1; break;
      case 'F': flat = 1;      break;
      case 'L': print_ext = 1;
	if (!*ss && --argc) ss = (++argv)[0];
	if (sscanf(ss,"%d,%d,%d",&max_r1,&min_r,&max_r2) != 3) 
	  halt(usage,NULL);
	ss = strchr(ss,0);
	break;
	
      case 't': getnum("%d",&thermo); break;

      case 'j': getnum("%d",&jack); break;

      case 'l': getnum("%d",&maxl); break;

      case 'r': getnum("%lg",&epsilon); break;

      case 'I':
	if (!*ss && --argc) ss = (++argv)[0];
	if (sscanf(ss,"%d:%d",&n1,&n2) != 2) halt(usage,NULL);
	ss = strchr(ss,0);
	break;

      default : halt(usage,NULL);

      } /* switch */
    }

  } /* while */

  data = dnum || fok || sok || nd;

  if (verbose) fprintf(stderr," -- ** -- FS analysis program:\n");

  if ((ishg || data) && betanumber <= 0) halt("Beta value(s) missing\n",NULL);
  if (deltaF && !is_cut) halt(usage,NULL);
  if (argc == 0) halt(usage,NULL);
  if (hg2 && nd != 2) halt("2-histogram requires 2 data vectors\n",NULL);
  if (cross && nd != 2) halt("Cross correlation requires 2 data vectors\n",NULL);


  for (i=0; i<nd; i++) dv[i] = (double **)calloc(MAXRUN,sizeof(double *));

  ss = argv[0];
  strcpy(name,ss);

  if ((dataf = fopen(name,"r")) == NULL) 
    halt("could not open datafile \'%s\'\n",name);

  if (npus[0]) ss = npus;

  {
    char bb[RW_STRLEN+1],*buf;
    int i;
    double g2weak;

    fgets(bb,RW_STRLEN-1,dataf);
    buf = bb; scan_whitespace(buf);

    /******* HERE ARE SECTIONS OF TAILORED CODE.  NOT THE MOST ELEGANT WAY ******/
    if (sscanf(buf,"higgs %lg %lg",&betag,&mh) == 2) {      
      if (verbose)fprintf(stderr," -- SU2-Higgs, Beta_g %g M_H %g\n",betag,mh);
      Higgs = 1;
    } else if (sscanf(buf,"susy23 %lg %lg %lg",&betag,&susy23_mh,&mass_U) == 3) {
      if (verbose) fprintf(stderr,
	  " -- SUSY2+3 -analysis, betag = %g, mh = %g, mU = %g\n",betag,susy23_mh,mass_U);
      susy23_amul = pow(4.0/(betag * 0.42),3.0);
      is_susy = is_susy23 = 1;
      reweight_mult = reweight_data = 1;
    } else if (sscanf(buf,"susy125 %lg",&betag) == 1) {
      if (verbose) fprintf(stderr,
          " -- SUSY125 -analysis, betag = %g\n",betag);
      susy23_mh = 125;
      susy23_amul = pow(4.0/(betag * 0.418),3.0);
      is_susy = is_susy23 = 1;
      reweight_mult = reweight_data = 1;
      fprintf(stderr,"Note: For NEW runs, measure labelling changed!\n");
      
    } else if (sscanf(buf,"susy %lg",&mass_U) == 1) {
      if (verbose) fprintf(stderr," -- SUSY -analysis, mU = %lg\n",mass_U);
      is_susy = 1;
      reweight_mult = reweight_data = 1;
    } else if (sscanf(buf,"su3h %lg",&betag) == 1) {
      if (verbose)
	fprintf(stderr," -- SU3+adjoint Higgs -analysis, betag = %g\n",betag);
      is_su3h = 1;
      isact2 = 1;
      reweight_mult = reweight_data = 1;
    } else if (sscanf(buf,"su2u1 %lg %lg",&betag,&p_x) == 2) {
      if (verbose) 
	fprintf(stderr," -- su2 x u1 -analysis, x = %lg  bg = %g\n",p_x,betag);
      is_u1su2 = 1;
      reweight_data = 9;
      reweight_mult = 1;
    } else {
      /* first, check if there are 2 beta's per file */
      scan_whitespace(buf);
      if (strncmp(buf,"b2",2) == 0) { 
	read_2_beta = 1; 
	buf += 2; 
	scan_whitespace(buf);
      } 
      else read_2_beta = 0;

      if (strncmp(buf,"ascii",5) == 0) {
	is_ascii = reweight_cmd = reweight_mult = 1; 
	scan_not_whitespace(buf); scan_whitespace(buf);
	if (verbose) fprintf(stderr," -- formatted data\n");
	h.lx = h.ly = h.lz = h.lt = 1;
      } else if (strncmp(buf,"reweight",8) == 0) {
	reweight_cmd = reweight_mult = 1;
	scan_not_whitespace(buf); scan_whitespace(buf);
      } else reweight_cmd = 0;

      if (reweight_cmd) {

	/* now read the reweight cmd strings */

	for (i=0; *buf && *buf != '\n' && *buf != ';'; i++, buf++) 
	  reweight_string[i] = *buf;
	reweight_string[i] = 0;

	if (read_2_beta) {
	  if (*buf != ';') 
	    halt("Expecting 2 reweight commands, separated by ';'\n");
	  buf++; scan_whitespace(buf);
	  for (i=0; *buf && *buf != '\n' && *buf != ';'; i++, buf++) 
	    reweight_string2[i] = *buf;
	  reweight_string2[i] = 0;
	  if (betanumber > 0 && !isact2) 
	    halt("Must give 2 beta-values (both -b and -B)\n");
	}
	
	if (verbose) {
	  fprintf(stderr," -- Reweighting with %s",reweight_string);
	  if (read_2_beta) fprintf(stderr," ; %s",reweight_string2);
	  fprintf(stderr,"\n");
	}

	if (is_ascii) {
	  if (sscanf(buf," volume %ld %ld %ld %ld",
		     &h.lx,&h.ly,&h.lz,&h.lt) != 4) {
	    h.lt = 1;
	    if (sscanf(buf," volume %ld %ld %ld",&h.lx,&h.ly,&h.lz) != 3)
	      h.lx = h.ly = h.lz = 1;
	  }
	}

      } else {
	sscanf(buf,"%d %lg",&reweight_data,&reweight_mult);
	if (verbose) fprintf(stderr," -- Reweighting data %d\n",reweight_data);
	if (verbose) fprintf(stderr," -- Multiplication factor %g\n",reweight_mult);
      }
    }
  }

  if (betanumber > 1 && ishg) {
    if (verbose) fprintf(stderr,"Calculating histogram vector\n");
    hgvec = 1;
  }

  if (is_susy) for (i=0; i<betanumber; i++) betav[i] = 1.0/sqr(betav[i]);
  if (is_u1su2 && !isact2) {
    for (i=0; i<betanumber; i++) beta2v[i] = p_x;
    p_x_global = p_x;
  }

  if (reweight_cmd && isact2 && !read_2_beta)
    halt("Cannot do -b and -B with this command file");

  /* if (Higgs && isact2) {
   *   for (i=0; i<betanumber; i++) beta2v[i] = betar(beta2v[i],betav[i]);
   * } 
   */

  i = 0;
  n = -1;
  while (1) {
    int weight;
    char buf[500];
    double b,b2;

    i++;

    if (i > n2) break;

    if (!read_2_beta) {
      if (fscanf(dataf,"%lg %d %s",&b,&nmeas,buf) != 3) break;
    } else {
      if (fscanf(dataf,"%lg %lg %d %s",&b,&b2,&nmeas,buf) != 4) break;
    }

    if (strcmp(buf,"w") == 0) {
      weight = 1; 
      fscanf(dataf," %s",hname);
    } else {
      weight = 0;
      strcpy(hname,buf);
    }

    chp = strrchr(hname,'/');
    if (chp == NULL) chp = hname;
    if (verbose) {
      if (is_susy) 
	fprintf(stderr,"%2d: T %lg measurements %d ",i,b,nmeas);
      else {
	if (!read_2_beta) fprintf(stderr,"%2d: beta %g measurements %d ",i,b,nmeas);
	else fprintf(stderr,"%2d: beta %g, %g measurements %d ",i,b,b2,nmeas);
	if (weight) fprintf(stderr,"+weight");
      }
      fprintf(stderr,"  ..%s\n",chp);
    }

    if (i < n1) {
      if (verbose) fprintf(stderr," skipping this\n");
      continue;
    }
    n++;

    /*    if (n > 0 && b < be[n-1])
      fprintf(stderr," * Warning: vector %d not in ordered place\n",n);
    */
    /* now open histogram file */

    datdim = inifil(hname,&inv_vol);
    if (verbose) fprintf(stderr," - Reading in %d variables\n",datdim);

    if (!Higgs && (datdim < reweight_data)) 
      halt("** illegal reweight-data value\n",NULL);

    /* thermalizing section */

    for (j=0; j < thermo; j++) getdbldat(tmp);

    if (is_susy) b = 1.0/sqr(b);  

    /* Convert T -> 1/T^2, according to
     * p_y = 0.548 - 0.379*sqr(100.0/T);
     * p_Y = 0.517 - 0.849*sqr(mU/T);
     */

    ntop = (nmeas-thermo);
    ap = act[n] = dblarr(ntop);
    if (isact2) a2p = act2[n] = dblarr(ntop);
    if (nd) {
      for (j=0; j<nd; j++) dv[j][n] = dblarr(ntop);
      dp = dv[0][n];
    }
    w[n] = dblarr(ntop);
    nv[n]   = ntop;
    beta[n] = b;
    if (read_2_beta) {
      beta2[n] = b2;
    } else if (Higgs && isact2) { 
      beta2[n] = beta2_Higgs(mh,b);
    } else if (is_u1su2) {
      beta[n] = beta_u1su2(b,p_x); 
      beta2[n] = p_x;
    } else if (is_su3h) {
      beta[n] = beta_su3h(b, &beta2[n]);
    } else beta2[n] = 0;

    /*    betar = sqr(betah)/betag *
     * ((1./8.)*sqr(mh/MW) - 3*G2/(128*pi*MD))/(1. - G2/(24*pi*MD));
     * The derivative
     */

    for (j=0; j<ntop; j++) {
      getdbldat(tmp);

      /*	if (verbose) fprintf(stderr,"%d %d %g\n",j,k,tmp[k-1]);
       */      
      /*      ap[j] = tmp[0];
       *      if (volmul) ap[j] = multiply*6*(1 - ap[j]/3.0);
       */

      /*** NOTE THE SIGN OF ap[j] :  it has to be exp(-b*ap[j])  ****/

      if (Higgs) {
	/*** NOTE: the old bug strikes again: divide tmp[1+w] by two ***/
	
	if (!isact2) ap[j] = tmp[1+weight]/b + (2/b)*tmp[3+weight];
	else {
	  ap[j]  = tmp[1+weight]/b;
	  a2p[j] = tmp[3+weight]/beta2[n];
	}
	sh = tmp[1+weight] + tmp[2+weight] + tmp[3+weight];

	if (sok) dp[j] = tmp[0+weight]+sh;

      } else if (is_susy) {
	
	if (is_susy23) {
	  if (susy23_mh == 105) {
	    ap[j] = susy23_amul*( 18384.1*tmp[3] - 3984.08*tmp[4] - 2*1191.72*tmp[7] 
				  - 2*96.6867*tmp[8] 
				  - sqr(mass_U)*tmp[14] )/inv_vol;
	  } else if (susy23_mh == 95) {
	    ap[j] = susy23_amul*( 13960*tmp[3] - 3960*tmp[4] - 2*990*tmp[7] + 2*97*tmp[8] 
				  - sqr(mass_U)*tmp[14] )/inv_vol;
	  } else if (susy23_mh == 125) {
            // New set! Note now H1=3+w, H2=6+w, Re=8+w, Im=9+w (not used)
            // U = 11+w
            ap[j] = susy23_amul*( 26504.3   * tmp[3+weight] 
                                  - 4004.3  * tmp[6+weight]
                                  - 2*1481.21 * tmp[8+weight]
                                  - 4957.71 * tmp[11+weight] )/inv_vol;
	  }
          
	} else {
	  ap[j] = -0.379*sqr(100.0)*tmp[8] - 0.849*sqr(mass_U)*tmp[16];
	  
	  /* now b ==- 1/T^2, so that have to multiply datavecs with 
	   * the appropriate factors:
	   *
	   * p_y = 0.548 - 0.379*sqr(100.0/T);
	   * p_Y = 0.517 - 0.849*sqr(mU/T);
	   */
	}

        
      } else if (is_u1su2) {
	
	ap[j] = tmp[9];                    /* requires mapping */
	if (isact2) a2p[j] = tmp[4]/p_x;   /* so mult. by x */

      } else if (is_su3h) {

	ap[j]  = tmp[4]/inv_vol;           /* total a2 and a4 */
	a2p[j] = tmp[6]/inv_vol;

      } else if (reweight_cmd) {      /* now reweighting with command */
	
	ap[j] = calclist(tmp,datdim,reweight_string);
	if (isact2) a2p[j] = calclist(tmp,datdim,reweight_string2);

      } else {             /* generic case -- everything should be so straightforward */

	ap[j] = tmp[reweight_data-1]*reweight_mult;
	if (isact2) a2p[j] = tmp[act2_data-1];

      }

      if (dnum) {
	dp[j] = tmp[dnum-1];
	if (voldiv) dp[j] *= inv_vol;
      } else if (nd) {
	int i;
	for (i=0; i<nd; i++) {
	  dv[i][n][j] = calclist(tmp,datdim,lists[i]);
	  if (voldiv) dv[i][n][j] *= inv_vol;
	}
      }
      if (weight) w[n][j] = tmp[0]; else w[n][j] = 0;

    } /* for over lines */

    /* close histogram file */

    if (raw && mom) {
      int i;
      for (i=0; i<nd; i++) {
	dp = dv[i][n];
	for (dum=j=0; j<ntop; j++) dum += dp[j];
	dum /= ntop;
	for (j=0; j<ntop; j++) dp[j] = pow(dp[j]-dum,mom);
      }
    }

  } /* while (1) */

  if (verbose && h.lx != 1) 
    fprintf(stderr," system size: %dx%dx%d\n",(int)h.lx,(int)h.ly,(int)h.lz);

  n++;

  if (raw) {

    /* now calculate only raw data */
    for (i=0; i<nd; i++) errors(n,nv,beta,dv[i],tc,jack,w);
    return(0);
  }

  /* initialize Ferrenberg-Swendsen package */
  
  if (ishg && b2search && !hgvec) betanumber = 1;
  FSinit(n,nv,beta,beta2,act,act2,epsilon,tc,maxl,init,ss,jack,irange,w);

  /* Initialization done, calculate observables */

  if (data || nd) {
    if (ishg && !hgvec) {
      if (hg2) FShgram2(dv[0],dv[1],betav[0],beta2v[0],bins,bins2);
      else {
	if (!is_cut && (eq_weight || eq_height || print_ext)) 
	  halt(" *** cut value not defined!",NULL);
	if (jack) FShgJack(dv,nd,betav[0],beta2v[0],bins,cut_value,jack,
			   jprint,js,ss,eq_weight,eq_height);
	else FShgram(dv,nd,betav[0],beta2v[0],
		     bins,cut_value,eq_weight,eq_height);
	/* FShgram(dv[0],betav[0],beta2v[0],
	   bins,cut_value,eq_weight,eq_height); */
      }
      return(0);
    }
    if (eig && jack) FSeigJack(dv,nd,jack,jprint,js,ss);
    else if (eig)    FSeig(dv,nd,eig<0,lists);
    else if (jack) {
      if (!hgvec) FSJack(dv,nd,jack,fok,jprint,js,ss,mom,binder,cross);
      else FShgvecJ(dv,jack,jprint,js,ss,bins,cut_value,eq_weight,eq_height);
    } else {
      if (!hgvec) FSCalc(dv,nd,fok,mom,binder,cross);
      else FShgvec(dv,bins,cut_value,eq_weight,eq_height);
    }
  }
  else if (ispart) FSCpart(imbetav,betanumber);
}

/*********************************************/

int
errors(int n,int nv[], double beta[], double *a[],int tc,int jack,
       double *w[])
{
  int j,i,k;
  double average,taverage,error,naive,norm,tnorm,timerel,bv, *dat;
  double t;

  bv = 0;
  if (jack) {
    for (j=0; j<n; j++) {
      dat = a[j];
      error = taverage = 0;
      tnorm = 0;
      for (k=0; k<nv[j]; k++) {
	t = exp(-w[j][k]);
	taverage += dat[k]*t;
	tnorm += t;
      }
      for (i=0; i<jack; i++) {
	average = norm = 0;
	for (k=i*nv[j]/jack; k<(i+1)*nv[j]/jack; k++) {
	  t = exp(-w[j][k]);
	  average += dat[k]*t;
	  norm += t;
	}
	error += sqr((taverage-average)/(tnorm-norm))/(jack*(jack-1));
      }
    }
    printf("%lg\t%.14lg\t%.14lg\n",p_beta(beta[j]),taverage/tnorm,sqrt(error));
    return(1);
  }
  for (j=0; j<n; j++) {
    if (a[j] != NULL) {
      errorcalc(nv[j],a[j],&average,&naive,&timerel,&bv,w[j]);
      if (tc) error = naive*(double)sqrt(2.0*fabs(timerel));
      else    error = naive;
      printf("%lg\t%.14lg\t%.14lg\n",p_beta(beta[j]),average,error);
    }
  }
  return(1);
}


/**************************************************************************
 *
 *     Ferrenberg-Swendsen package
 *
 *        Kari Rummukainen Apr 26 1990 - Apr 1993
 *        ref: Ferrenberg, Swendsen: Computers in Physics, sep/oct 1989
 *
 *     this subroutine package initializes the Ferrenberg-Swendsen
 *     multiple histogram analysis
 *
 *     int FSinit(num,nvec,beta,act,eps,tc,maxloop,init,name
 *            betaval,betanum,jack)
 *            returns the number of iterations
 *
 *     int num           number of vectors
 *     int nvec[num]     len of each vector
 *     double beta[num]  beta of each vector
 *     double *act[num]  energy values of the vectors
 *     double eps        relative accuracy of free-energy
 *     int tc            1/0 - calculate time correlations
 *     int maxloop       number of maximum loops; 0-infinity
 *     int init          1/0 - use data / calculate it
 *     double betaval[betanum]  betas for the values needed
 *     int jack          num/0 - use jackknife with blocks/do not use
 *
 *     void FSfree(free)        - calculates the free energy
 *
 *     double free[betan]        free energy values
 *
 *
 *     double FSval(values,sigma,obs)  - calculates the exp.value of observables
 *
 *     double values[betan]      return value
 *     double *obs[num]          values of observables
 *
 *
 *     void FSCalc(double data[betan],"echo","sym","xf.fs",isfree)
 *
 *     -- Jackknife routines 28 Aug 1990
 *
 */

void FSval(double val[],double **obs);
void FSfree(double free[]);
void FSmom(double * val,double **obs,double mom,double *average);
void FScross(double *val,double **obs1,double **obs2,double **average);
void FSBinder(double * val,double **obs,double *average);
void FShg(double *his,double **d[],int nd,double beta,double beta2,
	  double min1,double max1,int bins);
void FShg2(double *hg,double **d1,double **d2,double beta,double beta2,
	   double min1,double max1,int bins1,
	   double min2,double max2,int bins2);
double findhg(double hist[],double **d[],int nv,double beta,double beta2,
	      double min1,double max1,int bins,double cut,
	      int eq_weight,int eq_height);
double findb2hg(double hist[],double **d[],int nd,double beta,double *beta2,
		double min1,double max1,int bins,double cut,
		int eq_weight,int eq_height);
double FSBindersearch(double **dat,double *betain, double beta2);
void FShgvect(double *bv, double *val, double **obs[],
	      int bins,double cut_value,int eq_weight,int eq_height);
void get_minmax(double **d[],int nd, double *min1, double *max1);

void FSinidown();
void FSinimeas();
void FSsafeinit();
double FSlogsum(double beta,double beta2);
complex FSpartition(double beta,double beta2,double im);
int FSini1(int num,int *nvec,double *beta,double *beta2,
	   double *fp,double *gp,double **act,double **act2,
	   int maxl,int jack,double *jf,double *w[]);

int NewtonRaphson(int n,double v[],int maxl,double eps);

#define MAXJACK 100

static double *down[MAXRUN], *maxv[MAXRUN];
static double **a, **wt, **a2;
static double g[MAXRUN], f[MAXRUN], fc[MAXRUN];
static double *b, *b2, *ep, maxp[MAXRUN], logz[MAXRUN];
static int * nn, nv[MAXRUN], *nr;
static double * logsum;
static double b_ind[MAXRUN],b2_ind[MAXRUN],f_ind[MAXRUN];
static int nb_i[MAXRUN],pb_i[MAXRUN][10],ind_i[MAXRUN];
static int n,n_independent;
static double eps;
static int bnum;
static double *akt[MAXRUN],*wkt[MAXRUN],*akt2[MAXRUN];


int
FSinit(int num,int *nvec,double *beta,double *beta2,double **act,double **act2,
       double epsilon,int tc,int maxl,int init,char *ss,
       int jack, int irange,double **w) {

  int   i,j,k,nind,nbind[MAXRUN],ip,ir,nb;
  FILE  * fil;
  double tmp1, *ff;
  double oldact,bind[MAXRUN],b2ind[MAXRUN],fr[MAXRUN],gr[MAXRUN];
  double *jf,*jp;
  char  name[500];

  eps = epsilon;
  nind = 0;

  for (j=0; j<num; j++) {
    i = 1;
    for (k=0; k<nind && i; k++) {
      if (bind[k] == beta[j] && b2ind[k] == beta2[j]) {
	i = 0;
	nbind[k]++;
      }
    }
    if (i) {
      bind[nind]    = beta[j];
      b2ind[nind]   = beta2[j];
      nbind[nind++] = 1;
    }
  }

  /* for (j=0;j<nind; j++)
   * fprintf(stderr,"%lg\n",bind[j]);
   */

  if (verbose)
    fprintf(stderr," %d vectors, where %d independent\n",num,nind);

  if (init) {
    if (verbose)
      if (init > 0)
	fprintf(stderr,"- using old FS info for %d vectors\n",num);
      else fprintf(stderr,
	"- using old FS info as a starting point for %d vectors\n",num);

    strcpy(name,ss);  strcat(name,".values"); fil = fopen(name,"r");
    for (j=0; j<num; j++) {
      if (fscanf(fil,"%d %lg %lg %lg",&i,&tmp1,&fr[j],&gr[j]) != 4) 
	halt(" ** File error in %s\n",name);

      if (fabs(tmp1 - beta[j]) > 1e-13) 
	halt(" ** Beta error, b = %.14lg \n",tmp1);
    }
    fclose(fil);

  } else {
    double d_min,d_max,ave,naiv,integ;

    if (verbose) fprintf(stderr,"- initializing FS with %d vectors\n",num);

    /* calculate the time correlation */

    for (i=0; i<num; i++) {
      errorcalc(nvec[i],act[i],&ave,&naiv,&integ,&tmp1,w[i]);
      gr[i] = 1.0/(2.0*integ);
      if (!tc) gr[i] = 1.;

      if (i) fr[i] = fr[i-1] + 0.5*(ave+oldact)*(beta[i]-beta[i-1]);
      else   fr[i] = 0.0;

      oldact = ave;

      d_min = 1e90; d_max = -d_min;
      for (j=0; j<nvec[i]; j++) {
	d_min = (d_min < act[i][j]) ? d_min : act[i][j];
	d_max = (d_max > act[i][j]) ? d_max : act[i][j];
      }

      naiv *= sqrt((double)nvec[i]);
      if (verbose)
	fprintf(stderr,
		"%d: min: %g  w-ave-w: %g-%g-%g,  max %g  autocorr %g\n"
		,i+1,d_min,ave-naiv,ave,ave+naiv,d_max,integ);
    }
  }

  /** NORMALIZE WEIGHT **/

  for (i=0; i<num; i++) {
    double ws;
    for (ws=j=0; j<nvec[i]; j++) ws += exp(w[i][j]);
    ws = log(ws/nvec[i]);
    for (j=0; j<nvec[i]; j++) w[i][j] -= ws;
  }

  if (init <= 0) {
    if (jack) {
      jf = dblarr(jack*num);
      jp = dblarr(jack*num);
    }

    ff = dblarr(num);
    if (irange && nind > 2) {
      if (verbose) fprintf(stderr,"- calculating f in blocks::\n");
      ip = 0;
      for (ir=0; ir<nind-1; ir++) {
	nb = nbind[ir]+nbind[ir+1];
	if (verbose)
	  fprintf(stderr,"======== block %d, %d vectors ========\n",ir,nb);

	for (i=0; i<nb; i++) ff[i] = fr[i+ip] - fr[ip];

	FSini1(nb,nvec+ip,beta+ip,beta2+ip,ff,gr+ip,
	       act+ip,act2+ip,maxl,jack,jp,w+ip);

	for (i=0; i<nb; i++) fr[i+ip] = fr[ip] + ff[i];

	if (jack) {
	  for (i=0; i<jack; i++) for (j=0; j<nb; j++) {
	    jf[i*num+j+ip] = jp[i*nb+j] + jf[i*num+ip] - jp[i*nb];
	  }
	}
	ip += nbind[ir];
      }

    } else FSini1(num,nvec,beta,beta2,fr,gr,act,act2,maxl,jack,jf,w);

    /* write the files */

    strcpy(name,ss);  strcat(name,".values"); fil = fopen(name,"w");
    for (j=0; j<num; j++) {
      fprintf(fil,"%d %.18lg\t%.20lg \t%.20lg\n",
	      j+1,beta[j],fr[j],gr[j]);
    }
    fclose(fil);

    free(ff);

    if (jack) {
      strcpy(name,ss);  strcat(name,".jack"); fil = fopen(name,"w");
      for (i=0; i<jack; i++) {
	for (j=0; j<num; j++) {
	  fprintf(fil,"%d %.20lg %.20lg\n",j,beta[j],jf[i*num+j]);
	}
	fprintf(fil,"\n");
      }
      fclose(fil);

      free(jp); free(jf);
    }
  }

  if (jack) {
    jack = (jack > MAXJACK) ? MAXJACK : jack;
    for (i=0; i<num; i++) {
      akt[i] = dblarr(nvec[i] - nvec[i]/jack + 1);
      if (isact2) akt2[i] = dblarr(nvec[i] - nvec[i]/jack + 1);
      wkt[i] = dblarr(nvec[i] - nvec[i]/jack + 1);
    }
  }

  /* allocate workspaces */

  j = 0;
  for (i=0; i<num; i++) {
    maxv[i] = dblarr(nvec[i]);
    down[i] = dblarr(nvec[i]);
    j = max(j,nvec[i]);
  }
  ep   = dblarr(j);

  /* init the global variables */

  n   = num;
  nn  = nvec;
  b   = beta;  b2  = beta2;
  a   = act;   a2  = act2;
  wt  = w;
  for (i=0; i<n; i++) {
    nv[i] = nn[i]; fc[i] = f[i] = fr[i]; g[i] = gr[i];
  }

  n_independent = 0;
  for (j=0; j<n; j++) {
    i = 1;
    for (k=0; k<n_independent && i; k++) {
      if (b_ind[k] == beta[j] && b2_ind[k] == beta2[j]) {
	pb_i[k][nb_i[k]++] = j;
	ind_i[j] = k;
	i = 0;
      }
    }
    if (i) {
      b_ind[n_independent]   = beta[j];
      b2_ind[n_independent]  = beta2[j];
      pb_i[n_independent][0] = j;
      ind_i[j] = n_independent;
      nb_i[n_independent++]  = 1;
    }
  }

  FSsafeinit();

  if (betanumber > 0) {
    bnum    = betanumber;
    logsum  = dblarr(bnum+10);

    FSinimeas();

  } else bnum = 0;

  return(j);
}

/**********************************************************/

int
FSini1(int num,int *nvec,double *beta,double *beta2,double *fp,double *gp,
       double **act,double **act2,
       int maxl,int jack,double *jf,double **w) {
  
  int   i,j,k;
  double *ap,*wp,*w0p;
  int   ntop,nbot;
  double * pp;

  /* init global variables... temporarily */

  n   = num;
  nn  = nvec;
  b   = beta;
  b2  = beta2;
  a   = act;
  a2  = act2;
  wt  = w;

  /* allocate workspaces - also temporary */

  j = 0;
  for (i=0; i<num; i++) {
    maxv[i] = dblarr(nvec[i]);
    down[i] = dblarr(nvec[i]);
    j = max(j,nvec[i]);
  }
  ep   = dblarr(j);

  for (i=0; i<n; i++) {
    nv[i] = nn[i]; fc[i] = f[i] = fp[i]; g[i] = gp[i];
  }

  n_independent = 0;
  for (j=0; j<n; j++) {
    i = 1;
    for (k=0; k<n_independent && i; k++) {
      if (b_ind[k] == beta[j] && b2_ind[k] == beta2[j]) {
	pb_i[k][nb_i[k]++] = j;
	ind_i[j] = k;
	i = 0;
      }
    }
    if (i) {
      b_ind[n_independent]   = beta[j];
      b2_ind[n_independent]  = beta2[j];
      pb_i[n_independent][0] = j;
      ind_i[j] = n_independent;
      nb_i[n_independent++]  = 1;
    }
  }

  /* start jackknife routines .. */

  if (jack) {

    jack = (jack > MAXJACK) ? MAXJACK : jack;
    for (i=0; i<n; i++) {
      akt[i] = dblarr(nvec[i] - nvec[i]/jack + 1);
      if (isact2) akt2[i] = dblarr(nvec[i] - nvec[i]/jack + 1);
      wkt[i] = dblarr(nvec[i] - nvec[i]/jack + 1);
    }
    a = akt;
    a2= akt2;
    wt= wkt;

    for (j=0; j<n; j++) fc[j] = 0;

    for (k=0; k<jack; k++) {
      for (i=0; i<n; i++) {
	nbot = k*nn[i]/jack;
	ntop = (k+1)*nn[i]/jack;
	nv[i] = nn[i] - ntop + nbot;
	ap = act[i];
	pp = akt[i];
	for (j=0;    j<nbot;  j++) pp[j]           = ap[j];
	for (j=ntop; j<nn[i]; j++) pp[j-ntop+nbot] = ap[j];
	wp = wkt[i];
	w0p= w[i];
	for (j=0;    j<nbot;  j++) wp[j]           = w0p[j];
	for (j=ntop; j<nn[i]; j++) wp[j-ntop+nbot] = w0p[j];
	if (isact2) {
	  ap = act2[i];
	  pp = akt2[i];
	  for (j=0;    j<nbot;  j++) pp[j]           = ap[j];
	  for (j=ntop; j<nn[i]; j++) pp[j-ntop+nbot] = ap[j];
	}
      }

      FSsafeinit();

      for (j=0; j<n_independent; j++) f_ind[j] = f[pb_i[j][0]];
      j = NewtonRaphson(n_independent,f_ind,maxl,eps);
      if (verbose)
	fprintf(stderr," -- block %d, %d Newton-Raphson iterations\n",k,j);
      for (j=0;j<n_independent;j++)
	for (i=0;i<nb_i[j];i++) fp[ pb_i[j][i] ] = f[ pb_i[j][i] ] = f_ind[j];

      for (j=0; j<n; j++) {
	fc[j] += f[j];
	jf[k*n+j] = f[j];
      }
    }

    for (i=0; i<n; i++) {
      fp[i] = fc[i] /= jack;
      nv[i] = nn[i];
    }

    a = act;
    a2= act2;
    wt= w;
    for (i=n-1; i>=0; i--) {
      free(wkt[i]); free(akt[i]); if (isact2) free(akt2[i]);
    }

  } else {

    /* calculate the appropriate f */

    FSsafeinit();
    for (j=0; j<n_independent; j++) f_ind[j] = f[pb_i[j][0]];
    j = NewtonRaphson(n_independent,f_ind,maxl,eps);
    if (verbose) fprintf(stderr," .. %d Newton-Raphson iterations\n",j);
    for (j=0;j<n_independent;j++)
      for (i=0;i<nb_i[j];i++) fp[ pb_i[j][i] ] = f_ind[j];

    if (verbose) for (j=0; j<n; j++)
      fprintf(stderr,"  f:%g g:%g b:%g nv:%d\n",f[j],g[j],p_beta(b[j]),nv[j]);
  }

  free(ep);
  for (i=num-1; i>=0; i--) {
    free(down[i]); free(maxv[i]);
  }
  return(1);
}

/*******************************************************************/

void 
FSsafeinit() 
{
  int i,j,k;
  double *dp, *ap, *a2p, *wp;
  double *mxv;
  double val,aval;
  int end;

  if (verbose) { fprintf(stderr," + safe init .."); fflush(stderr); }

  for (k=0; k<n; k++) {
    ap  = a[k];
    a2p = a2[k];
    mxv = maxv[k];
    end = nv[k];
    dp  = down[k];
    wp  = wt[k];

    /* first phase - calculate the maximum of b*act */

    /* this just initializes ... */
#pragma ivdep
    for (j=0; j<end; j++) {
      dp[j]  = 0.0;
      mxv[j] = -b[k]*ap[j] + f[k];
    }
    if (isact2) for (j=0; j<end; j++) mxv[j] -= b2[k]*a2p[j];

    for (i=0; i<n; i++) if (i != k) {
      if (!isact2) 
#pragma ivdep
	for (j=0; j<end; j++) mxv[j] = max(mxv[j],-b[i]*ap[j] + f[i]);
      else 
#pragma ivdep
	for (j=0; j<end; j++) 
	  mxv[j] = max(mxv[j],-b[i]*ap[j] - b2[i]*a2p[j] + f[i]);
    }

    /* second - initialize down-array */
    
    for (i=0; i<n; i++) {
      val = log(nv[i]*g[i]) + f[i];
      if (!isact2)
#pragma ivdep
	for (j=0; j<end; j++) dp[j] += exp(-b[i]*ap[j] + wp[j] + val - mxv[j]);
      else
#pragma ivdep
 	for (j=0; j<end; j++) 
	  dp[j] += exp(-b[i]*ap[j] - b2[i]*a2p[j] + wp[j] + val - mxv[j]);
    }
#pragma ivdep
    for (j=0; j<end; j++) {
      if (dp[j] <= 0.0) 
	fprintf(stderr," ** down %lg, v %d  i %d\n",dp[j],k,j);
      dp[j] = log(dp[j]) + mxv[j];
    }

  }

  /* third - calculate the maximum of b*act - down : CHANGED TO n_ind */

  for (i=0; i<n_independent; i++) {
    maxp[i] = -1e100;

    for (k=0; k<n; k++) {
      end = nv[k];
      ap  = a[k]; 
      a2p = a2[k];
      dp  = down[k];
      wp  = wt[k];

      if (!isact2) 
#pragma ivdep
	for (j=0; j<end; j++) ep[j] = -b_ind[i]*ap[j] - dp[j];
      else 
#pragma ivdep
	for (j=0; j<end; j++) ep[j] = -b_ind[i]*ap[j] - b2_ind[i]*a2p[j] - dp[j];

      for (j=0; j<end; j++) maxp[i] = max(maxp[i],ep[j]);
    }
  }

  /* fourth - calculate logz-array : CHANGED TO n_ind */

  for (k=0; k<n_independent; k++) {
    val = 0.0;
    for (i=0; i<n; i++) {
      ap  = a[i];
      a2p = a2[i];
      wp  = wt[i];
      dp  = down[i];
      aval = log(g[i]) - maxp[k];

      if (!isact2)
#pragma ivdep
	for (j=0; j<nv[i]; j++) val += exp(-b_ind[k]*ap[j] - dp[j] + aval);
      else 
#pragma ivdep
	for (j=0; j<nv[i]; j++) 
	  val += exp(-b_ind[k]*ap[j] - b2_ind[k]*a2p[j] - dp[j] + aval);
	
      /*  fprintf(stderr,"k=%d i=%d val=%lg\n",k,i,val);
       */
    }
    logz[k] = log(val) + maxp[k];
    /* fprintf(stderr,"logz %d : %lg\n",k,logz[k]);
     */
  }

  if (verbose) fprintf(stderr,". done\n");

}


/**************/


void 
FSinidown() 
{
  int i,j,k;
  double *dp, *ap, *a2p, *wp;
  double *mxv;
  double val,aval;
  int end;


  for (k=0; k<n; k++) {
    ap  = a[k];
    a2p = a2[k];
    wp  = wt[k];
    mxv = maxv[k];
    end = nv[k];
    dp  = down[k];

#pragma ivdep
    for (j=0; j<end; j++) dp[j] = 0.0;

    for (i=0; i<n; i++) {
      val = log(nv[i]*g[i]) + f[i];

      if (!isact2)
#pragma ivdep
	for (j=0; j<end; j++) dp[j] += exp(-b[i]*ap[j] + val - mxv[j] + wp[j]);
      else 
#pragma ivdep
	for (j=0; j<end; j++) 
	  dp[j] += exp(-b[i]*ap[j] - b2[i]*a2p[j] + val - mxv[j] + wp[j]);
    }
#pragma ivdep
    for (j=0; j<end; j++) dp[j] = log(dp[j]) + mxv[j];

  }

  for (k=0; k<n_independent; k++) {
    val = 0.0;
    for (i=0; i<n; i++) {
      ap = a[i];
      a2p= a2[i];
      wp = wt[i];
      dp = down[i];
      aval = log(g[i]) - maxp[k];

      if (!isact2)
#pragma ivdep
	for (j=0; j<nv[i]; j++) 
	  val += exp(-b_ind[k]*ap[j] - dp[j] + aval);
      else
#pragma ivdep
	for (j=0; j<nv[i]; j++) 
	  val += exp(-b_ind[k]*ap[j] - b2_ind[k]*a2p[j] - dp[j] + aval);
    }
    logz[k] = log(val) + maxp[k];
  }
}


/**************/

void 
FSinimeas()
{
  int j;

  for (j=0; j<bnum; j++) logsum[j] = FSlogsum(betav[j],beta2v[j]);
}

/**************/

double 
FSlogsum(double beta,double beta2)
{
  int i,j,end;
  double *ap, *a2p, *dp, val, mval;

  val  = 0.0;
  mval = -1e100;

  map_beta(&beta,&beta2);

  for (i=0;i<n;i++) {
    ap  = a[i];
    a2p = a2[i];
    dp  = down[i];
    end = nv[i];

    if (!isact2)
#pragma ivdep
      for (j=0; j<end; j++) {
	ep[j] = -beta*ap[j] - dp[j];
	mval = max(mval,ep[j]);
      }
    else
#pragma ivdep
      for (j=0; j<end; j++) {
	ep[j] = -beta*ap[j] - beta2*a2p[j] - dp[j];
	mval = max(mval,ep[j]);
      }
  }

  for (i=0;i<n;i++) {
    ap = a[i];
    a2p= a2[i];
    dp = down[i];

    if (!isact2)
#pragma ivdep
      for (j=0;j<nv[i];j++) val += g[i]*exp(-beta*ap[j] - dp[j] - mval);
    else
#pragma ivdep
      for (j=0;j<nv[i];j++) 
	val += g[i]*exp(-beta*ap[j] - beta2*a2p[j] - dp[j] - mval);
  }
  return(log(val) + mval);
}


/**************/

void 
FSval(double * val,double ** obs)
{
  int i,j,k;
  double vp,*wsum,*wacc1,*wacc2,*val1,*val2;
  double *ap, *a2p, *dp, *op, dval, be, b2e, v1, v2;

#pragma ivdep
  for (k=0; k<bnum; k++) val[k] = 0.0;

  if (is_cut) {
    wsum  = dblarr(bnum);
    wacc1 = dblarr(bnum);
    wacc2 = dblarr(bnum);
    val1  = dblarr(bnum);
    val2  = dblarr(bnum);
  }

  for (i=0; i<n; i++) {
    ap = a[i];
    a2p= a2[i];
    dp = down[i];
    op = obs[i];

    for (k=0; k<bnum; k++) {
      be = betav[k];
      b2e= beta2v[k];
      map_beta(&be,&b2e);

      dval = logsum[k] - log(g[i]);
      v1 = v2 = 0.;

      if (!is_cut) {
	if (!isact2)
#pragma ivdep
	  for (j=0;j<nv[i];j++) v1 += op[j]*exp(-be*ap[j] - dp[j] - dval);
	else
#pragma ivdep
	  for (j=0;j<nv[i];j++) 
	    v1 += op[j]*exp(-be*ap[j] - b2e*a2p[j] - dp[j] - dval);
	val[k] += v1;
      } else {
	for (j=0; j<nv[i]; j++) {
	  if (!isact2) vp = exp(-be*ap[j] - dp[j] -dval);
	  else vp = exp(-be*ap[j] - b2e*a2p[j] - dp[j] -dval);
	  wsum[k] += vp;
	  if (op[j] > cut_v2) {
	    v1 += op[j]*vp;
	    wacc1[k] += vp;
	  } else if (op[j] <= cut_v1) {
	    v2 += op[j]*vp;
	    wacc2[k] += vp;
	  }
	}
	val1[k] += v1;
	val2[k] += v2;
      }
    }
  }

  if (is_cut) {
    if (deltaF) for (k=0; k<bnum; k++)
      val[k] = (log(wacc1[k]) - log(wacc2[k]))*inv_vol;
    else if (delta) for (k=0; k<bnum; k++) 
      val[k] = val1[k]*wsum[k]/wacc1[k] - val2[k]*wsum[k]/wacc2[k];
    else if (is_cut > 0) 
      for (k=0; k<bnum; k++) val[k] = val1[k]*wsum[k]/wacc1[k];
    else 
      for (k=0; k<bnum; k++) val[k] = val2[k]*wsum[k]/wacc2[k];

    free(wsum);
    free(wacc1);
    free(wacc2);
    free(val1);
    free(val2);
  }
}

/**************/

void 
FSCalc(double **dv[],int nd,int isfree,double mom,int binder,int cross)
{
  int i,j,k,l,nv,ntop,nbot;
  double **val;

  if (isfree) nd = nv = 1;
  else if (cross) nv = 1;
  else nv = nd;
  val = (double **)calloc(nv,sizeof(double *));
  for (i=0; i<nv; i++) val[i] = dblarr(bnum);
  
  if (cross) FScross(val[0],dv[0],dv[1],(double **)NULL);
  else for (i=0; i<nd; i++) {
    if (mom) FSmom(val[i],dv[i],mom,(double *)NULL);
    else if (binder) FSBinder(val[i],dv[i],(double *)NULL);
    else if (isfree) FSfree(val[i]);
    else {FSval(val[i],dv[i]); /* FSerr(err,dat); */}
  }
  for (i=0; i<bnum; i++) {
    fprintf(stdout,"%.12lg  ",p_beta(betav[i]));
    if (isact2 && !is_su3h) fprintf(stdout,"%.12lg  ",(beta2v[i]));
    for (j=0; j<nv; j++) fprintf(stdout,"%.14lg  ",val[j][i]);
    fprintf(stdout,"\n");
  }

  for (i=0; i<nv; i++) free(val[i]);
  free(val);
}

/**************/

void 
FSmom(double *val,double **obs,double mom,double *average) 
{
  int i,j,k,nn,cutbuf;
  double *ap,*a2p,*dp,*op,beta,beta2,dval,s1,*ave;

  for (k=0; k<bnum; k++) val[k] = 0.0;

  cutbuf = is_cut; is_cut = 0;
  if (average == NULL) {
    ave = dblarr(bnum);
    FSval(ave,obs);
  } else ave = average;

  is_cut = cutbuf;

  for (i=nn=0; i<n; i++) {
    ap = a[i];
    a2p= a2[i];
    dp = down[i];
    op = obs[i];

    nn += nv[i];
    for (k=0; k<bnum; k++) {
      beta = betav[k];
      beta2= beta2v[k];
      map_beta(&beta,&beta2);
      dval = logsum[k] - log(g[i]);

      s1 = 0.;

      if (mom == 2.0) {
	if (!isact2) {
#pragma ivdep
	  for (j=0;j<nv[i];j++) {
	    s1 += sqr(op[j] - ave[k])*exp(-beta*ap[j] - dp[j] - dval);
	  } 
	} else {
#pragma ivdep
	  for (j=0;j<nv[i];j++) {
	    s1 += sqr(op[j] - ave[k])*
	      exp(-beta*ap[j] -beta2*a2p[j] -dp[j] - dval);
	  }
	}
      } else {
	if (!isact2) { 
#pragma ivdep
	  for (j=0;j<nv[i];j++) {
	    s1 += pow(op[j] - ave[k],mom)*
	      exp(-beta*ap[j] - dp[j] - dval);
	  } 
	} else { 
#pragma ivdep
	  for (j=0;j<nv[i];j++) {
	    s1 += pow(op[j] - ave[k],mom)*
	      exp(-beta*ap[j] -beta2*a2p[j] -dp[j] - dval);
	  }
	}
      }
      val[k] += s1;
    }
  }
  if (average == NULL) free(ave);
}

/**************/

void 
FScross(double *val,double **obs1,double **obs2,double **average)
{
  int i,j,k,nn,cutbuf;
  double *ap,*a2p,*dp,beta,beta2,dval,s1,*ave1,*ave2;

  for (k=0; k<bnum; k++) val[k] = 0.0;

  cutbuf = is_cut; is_cut = 0;

  /* calculate average here only if desired */
  if (average == NULL) {
    ave1 = dblarr(bnum);
    ave2 = dblarr(bnum);
    FSval(ave1,obs1);
    FSval(ave2,obs2);
  } else {
    ave1 = average[0];
    ave2 = average[1];
  }
  is_cut = cutbuf;

  for (i=nn=0; i<n; i++) {
    ap = a[i];
    a2p= a2[i];
    dp = down[i];

    nn += nv[i];
    for (k=0; k<bnum; k++) {
      beta = betav[k];
      beta2= beta2v[k];
      map_beta(&beta,&beta2);
      dval = logsum[k] - log(g[i]);

      s1 = 0.;

      if (!isact2) {
#pragma ivdep
	for (j=0;j<nv[i];j++) {
	  s1 += (obs1[i][j] - ave1[k])*(obs2[i][j] - ave2[k])*
	    exp(-beta*ap[j] - dp[j] - dval);
	} 
      } else {
#pragma ivdep
	for (j=0;j<nv[i];j++) {
	  s1 += (obs1[i][j] - ave1[k])*(obs2[i][j] - ave2[k])*
	    exp(-beta*ap[j] -beta2*a2p[j] -dp[j] - dval);
	}
      }
      val[k] += s1;
    }
  }
  if (average == NULL) {
    free(ave1); free(ave2);
  }
}

/**************/

void 
FSBinder(double * val,double **obs,double *average) 
{
  int i,j,k,nn,cutbuf;
  double *ap,*ave,*a2p,*dp,*op,beta,beta2,dval,s2,s4,od,ob,*a2t;

  a2t = dblarr(bnum);

  for (k=0; k<bnum; k++) val[k] = 0.0;

  cutbuf = is_cut; is_cut = 0;

  if (average == NULL) {
    ave = dblarr(bnum);
    FSval(ave,obs);
  } else ave = average;

  is_cut = cutbuf;

  for (i=nn=0; i<n; i++) {
    ap = a[i];
    a2p= a2[i];
    dp = down[i];
    op = obs[i];

    nn += nv[i];
    for (k=0; k<bnum; k++) {
      beta = betav[k];
      beta2= beta2v[k];
      map_beta(&beta,&beta2);
      dval = logsum[k] - log(g[i]);

      s2 = s4 = 0.;

#pragma ivdep
      for (j=0;j<nv[i];j++) {
	od = sqr(op[j]-ave[k]);
	if (!isact2) ob = exp(-beta*ap[j] - dp[j] - dval);
        else ob = exp(-beta*ap[j] - beta2*a2p[j] - dp[j] - dval);
	s2 += od*ob;
	s4 += sqr(od)*ob;
      }
      val[k] += s4;
      a2t[k] += s2;
    }
  }

  for (k=0; k<bnum; k++) val[k] = 1.0 - val[k]/(3.0*sqr(a2t[k]));

  free(a2t);
  if (average == NULL) free(ave);
}


/**************
 *
 * calculate the partition function with a complex
 * coupling
 */

complex
FSpartition(double beta,double beta2,double im)
{
  int i,j,end;
  double *ap, *a2p, *dp, mval,t;
  complex val;

  val.r = val.i  = 0.0;
  mval = -1e100;

  map_beta(&beta,&beta2);

  for (i=0;i<n;i++) {
    ap  = a[i];
    a2p = a2[i];
    dp  = down[i];
    end = nv[i];

    if (!isact2)
#pragma ivdep
      for (j=0; j<end; j++) {
	ep[j] = -beta*ap[j] - dp[j];
	mval = max(mval,ep[j]);
      }
    else
#pragma ivdep
      for (j=0; j<end; j++) {
	ep[j] = -beta*ap[j] - beta2*a2p[j] - dp[j];
	mval = max(mval,ep[j]);
      }
  }

  for (i=0;i<n;i++) {
    ap = a[i];
    a2p= a2[i];
    dp = down[i];

    if (!isact2)
#pragma ivdep
      for (j=0;j<nv[i];j++) {
	t = g[i]*exp(-beta*ap[j] - dp[j] - mval);
	val.r += t * cos( im * ap[j] );
	val.i += t * sin( im * ap[j] );
      }
    else
#pragma ivdep
      for (j=0;j<nv[i];j++) {
	t = g[i]*exp(-beta*ap[j] - beta2*a2p[j] - dp[j] - mval);
	val.r += t * cos( im * ap[j] );
	val.i += t * sin( im * ap[j] );
      }
  }

  if (mval > 200) mval = 200;
  val.r *= exp(mval);
  val.i *= exp(mval);

  return(val);
}


/*********************************************/

void
FShgvect(double *bv, double *val, double **dv[],
	 int bins,double cut_value,int eq_weight,int eq_height)
{
  int i;
  double *his;
  double min1,max1,estep;
  double maxv1,maxv2,minv,maxl1,maxl2,minl,cu1,cu2;

  if (bins == 0) bins = 100;

  his = dblarr(bins);

  get_minmax(dv,1,&min1,&max1);
  if (verbose) fprintf(stderr,"  Min: %lg  Max: %lg  Bins: %d\n",min1,max1,bins);

  /* modify min and max a bit ... */
  
  estep = (max1 - min1)/bins;
  min1 -= 0.5*estep;
  max1 += 0.5*estep;

  for (i=0; i<bnum; i++) {
    
    bv[i] = findb2hg(his,dv,1,betav[i],&beta2v[i],min1,max1,
		     bins,cut_value,eq_weight,eq_height);
    getextrema(his,bins,min1,max1,cut_value,
	       &maxv1,&minv,&maxv2,&maxl1,&minl,&maxl2,&cu1,&cu2);
    val[i] = (exp(maxv1) + exp(maxv2))/(2.0*exp(minv));
    if (verbose) fprintf(stderr,"+"); fflush(stderr);
  }
  free(his);
  return;
}
  
/*********************************************/

void 
FShgram(double **d[],int nd, double beta, double beta2, int bins, double cut,
	int eq_weight, int eq_height) {
  int j,i,k;
  double *his, *op;
  double min1, max1, ee, estep, sleft, sright, a1, a2;
  double minv,maxv1,maxv2,minl,maxl1,maxl2,cu1,cu2;

  if (bins == 0) bins = 100;

  if (verbose) {
    if (eq_weight || eq_height) {
      fprintf(stderr,
	      "- Pseudocritical histogram, starting b = %lg ...\n",
	      p_beta(beta));
    } else
      fprintf(stderr,"- Histogram for b = %lg ...\n",p_beta(beta));
    if (nd > 1) fprintf(stderr,"- Combined hg for %d elements\n",nd);
  }

  his = dblarr(bins);

  get_minmax(d,nd,&min1,&max1);

  if (verbose)
    fprintf(stderr,"  Min: %lg  Max: %lg  Bins: %d\n",min1,max1,bins);

  /* modify min and max a bit ... */

  estep = (max1 - min1)/bins;
  min1 -= 0.5*estep;
  max1 += 0.5*estep;

  beta = findb2hg(his,d,nd,beta,&beta2,min1,max1,bins,cut,eq_weight,eq_height);

  if (verbose && (eq_weight || eq_height)) {
    fprintf(stderr,"critical b = %.12lg \n",p_beta(beta));
    if (b2search)
      fprintf(stderr," - critical b2:  %.12lg \n",beta2);
  }

  if (print_ext) {
    int area;
    getextrema(his,bins,min1,max1,cut,
	       &maxv1,&minv,&maxv2,&maxl1,&minl,&maxl2,&cu1,&cu2);

    fprintf(stderr," - location  --  curvature\n");
    fprintf(stderr," max1  %.10lg    %.10lg\n",maxl1,cu1);
    fprintf(stderr," min   %.10lg\n",minl);
    fprintf(stderr," max2  %.10lg    %.10lg\n",maxl2,cu2);

    fprintf(stderr," - values of max - min - max:\n");
    fprintf(stderr,"   %.10lg    %.10lg    %.10lg\n",maxv1,minv,maxv2);

    ee = (maxl2 - minl)/(maxl2 - maxl1);

    if (h.lx <= h.ly) {
      if (h.lz <= h.ly) area = h.lx * h.lz; else area = h.lx * h.ly;
    } else {
      if (h.lz <= h.lx) area = h.ly * h.lz; else area = h.ly * h.lx;
    }
    fprintf(stderr," - surface tension:  %.10lg\n",
	    (ee*maxv1 + (1-ee)*maxv2 - minv)/(2*area));

  }

  estep = (max1-min1)/bins;
  ee    = min1+0.5*estep;
  for (j=0; j<bins; j++) {
    if ((!hp_low || ee <= hprint_low) && (!hp_high || ee >= hprint_high))
      fprintf(stdout,"%.8lg\t%.12lg\n",ee,his[j]);
    ee += estep;
  }

  free(his);
  return;
}

/**********************************************/

void
FShg(double *his,double **d[],int nd,double beta,double beta2,
     double min1,double max1,int bins) {

  int i,id,k,j,nn;
  double *ap,*a2p,*wp,*dp, v1, dval, logs;

  logs = FSlogsum(beta,beta2);
  map_beta(&beta,&beta2);

  for (i=0; i<bins; i++) his[i] = 0;
  for (i=0; i<n; i++) if (a[i] != NULL) {
    ap = a[i];
    a2p= a2[i];
    wp = wt[i];
    dp = down[i];
    nn = nv[i];

    dval = logs - log(g[i]);

    for (j=0;j<nn;j++) {
      if (!isact2) v1 = exp(-beta*ap[j] - dp[j] - dval);
      else v1 = exp(-beta*ap[j] -beta2*a2p[j] - dp[j] - dval);
      for (id=0; id<nd; id++) {
	k  = floor(((d[id][i][j]-min1)/(max1-min1))*bins);
	if (k == bins) k--;
	his[k] += v1;
      }
    }
  }
}

/**********************************************/

#ifdef AWAY
double
findhg(double hist[],double **d,double beta,double beta2,
       double min1,double max1,int bins,double cut,
       int eq_weight,int eq_height)
{
  int i,k;
  double sleft,sright,a1,a2;
  double minv,maxv1,maxv2,minl,maxl1,maxl2,cu1,cu2;

  if (!eq_weight && !eq_height) {
    FShg(hist,d,beta,beta2,min1,max1,bins);
  } 

  k = 0;
  if (eq_weight || eq_height) {
    double deltab;

    deltab = 1e-7*beta;
    /* now have to solve the critical point */
	
    do {
      k++;

      FShg(hist,d,beta+deltab,beta2,min1,max1,bins);

      sleft = sright = 0;
      for (i=0; i <= (((cut-min1)/(max1-min1))*bins); i++) sleft += hist[i];
      for ( ; i<bins; i++) sright += hist[i];
      a2 = log(sleft/sright);

      FShg(hist,d,beta,beta2,min1,max1,bins);

      sleft = sright = 0;
      for (i=0; i <= (((cut-min1)/(max1-min1))*bins); i++) sleft += hist[i];
      for ( ; i<bins; i++) sright += hist[i];
      a1 = log(sleft/sright);

      beta = beta + (a1*(deltab)/(a1 - a2));
      
    } while ((!eq_height && fabs(a1) > 1e-6 && k < 30) || 
	     (eq_height && k < 10 && fabs(a1) > 1e-2));

  } 
  
  if (eq_height) {
    double delta;

    delta = beta*1e-10;
    k = 0;
    do {
      k++;

      FShg(hist,d,beta+delta,beta2,min1,max1,bins);
      getextrema(hist,bins,min1,max1,cut,
		 &maxv1,&minv,&maxv2,&maxl1,&minl,&maxl2,&cu1,&cu2);
      a2 = maxv1 - maxv2;

      FShg(hist,d,beta,beta2,min1,max1,bins);
      getextrema(hist,bins,min1,max1,cut,
		 &maxv1,&minv,&maxv2,&maxl1,&minl,&maxl2,&cu1,&cu2);
      a1 = maxv1 - maxv2;

      beta = beta + (a1*(delta)/(a1 - a2));
	
    } while (fabs(a1) > 1e-5*fabs(maxv1+maxv2) && k < 30);
  }

  if (k>=30) {
    fprintf(stderr," **** Warning: bc-search did not converge!\n");
    fprintf(stderr,"      norm %g, bc-value %.10g\n",fabs(a1),
	    p_beta(beta));
  }

  return(beta);
}

#endif

double
findhg(double hist[],double **d[],int nd,double beta,double beta2,
       double min1,double max1,int bins,double cut,
       int eq_weight,int eq_height)
{
  int i,k;
  double sleft,sright,a1,a2;
  double minv,maxv1,maxv2,minl,maxl1,maxl2,cu1,cu2;

  if (!eq_weight && !eq_height) {
    FShg(hist,d,nd,beta,beta2,min1,max1,bins);
  } 

  k = 0;
  if (eq_weight || eq_height) {
    double deltab;

    deltab = 1e-7*beta;
    /* now have to solve the critical point */
	
    do {
      k++;

      FShg(hist,d,nd,beta+deltab,beta2,min1,max1,bins);

      sleft = sright = 0;
      for (i=0; i <= (((cut-min1)/(max1-min1))*bins); i++) sleft += hist[i];
      for ( ; i<bins; i++) sright += hist[i];
      a2 = log(sleft/sright);

      FShg(hist,d,nd,beta,beta2,min1,max1,bins);

      sleft = sright = 0;
      for (i=0; i <= (((cut-min1)/(max1-min1))*bins); i++) sleft += hist[i];
      for ( ; i<bins; i++) sright += hist[i];
      a1 = log(sleft/sright);

      beta = beta + (a1*(deltab)/(a1 - a2));
      
    } while ((!eq_height && fabs(a1) > 1e-6 && k < 30) || 
	     (eq_height && k < 10 && fabs(a1) > 1e-2));

  } 
  
  if (eq_height) {
    double delta;

    delta = beta*1e-10;
    k = 0;
    do {
      k++;

      FShg(hist,d,nd,beta+delta,beta2,min1,max1,bins);
      getextrema(hist,bins,min1,max1,cut,
		 &maxv1,&minv,&maxv2,&maxl1,&minl,&maxl2,&cu1,&cu2);
      a2 = maxv1 - maxv2;

      FShg(hist,d,nd,beta,beta2,min1,max1,bins);
      getextrema(hist,bins,min1,max1,cut,
		 &maxv1,&minv,&maxv2,&maxl1,&minl,&maxl2,&cu1,&cu2);
      a1 = maxv1 - maxv2;

      beta = beta + (a1*(delta)/(a1 - a2));
	
    } while (fabs(a1) > 1e-5*fabs(maxv1+maxv2) && k < 30);
  }

  if (k>=30) {
    fprintf(stderr," **** Warning: bc-search did not converge!\n");
    fprintf(stderr,"      norm %g, bc-value %.10g\n",fabs(a1),
	    p_beta(beta));
  }

  return(beta);
}


/**********************************************/

#ifdef AWAY
double
findb2hg(double hist[],double **d,double beta,double *beta2p,
	 double min1,double max1,int bins,double cut,
	 int eq_weight,int eq_height)
{
  double delta,oldbeta2,olda1,prevbeta2,beta2,b1,b2,n1,n2;
  double a1,a2,maxv1,maxv2,minv,maxl1,minl,maxl2,cu1,cu2;
  int k,bracket;

  prevbeta2 = beta2 = *beta2p;
  if (!b2search) {
    beta = findhg(hist,d,beta,beta2,min1,max1,bins,cut,
		  eq_weight,eq_height);
  } else {
    k = 0;
    n1 = n2 = 0;
    delta = beta2 * 1e-6;
    a1 = 1e100;
    do {
      k++;
      olda1 = a1;
      oldbeta2 = prevbeta2;
      prevbeta2 = beta2;

      beta = findhg(hist,d,beta,beta2+delta,min1,max1,bins,cut,
		    eq_weight,eq_height);
      getextrema(hist,bins,min1,max1,cut,
		 &maxv1,&minv,&maxv2,&maxl1,&minl,&maxl2,&cu1,&cu2);
      if (b2search > 0) {
	a2 = (exp(maxv1) + exp(maxv2))/(2.0*exp(minv)) - 2.173;
      } else if (b2search == -2) {
	a2 = FSBindersearch(d,&beta,(beta2+delta));
	a2 -= (1-1.604/3);
      } else {
	betav[0] = beta; beta2v[0] = (beta2+delta); FSinimeas();
	FSBinder(&a2,d,(double *)NULL);
	a2 -= (1-1.604/3);
      }
      beta = findhg(hist,d,beta,beta2,min1,max1,bins,cut,
		    eq_weight,eq_height);
      getextrema(hist,bins,min1,max1,cut,
		 &maxv1,&minv,&maxv2,&maxl1,&minl,&maxl2,&cu1,&cu2);
      if (b2search > 0) {
	a1 = (exp(maxv1) + exp(maxv2))/(2.0*exp(minv)) - 2.173;
      } else if (b2search == -2) {
	a1 = FSBindersearch(d,&beta,beta2);
	a1 -= (1-1.604/3);
      } else {
	betav[0] = beta; beta2v[0] = beta2; FSinimeas();
	FSBinder(&a1,d,(double *)NULL);
	a1 -= (1.0-1.604/3.0);
      }
      beta2 = beta2 + (a1*(delta)/(a1 - a2));	
      
      if (verbose) {
	if (Higgs) fprintf(stderr," M_H: %.10lg  beta: %.10lg  norm %g\n",
			   beta2,beta,a1);
	else fprintf(stderr," beta2: %.10lg  beta: %.10lg  norm %g\n",
		     beta2,beta,a1);
      }
      
      bracket = 0;
      if (a1 > 0 && (n1 == 0.0 || a1 < n1)) { n1 = a1, b1 = prevbeta2; }
      else if (a1 < 0 && (n2 == 0.0 || a1 > n2)) { n2 = a1, b2 = prevbeta2; }
      else bracket = 1;

      if (fabs(a1) > fabs(olda1) || bracket) {
	
	beta2 = (b1 + b2)/2;
	if (a1*olda1 < 0) beta2 = (prevbeta2*olda1 - oldbeta2*a1)/(olda1-a1);
	if (verbose) fprintf(stderr,"Using beta -> %.8lg\n", beta2);
	k += 20;
	
      }
    } while (fabs(a1) > 1e-5 && k < 30);
    if (k>=30) {
      fprintf(stderr," **** Warning: beta2-search did not converge!\n");
      fprintf(stderr,"      norm %g, bc-value %.10g\n",fabs(a1),beta2);
    }
    *beta2p = beta2;
  }
  return(beta);
}

#endif

/**********************************************/

double
findb2hg(double hist[],double **d[],int nd,double beta,double *beta2p,
	 double min1,double max1,int bins,double cut,
	 int eq_weight,int eq_height)
{
  double delta,oldbeta2,olda1,prevbeta2,beta2,b1,b2,n1,n2;
  double a1,a2,maxv1,maxv2,minv,maxl1,minl,maxl2,cu1,cu2;
  int k,bracket;

  prevbeta2 = beta2 = *beta2p;
  if (!b2search) {
    beta = findhg(hist,d,nd,beta,beta2,min1,max1,bins,cut,
		  eq_weight,eq_height);
  } else {
    if (b2search < 0 && nd > 1) 
      fprintf(stderr,
	      "-WARNING: combined hg beta2-search not fully implemented with Binder\n");
    k = 0;
    n1 = n2 = 0;
    delta = beta2 * 1e-6;
    a1 = 1e100;
    do {
      k++;
      olda1 = a1;
      oldbeta2 = prevbeta2;
      prevbeta2 = beta2;

      beta = findhg(hist,d,nd,beta,beta2+delta,min1,max1,bins,cut,
		    eq_weight,eq_height);
      getextrema(hist,bins,min1,max1,cut,
		 &maxv1,&minv,&maxv2,&maxl1,&minl,&maxl2,&cu1,&cu2);
      if (b2search > 0) {
	a2 = (exp(maxv1) + exp(maxv2))/(2.0*exp(minv)) - 2.173;
      } else if (b2search == -2) {
	a2 = FSBindersearch(d[0],&beta,(beta2+delta));
	a2 -= (1-1.604/3);
      } else {
	betav[0] = beta; beta2v[0] = (beta2+delta); FSinimeas();
	FSBinder(&a2,d[0],(double *)NULL);
	a2 -= (1-1.604/3);
      }
      beta = findhg(hist,d,nd,beta,beta2,min1,max1,bins,cut,
		    eq_weight,eq_height);
      getextrema(hist,bins,min1,max1,cut,
		 &maxv1,&minv,&maxv2,&maxl1,&minl,&maxl2,&cu1,&cu2);
      if (b2search > 0) {
	a1 = (exp(maxv1) + exp(maxv2))/(2.0*exp(minv)) - 2.173;
      } else if (b2search == -2) {
	a1 = FSBindersearch(d[0],&beta,beta2);
	a1 -= (1-1.604/3);
      } else {
	betav[0] = beta; beta2v[0] = beta2; FSinimeas();
	FSBinder(&a1,d[0],(double *)NULL);
	a1 -= (1.0-1.604/3.0);
      }
      beta2 = beta2 + (a1*(delta)/(a1 - a2));	
      
      if (verbose) {
	if (Higgs) fprintf(stderr," M_H: %.10lg  beta: %.10lg  norm %g\n",
			   beta2,beta,a1);
	else fprintf(stderr," beta2: %.10lg  beta: %.10lg  norm %g\n",
		     beta2,beta,a1);
      }
      
      bracket = 0;
      if (a1 > 0 && (n1 == 0.0 || a1 < n1)) { n1 = a1, b1 = prevbeta2; }
      else if (a1 < 0 && (n2 == 0.0 || a1 > n2)) { n2 = a1, b2 = prevbeta2; }
      else bracket = 1;

      if (fabs(a1) > fabs(olda1) || bracket) {
	
	beta2 = (b1 + b2)/2;
	if (a1*olda1 < 0) beta2 = (prevbeta2*olda1 - oldbeta2*a1)/(olda1-a1);
	if (verbose) fprintf(stderr,"Using beta -> %.8lg\n", beta2);
	k += 20;
	
      }
    } while (fabs(a1) > 1e-5 && k < 30);
    if (k>=30) {
      fprintf(stderr," **** Warning: beta2-search did not converge!\n");
      fprintf(stderr,"      norm %g, bc-value %.10g\n",fabs(a1),beta2);
    }
    *beta2p = beta2;
  }
  return(beta);
}


/*********************************************/

void 
FShgram2(double **d1,double **d2,double beta, double beta2, 
	 int bins1, int bins2) {

  int j,i,k;
  double *his, *op;
  double min1, max1, min2, max2, ee1, es1, ee2, es2;

  if (bins1 == 0) bins1 = 20;
  if (bins2 == 0) bins2 = 20;

  if (verbose)
    fprintf(stderr,"- Double histogram for b = %lg ...\n",p_beta(beta));

  his = dblarr(bins1*bins2);

  min1 = 1e100;
  max1 = -min1;
  for (i=0; i<n; i++) {
    op = d1[i];
    for (j=0; j<nv[i]; j++) {
      min1 = (min1 < op[j]) ? min1 : op[j];
      max1 = (max1 > op[j]) ? max1 : op[j];
    }
  }

  min2 = 1e100;
  max2 = -min2;
  for (i=0; i<n; i++) {
    op = d2[i];
    for (j=0; j<nv[i]; j++) {
      min2 = (min2 < op[j]) ? min2 : op[j];
      max2 = (max2 > op[j]) ? max2 : op[j];
    }
  }

  if (verbose) {
    fprintf(stderr,"  Min1: %lg  Max1: %lg  Bins: %d\n",min1,max1,bins1);
    fprintf(stderr,"  Min2: %lg  Max2: %lg  Bins: %d\n",min2,max2,bins2);
  }

  FShg2(his,d1,d2,beta,beta2,min1,max1,bins1,min2,max2,bins2);

  es1 = (max1-min1)/bins1;
  es2 = (max2-min2)/bins2;

  ee2 = min2+0.5*es2;
  for (k=0; k<bins2; k++) {   
    ee1 = min1+0.5*es1;
    for (j=0; j<bins1; j++) {
      fprintf(stdout,"%.8lg\t%.8lg\t%.12lg\n",ee1,ee2,his[j+bins1*k]);
      ee1 += es1;
    }
    ee2 += es2;
    fprintf(stdout,"\n");
  }

  free(his);
  return;
}

/**********************************************/

void
FShg2(double *hg,double **d1,double **d2,double beta,double beta2,
      double min1,double max1,int bins1,
      double min2,double max2,int bins2) {

  int i,k1,k2,j,nn;
  double *ap,*a2p,*wp,*dp,*op1,*op2, v1, dval, logs;

  logs = FSlogsum(beta,beta2);
  map_beta(&beta,&beta2);

  for (i=0; i<bins1*bins2; i++) hg[i] = 0;
  for (i=0; i<n; i++) if (a[i] != NULL) {
    ap = a[i];
    a2p= a2[i];
    wp = wt[i];
    dp = down[i];
    op1 = d1[i]; op2 = d2[i];
    nn = nv[i];

    dval = logs - log(g[i]);

    for (j=0;j<nn;j++) {
      if (!isact2) v1 = exp(-beta*ap[j] - dp[j] - dval);
      else v1 = exp(-beta*ap[j] -beta2*a2p[j] - dp[j] - dval);
      k1 = floor(((op1[j]-min1)/(max1-min1))*bins1);
      if (k1 == bins1) k1--;
      k2 = floor(((op2[j]-min2)/(max2-min2))*bins2);
      if (k2 == bins2) k2--;
      hg[k1 + k2*bins1] += v1;
    }
  }
}



/*****************************/

double 
jack_c(double d[],int jack, double *ee) {

  int k;
  double v,e;
  
  v = e = 0;
  for (k=0; k<jack; k++) v += d[k]/jack;
  for (k=0; k<jack; k++) e += sqr(v - d[k]);
  *ee = sqrt((jack-1)*e/jack);
  return(v);
}


#ifdef AWAY
void 
FShgJack(double **dat,double beta, double beta2, int bins, double cut,
	 int jack,int jprint,char *js,char *ss,
	 int eq_weight, int eq_height) {

  int i,j,k,ntop,nbot,iter;
  double min1, max1, ee, estep, sleft, sright, a1, a2;
  FILE *out;
  double *op, *jp, *dp, *his, *hist, *sig, be[200], be2[200];
  double *jdat[MAXRUN];
  double mv1[MAXJACK],mv2[MAXJACK],miv[MAXJACK];
  double ml1[MAXJACK],ml2[MAXJACK],mil[MAXJACK];
  double cv1[MAXJACK],cv2[MAXJACK];
  double low_s[MAXJACK],high_s[MAXJACK];
  double minv,maxv1,maxv2,minl,maxl1,maxl2,cu1,cu2;

  if (bins == 0) bins = 100;

  if (verbose) if (eq_weight || eq_height)
    fprintf(stderr,
	    "- Pseudocritical histogram, starting b = %lg ...\n",
	    p_beta(beta));
  else
    fprintf(stderr,"- Histogram for b = %lg ...\n",p_beta(beta));

  his = dblarr(bins);
  hist = dblarr(bins);
  sig = dblarr(bins);
  
  min1 = 1e100;
  max1 = -min1;
  for (i=0; i<n; i++) {
    op = dat[i];
    for (j=0; j<nv[i]; j++) {
      min1 = (min1 < op[j]) ? min1 : op[j];
      max1 = (max1 > op[j]) ? max1 : op[j];
    }
  }

  if (verbose)
    fprintf(stderr,"  Min: %lg  Max: %lg  Bins: %d\n",min1,max1,bins);

  /* modify min and max ... */

  estep = (max1 - min1)/bins;
  min1 -= 0.5*estep;
  max1 += 0.5*estep;


  for (i=0; i<n; i++) jdat[i] = dblarr(nn[i] - nn[i]/jack + 1);

  if (jprint) if (jstdout) out = stdout; else out = fopen(js,"w");

  for (iter=0; iter<jack; iter++) {
    FSiniJack(iter,jack,ss);

    for (i=0; i<n; i++) {
      nbot  = iter*nn[i]/jack;
      ntop  = (iter+1)*nn[i]/jack;

      jp = jdat[i];
      dp = dat[i];

#pragma ivdep
      for (j=0; j<nbot; j++) jp[j] = dp[j];
#pragma ivdep
      for (j=ntop; j<nn[i]; j++) jp[j-ntop+nbot] = dp[j];
    }

    beta = findb2hg(hist,jdat,beta,&beta2,min1,max1,bins,cut,
		    eq_weight,eq_height);

    be[iter]  = beta;
    be2[iter] = beta2;
    
    for (i=0; i<bins; i++) {
      his[i] += hist[i];
      sig[i] += sqr(hist[i]);
    }

    if (print_ext) {
      double pt;

      getextrema(hist,bins,min1,max1,cut,
		 &mv1[iter],&miv[iter],&mv2[iter],
		 &ml1[iter],&mil[iter],&ml2[iter],
		 &cv1[iter],&cv2[iter]);

      /* get the exp. values too --- */
      sleft = sright = 0;
      estep = (max1-min1)/bins;
      ee    = min1+0.5*estep;
      pt    = 0;
      for (i=0; i <= (((cut-min1)/(max1-min1))*bins); i++) {
	sleft += hist[i]*ee;
	pt += hist[i];
	ee += estep;
      }
      low_s[iter]  = sleft/pt;
      pt = 0;
      for ( ; i<bins; i++) {
	sright += hist[i]*ee;
	pt += hist[i];
	ee += estep;
      }
      high_s[iter] = sright/pt;
    }

    /* write stuff */

    if (jprint) {
      estep = (max1-min1)/bins;
      ee    = min1+0.5*estep;
      for (j=0; j<bins; j++) {
	if ((!hp_low || ee <= hprint_low) && (!hp_high || ee >= hprint_high))
	  /* fprintf(out,"%d %d  %.10lg  %.12lg\n",iter,j,ee,hist[j]); */
	  fprintf(out,"%.10lg  %.12lg\n",ee,hist[j]);
	ee += estep;
      }
    }
  }
  if (jprint && !jstdout) fclose(out);

#pragma ivdep
  for (k=0; k<bins; k++) his[k] /= jack;

#pragma ivdep
  for (k=0; k<bins; k++)
    sig[k] = sqrt((jack-1)*(sig[k]/jack - sqr(his[k])));

  estep = (max1-min1)/bins;
  ee    = min1+0.5*estep;
  if (!jstdout) for (j=0; j<bins; j++) {
    if ((!hp_low || ee <= hprint_low) && (!hp_high || ee >= hprint_high))
      fprintf(stdout,"%.8lg\t%.12lg\t%.12lg\n",ee,his[j],sig[j]);
    ee += estep;
  }

  for (iter=0; iter<jack; iter++) be[iter] = p_beta(be[iter]);
  beta = jack_c(be,jack,&ee);
  fprintf(stderr," - critical b:   %.12lg   %.12lg\n",beta,ee);
  if (b2search) {
    beta2 = jack_c(be2,jack,&ee);
    fprintf(stderr," - critical b2:  %.12lg   %.12lg\n",beta2,ee);
    fprintf(stderr," - block by block: \n");
    for (iter=0; iter<jack; iter++) fprintf(stderr,"%d %.12g  %.12g\n",
					    iter,be[iter],be2[iter]);
  }

  if (print_ext) {
    double val,ee,tmp[MAXJACK],area;

    fprintf(stderr," -- extrema locations --\n");
    val = jack_c(ml1,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = maxloc1\n",val,ee);
    val = jack_c(mil,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = minloc\n",val,ee);
    val = jack_c(ml2,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = maxloc2\n",val,ee);

    for (k=0; k<jack; k++) tmp[k] = ml2[k] - ml1[k];
    val = jack_c(tmp,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = loc2-loc1\n",val,ee);

    fprintf(stderr," -- extrema values -- \n");
    val = jack_c(mv1,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = maxval1\n",val,ee);
    val = jack_c(miv,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = minval\n",val,ee);
    val = jack_c(mv2,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = maxval2\n",val,ee);

    fprintf(stderr," -- curvatures -- \n");
    val = jack_c(cv1,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = curvature1\n",val,ee);
    val = jack_c(cv2,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = curvature2\n",val,ee);
    
    if (h.lx <= h.ly) {
      if (h.lz <= h.ly) area = h.lx * h.lz; else area = h.lx * h.ly;
    } else {
      if (h.lz <= h.lx) area = h.ly * h.lz; else area = h.ly * h.lx;
    }
    for (k=0; k<jack; k++) {
      ee = (ml2[k] - mil[k])/(ml2[k] - ml1[k]);
      tmp[k] = (ee*mv1[k] + (1-ee)*mv2[k] - miv[k])/(2*area);
    }
    val = jack_c(tmp,jack,&ee);
    fprintf(stderr," -- surface tension value --\n");
    fprintf(stderr," %.10lg  %.10lg    = tension\n",val,ee);

    fprintf(stderr," -- low and high expectation values --\n");
    val = jack_c(low_s,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = low\n",val,ee);
    val = jack_c(high_s,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = high\n",val,ee);
    for (k=0; k<jack; k++) tmp[k] = high_s[k] - low_s[k];
    val = jack_c(tmp,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = h-l\n",val,ee);
  }

}

#endif

void 
FShgJack(double **dat[],int nd,double beta, double beta2, int bins, double cut,
	 int jack,int jprint,char *js,char *ss,
	 int eq_weight, int eq_height) {

  int i,j,k,ntop,nbot,iter;
  double min1, max1, ee, estep, sleft, sright, a1, a2;
  FILE *out;
  double *op, *jp, *dp, *his, *hist, *sig, be[200], be2[200];
  double ***jdat;
  double mv1[MAXJACK],mv2[MAXJACK],miv[MAXJACK];
  double ml1[MAXJACK],ml2[MAXJACK],mil[MAXJACK];
  double cv1[MAXJACK],cv2[MAXJACK];
  double low_s[MAXJACK],high_s[MAXJACK];
  double minv,maxv1,maxv2,minl,maxl1,maxl2,cu1,cu2;


  if (bins == 0) bins = 100;

  if (verbose) {
    if (eq_weight || eq_height)
      fprintf(stderr,
	      "- Pseudocritical histogram, starting b = %lg ...\n",
	      p_beta(beta));
    else
      fprintf(stderr,"- Histogram for b = %lg ...\n",p_beta(beta));
    if (nd > 1) fprintf(stderr,"- Combined hg for %d data\n",nd);
  }

  his = dblarr(bins);
  hist = dblarr(bins);
  sig = dblarr(bins);
  
  get_minmax(dat,nd,&min1,&max1);

  if (verbose)
    fprintf(stderr,"  Min: %lg  Max: %lg  Bins: %d\n",min1,max1,bins);

  /* modify min and max ... */

  estep = (max1 - min1)/bins;
  min1 -= 0.5*estep;
  max1 += 0.5*estep;

  jdat = (double ***)calloc(nd,sizeof(double **));
  for (j=0; j<nd; j++) jdat[j] = (double **)calloc(n,sizeof(double *));
  for (j=0; j<nd; j++) for (i=0; i<n; i++) 
    jdat[j][i] = dblarr(nn[i] - nn[i]/jack + 1);

  if (jprint) if (jstdout) out = stdout; else out = fopen(js,"w");

  for (iter=0; iter<jack; iter++) {
    FSiniJack(iter,jack,ss);

    for (i=0; i<n; i++) {
      nbot  = iter*nn[i]/jack;
      ntop  = (iter+1)*nn[i]/jack;

      for (k=0; k<nd; k++) {
	jp = jdat[k][i];
	dp = dat[k][i];
#pragma ivdep
	for (j=0; j<nbot; j++) jp[j] = dp[j];
#pragma ivdep
	for (j=ntop; j<nn[i]; j++) jp[j-ntop+nbot] = dp[j];
      }
    }

    beta = findb2hg(hist,jdat,nd,beta,&beta2,min1,max1,bins,cut,
		    eq_weight,eq_height);

    be[iter]  = beta;
    be2[iter] = beta2;
    
    for (i=0; i<bins; i++) {
      his[i] += hist[i];
      sig[i] += sqr(hist[i]);
    }

    if (print_ext) {
      double pt;

      getextrema(hist,bins,min1,max1,cut,
		 &mv1[iter],&miv[iter],&mv2[iter],
		 &ml1[iter],&mil[iter],&ml2[iter],
		 &cv1[iter],&cv2[iter]);

      /* get the exp. values too --- */
      sleft = sright = 0;
      estep = (max1-min1)/bins;
      ee    = min1+0.5*estep;
      pt    = 0;
      for (i=0; i <= (((cut-min1)/(max1-min1))*bins); i++) {
	sleft += hist[i]*ee;
	pt += hist[i];
	ee += estep;
      }
      low_s[iter]  = sleft/pt;
      pt = 0;
      for ( ; i<bins; i++) {
	sright += hist[i]*ee;
	pt += hist[i];
	ee += estep;
      }
      high_s[iter] = sright/pt;
    }

    /* write stuff */

    if (jprint) {
      estep = (max1-min1)/bins;
      ee    = min1+0.5*estep;
      for (j=0; j<bins; j++) {
	if ((!hp_low || ee <= hprint_low) && (!hp_high || ee >= hprint_high))
	  /* fprintf(out,"%d %d  %.10lg  %.12lg\n",iter,j,ee,hist[j]); */
	  fprintf(out,"%.10lg  %.12lg\n",ee,hist[j]);
	ee += estep;
      }
    }
  }
  if (jprint && !jstdout) fclose(out);

#pragma ivdep
  for (k=0; k<bins; k++) his[k] /= jack;

#pragma ivdep
  for (k=0; k<bins; k++)
    sig[k] = sqrt((jack-1)*(sig[k]/jack - sqr(his[k])));

  estep = (max1-min1)/bins;
  ee    = min1+0.5*estep;
  if (!jstdout) for (j=0; j<bins; j++) {
    if ((!hp_low || ee <= hprint_low) && (!hp_high || ee >= hprint_high))
      fprintf(stdout,"%.8lg\t%.12lg\t%.12lg\n",ee,his[j],sig[j]);
    ee += estep;
  }

  for (iter=0; iter<jack; iter++) be[iter] = p_beta(be[iter]);
  beta = jack_c(be,jack,&ee);
  fprintf(stderr," - critical b:   %.12lg   %.12lg\n",beta,ee);
  if (b2search) {
    beta2 = jack_c(be2,jack,&ee);
    fprintf(stderr," - critical b2:  %.12lg   %.12lg\n",beta2,ee);
    fprintf(stderr," - block by block: \n");
    for (iter=0; iter<jack; iter++) fprintf(stderr,"%d %.12g  %.12g\n",
					    iter,be[iter],be2[iter]);
  }

  if (print_ext) {
    double area,val,ee,tmp[MAXJACK];

    fprintf(stderr," -- extrema locations --\n");
    val = jack_c(ml1,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = maxloc1\n",val,ee);
    val = jack_c(mil,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = minloc\n",val,ee);
    val = jack_c(ml2,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = maxloc2\n",val,ee);

    for (k=0; k<jack; k++) tmp[k] = ml2[k] - ml1[k];
    val = jack_c(tmp,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = loc2-loc1\n",val,ee);

    fprintf(stderr," -- extrema values -- \n");
    val = jack_c(mv1,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = maxval1\n",val,ee);
    val = jack_c(miv,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = minval\n",val,ee);
    val = jack_c(mv2,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = maxval2\n",val,ee);

    fprintf(stderr," -- curvatures -- \n");
    val = jack_c(cv1,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = curvature1\n",val,ee);
    val = jack_c(cv2,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = curvature2\n",val,ee);
    
    if (h.lx <= h.ly) {
      if (h.lz <= h.ly) area = h.lx * h.lz; else area = h.lx * h.ly;
    } else {
      if (h.lz <= h.lx) area = h.ly * h.lz; else area = h.ly * h.lx;
    }
    for (k=0; k<jack; k++) {
      ee = (ml2[k] - mil[k])/(ml2[k] - ml1[k]);
      tmp[k] = (ee*mv1[k] + (1-ee)*mv2[k] - miv[k])/(2*area);
    }
    val = jack_c(tmp,jack,&ee);
    fprintf(stderr," -- surface tension value --\n");
    fprintf(stderr," %.10lg  %.10lg    = tension\n",val,ee);

    fprintf(stderr," -- low and high expectation values --\n");
    val = jack_c(low_s,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = low\n",val,ee);
    val = jack_c(high_s,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = high\n",val,ee);
    for (k=0; k<jack; k++) tmp[k] = high_s[k] - low_s[k];
    val = jack_c(tmp,jack,&ee);
    fprintf(stderr," %.10lg  %.10lg    = h-l\n",val,ee);
  }

  for (j=0; j<nd; j++) {
    for (i=0; i<n; i++) free(jdat[j][i]);
    free(jdat[j]);
  }
}

/***************************************************************/

void get_minmax(double **d[],int nd, double *min1, double *max1)
{
  double mi,ma,*op;
  int i,k,j;

  mi = 1e100;
  ma = -mi;
  for (k=0; k<nd; k++) for (i=0; i<n; i++) {
    op = d[k][i];
    for (j=0; j<nv[i]; j++) {
      mi = (mi < op[j]) ? mi : op[j];
      ma = (ma > op[j]) ? ma : op[j];
    }
  }
  *min1 = mi; *max1 = ma;
}

/***************************************************************/

void fitconst(double h[],int bins, int ind, int len, double *val);
void fitextrema(double h[],int bins,double min1,double max1,int ind,int len,
		double *val,double *loc,double *curv);

void 
getextrema(double hist[],int bins,double min1,double max1, double cut,
	   double *maxv1,double *minv,double *maxv2,
	   double *maxl1,double *minl,double *maxl2,
	   double *cu1,double *cu2) {

  int i,mxl1,mxl2,mnl;
  double mx1,mx2,mn;

  mx1 = mx2 = 0;
  for (i=0; i <= (((cut-min1)/(max1-min1))*bins); i++) 
    if (hist[i] > mx1) mx1 = hist[mxl1 = i];
  for ( ; i<bins; i++) if (hist[i] > mx2) mx2 = hist[mxl2 = i];
  mn = mx1;
  for (i=mxl1; i<=mxl2; i++) if (hist[i] < mn) mn = hist[mnl = i];

  /*  fprintf(stderr,"Maxl: %d -- minl %d -- maxl %d\n",mxl1,mnl,mxl2); */

  /* now maximum 1 - minimum - maximum 2 */
  fitextrema(hist,bins,min1,max1,mxl1,max_r1,maxv1,maxl1,cu1);

  if (!flat) fitextrema(hist,bins,min1,max1,mnl, min_r, minv, minl, cu2);
  else {
    fitconst(hist,bins,(int)(((cut-min1)/(max1-min1))*bins),min_r,minv);
    *minl = cut;
  }

  fitextrema(hist,bins,min1,max1,mxl2,max_r2,maxv2,maxl2,cu2);
  
}

void 
fitconst(double h[],int bins, int ind, int len, double *val) {

  int i,j;
  double pe[1],x[300],y[300];

  if (ind-len/2 < 0 || ind+len/2 >= bins || len > 300) {
    halt(" *** point %d too close to end\n",ind);
  }
  for (j=i=0; i<len; i++) if (h[ind+i-len/2] > 0) {
    x[j] = i;
    y[j] = log(h[ind + i - len/2]);
    j++;
  }

  polyfit(j,x,y,NULL,0,val,pe);

}

void 
fitextrema(double h[],int bins,double min1,double max1,int ind,int len,
		double *val,double *loc,double *curv) {

  int i,j;
  double x[300],y[300],e[300],sig[300],p[3],pe[3],step;

  if (ind-len/2 < 0 || ind+len/2 >= bins) 
    halt(" *** point %d too close to end\n",ind);
    
  j=0;
  step = (max1-min1)/bins;
  for (i=ind-len/2; i<=ind+len/2; i++) if (h[i] > 0) {
    x[j] = min1 + 0.5*step + i*step;
    y[j] = log(h[i]);
    e[j] = exp(sqr((2.0*(i-ind))/len)*2.0);
    j++;
  }

  if (j > 3) polyfit(j,x,y,e,2,p,pe);
  else fprintf(stderr,"too small range\n");
  *loc  = -p[1]/(2*p[2]);
  *val  = p[0] + p[1]*(*loc) + p[2]*sqr(*loc);
  *curv = 2*p[2];
  if (*loc < x[0] || *loc > x[j-1])
    fprintf(stderr,"*Warning: extremum range %g - %g, value %g\n",
	    x[0],x[j-1],*loc);
  
}


/****************************************************************
 *                                                              *
 *   the NR-minimizing function                                 *
 *   Kari Rummukainen 1991                                      *
 *                                                              *
 *                                                              *
 *   global variable bsum contains derivative hash table        *
 *                                                              *
 ***************************************************************/

void NRFuncInit();
void NRFunc(double val[],double vec[]);
void NRDeriv(double matrix[]);

static double fnr[MAXRUN];
static int n_ind;
static int *bsum;

void 
NRFuncInit() {

  int i,j,k,l,m,new;
  double bet,bet2;

  n_ind = n_independent;

  for (j=0; j<n_ind; j++) fnr[j] = f[ pb_i[j][0] ];

  for (j=0; j<n_ind; j++) for (i=0; i<n_ind; i++) {
    bet = b_ind[i]+b_ind[j];
    bet2= b2_ind[i]+b2_ind[j];
    new = 1;
    for (k=0; k<n_ind*j+i && new; k++) {
      m = k % n_ind;
      l = (k - m)/n_ind;
      if (bet == (b_ind[l]+b_ind[m]) && bet2 == (b2_ind[l]+b2_ind[m])) {
	bsum[n_ind*j+i] = k;
	new = 0;
      }
    }
    if (new) bsum[n_ind*j+i] = -1;
  }
}


/********************************/

#define NLIM 1.0

void 
NRFunc(double val[],double v[]) {

  int i,j,k;
  double norm;

  norm = 0.0;
  for (k=1; k<n_ind; k++) {
    norm += sqr(-logz[k] - fnr[k] + logz[0]);
    fnr[k] = v[k]- v[0];
    for (j=0; j<nb_i[k]; j++) f[ pb_i[k][j] ] = fnr[k];
  }
  norm = sqrt(norm);

  fnr[0] = 0.0;
  for (i=0; i<nb_i[0]; i++) f[i] = 0.;
  if (norm < NLIM*n_ind) FSinidown(); else FSsafeinit();

  for (i=0; i<n_ind; i++) val[i] = -logz[i] - fnr[i] + logz[0];
}

/*******************************/

#define NOINVERT 1.0e-100

void 
NRDeriv(double * dmatr) {

  int i,j,k,m,l,ij;
  double val,bet,bet2,mp,mpt;
  double *ap,*a2p,*wp,*dp,dm,dd;

  for (j=0; j<n_ind; j++) for (i=0; i<n_ind; i++) {

    /* check if have to calculate */

    k = bsum[n_ind*j+i];
    m = k % n_ind;
    l = (k-m)/n_ind;

    if (k < 0 || (dmatr[n_ind*l+m] < NOINVERT)) {

      /* now let us calculate */

      bet = b_ind[i] + b_ind[j];
      bet2= b2_ind[i] + b2_ind[j];
      val = 0.0;
      mpt = 0.0;
      for (ij=0; ij < nb_i[j]; ij++)
	mpt += nv[ pb_i[j][ij] ]*g[ pb_i[j][ij] ];

      for (k=0; k<n; k++) {
	ap = a[k];
	a2p= a2[k];
	wp = wt[k];
	dp = down[k];
	mp = log(g[k]*mpt) + fnr[j] - logz[i];

	if (!isact2)
#pragma ivdep
	  for (m=0; m<nv[k]; m++) val += exp(-bet*ap[m] - 2.0*dp[m] + mp);
	else
#pragma ivdep
	  for (m=0; m<nv[k]; m++) 
	    val += exp(-bet*ap[m] - bet2*a2p[m] - 2.0*dp[m] + mp);
      }
      dmatr[n_ind*j+i] = val;

    }

    else {

      /* now can use old one */

      if (dmatr[n_ind*l+m] == 0.0) dmatr[n_ind*j+i] = 0.0;
      else {
	dm = dd = 0.;
	for (ij=0; ij<nb_i[j]; ij++) dm += nv[ pb_i[j][ij] ] * g[ pb_i[j][ij] ];
	for (ij=0; ij<nb_i[l]; ij++) dd += nv[ pb_i[l][ij] ] * g[ pb_i[l][ij] ];
	dmatr[n_ind*j+i] = dmatr[n_ind*l+m]
	  * dm/dd * exp(fnr[j] - fnr[l] - logz[i] + logz[m]);
      }
    }
  }

  for (j=0; j<n_ind; j++) {
    for (i=1; i<n_ind; i++) dmatr[n_ind*j+i] -= dmatr[n_ind*j];
    dmatr[n_ind*j] = dmatr[j] = 0.0;
    dmatr[n_ind*j+j] -= 1.0;
  }

}


/****************************************************************************
*
*   The Newton-Raphson multidimensional solving system
*
*     Kari Rummukainen 20 Aug 1990
*
*   parameters:
*
*   int NewtonRaphson(int n,double vec[],int maxl,double eps)
*        nv    - number of variables
*        vec   - solution vector + initializing vector!
*        maxl  - maximum number of iterations
*        eps   - required accuracy
*
*   Returns the number of iterations
*/

#define TOLER 1.0e-8
#define MVEC 50
#define DELTA 1.0e-9


int 
NewtonRaphson(int n,double vec[],int maxl,double eps) 
{
  double *v,*val,*val2,*vv;
  double *matv,det,tol;
  double norm;
  double mult,onorm;
  int   i,j,loop;

  loop    = 0;
  tol     = TOLER;
  matv    = dblarr(n*n);
  bsum    = (int *)calloc(n*n,sizeof(int));
  v       = dblarr(n);
  val     = dblarr(n);
  val2    = dblarr(n);
  vv      = dblarr(n);

  /* initialize functions */

  for (j=0; j<n; j++) v[j] = vec[j];
  NRFuncInit();

  norm = 0.0;

  NRFunc(val,v);

  for (j=0; j<n; j++) {
    norm += sqr(val[j]);
  }
  norm = sqrt(norm);

  mult = 1;
  while (norm > eps) {
    loop++;

    /* calculate the derivative matrix (dF_i/df_j) */

    NRDeriv(matv);
    for (i=0; i<n; i++) vv[i] = -val[i];

/*
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) fprintf(stderr,"%8.4lg ",matv[n*i+j]);
      fprintf(stderr,"   %8.4lg\n",vv[i]);
    }
*/
    svdecomp(matv,n,vv,1);

/*    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) fprintf(stderr,"%8.4lg ",matv[n*i+j]);
      fprintf(stderr,"   %8.4lg\n",vv[i]);
    }
*/

    for (i=0; i<n; i++) val2[i] = v[i];

    onorm = norm;
    j = 0;
    mult = 1;

    do {
      j++;
      norm = 0.0;

      for (i=0; i<n; i++) v[i] = val2[i] + mult*vv[i];

      NRFunc(val,v);

      for (i=0; i<n; i++) norm += sqr(val[i]);
      norm = sqrt(norm);
      mult *= -0.1;
    } while (onorm < norm && j < 20);

    if (verbose) fprintf(stderr," - NR - %d \tnorm: %lg fork: %d\n",loop,norm,j);
    if (onorm < norm) {
      /* fprintf(stderr," - NR - Newton-Raphson failure, sorry => norm <- 0\n"); */
      fprintf(stderr," - NR - Newton-Raphson failure, sorry => naive iteration\n");

      j = 0;
      do {
	j++;
	for (i=0; i<n; i++) v[i] = 0.9*v[i] + 0.1*val[i];
	NRFunc(val,v);
	norm = 0.0;
	for (i=0; i<n; i++) norm += sqr(val[i]);
	norm = sqrt(norm);
	fprintf(stderr," - Naive - %d norm %g\n",j,norm);
      } while (norm > onorm && j < 20);
    }
  }

  for (j=0; j<n; j++) vec[j] = v[j];

  free(v); free(vv); free(val); free(val2);
  free(matv); free(bsum);
  return(loop);
}


/*************************/

void 
FSCpart(double *imbeta,int ni)
{
  int i,j,k,l,ntop,nbot;
  complex *val;
  double bi,b2;

  val = (complex *)calloc(bnum,sizeof(complex));
  if (ni > 1 && ni < bnum) 
    halt(" ** wrong number of imaginary parts\n",NULL);

  for (i=0; i<bnum; i++) {
    if (ni > 1) bi = imbeta[i]; else bi = imbeta[0];
    if (isact2) b2 = beta2v[i];
    val[i] = FSpartition(betav[i],b2,bi);
  }
  for (i=0; i<bnum; i++) {
    fprintf(stdout,"%.12lg  ",p_beta(betav[i]));
    if (isact2 && !is_su3h) fprintf(stdout,"%.12lg  ",(beta2v[i]));
    fprintf(stdout,"%.14lg  %.14lg\n",val[i].r, val[i].i);
  }
  free(val);
}


/*************************/

void 
FShgvec(double **dv[],int bins,double cut_value,int eq_weight,int eq_height)
{
  int i,j,k,l,ntop,nbot;
  double *val,*bv;

  val = dblarr(bnum);
  bv  = dblarr(bnum);
  
  FShgvect(bv,val,dv,bins,cut_value,eq_weight,eq_height);

  for (i=0; i<bnum; i++) {
    fprintf(stdout,"%.12lg  ",p_beta(bv[i]));
    if (isact2 && !is_su3h) fprintf(stdout,"%.12lg  ",(beta2v[i]));
    fprintf(stdout,"%.14lg\n",val[i]);
  }

  free(val);
  free(bv);
}

/**********************************************************
 *  calculate the second moment eigenmatrix
 */

void
FSeig(double **dv[],int nd,int edprint,char *lists[])
{
  int i,j,k;
  double beta, *jp, *dp;
  double *resm,*eval,*evec;

  if (bnum != 1) halt("One and only one beta-value!",NULL);

  /* allocate workspace, if needed */
  
  eval = dblarr(nd);
  evec = dblarr(nd*nd);
  resm = dblarr(nd*nd);

  FSm2mat(dv,nd,resm);
  jacobi(nd,resm,eval,evec,1);

  if (!edprint) for (i=0; i<nd; i++) {
    printf("%d eigenvalue %.14lg \n",i,eval[i]);
    printf("%d vector  ",i);
    for (j=0; j<nd; j++) printf("%.8lg ",evec[i*nd+j]);
    printf("\n");
  } else {
    for (i=0; i<nd; i++) {
      /* printf("# %d eigenvalue %.14lg \n",i,eval[i]); */
      printf("set vect%d = \'",i);
      for (j=0; j<nd; j++) {
	printf("%.8lg*(%s)",evec[i*nd+j],lists[j+nd*i]);
	if (j<nd-1) printf(" + ");
	else printf("\'\n");
      }
    }
  }
}

/**************/

void 
FSiniJack(int iter, int jack, char *pre) {

  int i,j,k,ntop,nbot;
  static FILE * fil;
  static double ** ac, ** a2c, ** wc;
  char name[500];
  double beta, *ap, *a2p, *pp, *p2p, *wp, *w0p;

  if (verbose) fprintf(stderr," - %2d: ",iter+1); fflush(stdout);

  if (iter == 0) {
    strcpy(name,pre);  strcat(name,".jack");
    if ((fil = fopen(name,"r")) == NULL) 
      halt(" ++++ error in jackknife file %s",name);

    for (i=0; i<n; i++) {
      f[i]  = fc[i];
      nv[i] = nn[i];
    }
    ac = a;   a2c = a2; 
    a  = akt; a2  = akt2;
    wc = wt;
    wt = wkt;
  }

  for (i=0; i<n; i++) {
    fscanf(fil,"%d %lg %lg",
	   &k,&b[i],&f[i]);
    nbot  = iter*nn[i]/jack;
    ntop  = (iter+1)*nn[i]/jack;
    nv[i] = nn[i] - ntop + nbot;
    ap = ac[i];
    a2p= a2c[i];
    pp = akt[i];
    p2p= akt2[i];
    wp = wkt[i];
    w0p= wc[i];

#pragma ivdep
    for (j=0; j<nbot; j++) pp[j] = ap[j];
#pragma ivdep
    for (j=ntop; j<nn[i]; j++) pp[j-ntop+nbot] = ap[j];

    if (isact2) {
#pragma ivdep
      for (j=0; j<nbot; j++) p2p[j] = a2p[j];
#pragma ivdep
      for (j=ntop; j<nn[i]; j++) p2p[j-ntop+nbot] = a2p[j];
    }

#pragma ivdep
    for (j=0; j<nbot; j++) wp[j] = w0p[j];
#pragma ivdep
    for (j=ntop; j<nn[i]; j++) wp[j-ntop+nbot] = w0p[j];
  }

  FSsafeinit();
  FSinimeas();

  if (iter+1 == jack) fclose(fil);
}


/*****************************/

void
FSJack(double **dv[],int nd,int jack,int isfree,int jprint,char *js,
       char *ss,double mom,int binder,int cross)
{
  int i,j,k,nv,ntop,nbot,iter;
  FILE *out;
  double *jp, *dp;
  double **res, **sig, **resl,***jdat;
  double **average;

  /* in moments, calculate the average globally, without jackknife */

  average = (double **)calloc(nd,sizeof(double *));
  if (is_global_ave && (mom || cross || binder)) {
    fprintf(stderr," - calculating global averages .."); fflush(stdout);
    for (j=0; j<nd; j++) {
      average[j] = dblarr(bnum);
      FSval(average[j],dv[j]);
    }
    fprintf(stderr,". done\n");
  } else {
    for (j=0; j<nd; j++) average[j] = (double *)NULL;
  }

  /* allocate workspace, if needed */

  if (isfree) nd = nv = 1;
  else if (cross) nv = 1;
  else nv = nd;
  if (!isfree) {
    jdat = (double ***)calloc(nd,sizeof(double **));
    for (j=0; j<nd; j++) {
      jdat[j] = (double **)calloc(n,sizeof(double *));
      for (i=0; i<n; i++) jdat[j][i] = dblarr(nn[i] - nn[i]/jack + 1);
    }
  }
  res = (double **)calloc(nv,sizeof(double *));
  sig = (double **)calloc(nv,sizeof(double *));
  resl = (double **)calloc(nv,sizeof(double *));

  for (j=0; j<nv; j++) {
    res[j] = dblarr(bnum);
    sig[j] = dblarr(bnum);
    resl[j] = dblarr(bnum);
  }

  if (jprint) if (jstdout) out = stdout; else out = fopen(js,"w");

  for (iter=0; iter<jack; iter++) {
    FSiniJack(iter,jack,ss);

    if (!isfree) for (i=0; i<n; i++) {
      nbot  = iter*nn[i]/jack;
      ntop  = (iter+1)*nn[i]/jack;

      for (k=0; k<nd; k++) {
	jp = jdat[k][i];
	dp = dv[k][i];
#pragma ivdep
	for (j=0; j<nbot; j++) jp[j] = dp[j];
#pragma ivdep
	for (j=ntop; j<nn[i]; j++) jp[j-ntop+nbot] = dp[j];
      }
    }

    if (cross) FScross(resl[0],jdat[0],jdat[1],average);
    else for (j=0; j<nd; j++) {
      if (mom) FSmom(resl[j],jdat[j],mom,average[j]);
      else if (binder) FSBinder(resl[j],jdat[j],average[j]);
      else if (isfree) FSfree(resl[j]);
      else FSval(resl[j],jdat[j]);
    }

    for (j=0; j<nv; j++) for (k=0; k<bnum; k++) {
      res[j][k] += resl[j][k];
      sig[j][k] += sqr(resl[j][k]);
    }

    /* write stuff */

    if (jprint) for (i=0; i<bnum; i++) {
      fprintf(out,"%d %d  %.10lg  ",iter,i,p_beta(betav[i]));
      if (isact2 && !is_su3h) fprintf(out,"%.12lg  ",(beta2v[i]));
      for (j=0; j<nv; j++) fprintf(out,"%.13lg ",resl[j][i]);
      fprintf(out,"\n");
    }
  }
  if (jprint && !jstdout) fclose(out);

  for (j=0; j<nv; j++) for (k=0; k<bnum; k++) res[j][k] /= jack;

  for (j=0; j<nv; j++) for (k=0; k<bnum; k++)
    sig[j][k] = sqrt((jack-1)*(sig[j][k]/jack - sqr(res[j][k])));

  if (!jstdout) for (i=0; i<bnum; i++) {
    fprintf(stdout,"%.10lg  ",p_beta(betav[i]));
    if (isact2 && !is_su3h) fprintf(stdout,"%.12lg  ",(beta2v[i]));
    for (j=0; j<nv; j++) 
      fprintf(stdout,"%.14lg %.14lg   ",res[j][i],sig[j][i]);
    fprintf(stdout,"\n");
  }

  if (!isfree) {
    for (j=0; j<nd; j++) {
      for (i=0; i<n; i++) free(jdat[j][i]);
      free(jdat[j]);
    }
    free(jdat);
  }

  for (j=0; j<nv; j++) {
    free(res[j]); free(sig[j]); free(resl[j]);
  }
  free(res); free(sig); free(resl);
  for (i=0; i<nd; i++) if (average[i] != NULL) free(average[i]);
}

/*****************************/

void
FShgvecJ(double **dv[],int jack,int jprint,char *js,
	 char *ss,int bins,double cut_value,int eq_weight,int eq_height)
{

  int i,j,k,ntop,nbot,iter;
  FILE *out;
  double *jp, *dp;
  double ave, sig, **resl,**jdat,**bv;

  /* allocate workspace, if needed */

  jdat = (double **)calloc(n,sizeof(double *));
  for (i=0; i<n; i++) jdat[i] = dblarr(nn[i] - nn[i]/jack + 1);
  bv   = (double **)calloc(jack,sizeof(double *));
  resl = (double **)calloc(jack,sizeof(double *));
  for (i=0; i<jack; i++) bv[i] = dblarr(bnum);
  for (i=0; i<jack; i++) resl[i] = dblarr(bnum);

  if (jprint) if (jstdout) out = stdout; else out = fopen(js,"w");

  for (iter=0; iter<jack; iter++) {
    FSiniJack(iter,jack,ss);

    for (i=0; i<n; i++) {
      nbot  = iter*nn[i]/jack;
      ntop  = (iter+1)*nn[i]/jack;

      jp = jdat[i];
      dp = dv[0][i];
#pragma ivdep
      for (j=0; j<nbot; j++) jp[j] = dp[j];
#pragma ivdep
      for (j=ntop; j<nn[i]; j++) jp[j-ntop+nbot] = dp[j];
    }

    FShgvect(bv[iter],resl[iter],&jdat,bins,cut_value,eq_weight,eq_height);

    /* write stuff */

    if (jprint) for (i=0; i<bnum; i++) {
      fprintf(out,"%d %d  %.10lg  ",iter,i,p_beta(bv[iter][i]));
      if (isact2 && !is_su3h) fprintf(out,"%.12lg  ",(beta2v[i]));
      fprintf(out,"%.13lg\n",resl[iter][i]);
    }
  }
  if (jprint && !jstdout) fclose(out);

  if (!jstdout) for (i=0; i<bnum; i++) {
    for (j=0; j<jack; j++) ave += bv[j][i];
    ave /= jack;
    for (j=0; j<jack; j++) sig += sqr(bv[j][i] - ave);
    sig = sqrt((jack-1)*sig/jack);
    fprintf(stdout,"%.12lg  %.12lg ",p_beta(ave),sig);

    if (isact2 && !is_su3h) fprintf(stdout,"%.12lg  ",(beta2v[i]));

    for (j=0; j<jack; j++) ave += resl[j][i];
    ave /= jack;
    for (j=0; j<jack; j++) sig += sqr(resl[j][i] - ave);
    sig = sqrt((jack-1)*sig/jack);
    fprintf(stdout,"%.14lg  %.14lg\n",ave,sig);
  }

  for (i=0; i<n; i++) free(jdat[j]);
  free(jdat);
  for (j=0; j<jack; j++) free(bv[j]);   free(bv);
  for (j=0; j<jack; j++) free(resl[j]); free(resl);
}


/**************/

void 
FSfree(double free[]) {

  int i;

  for (i=0; i<bnum; i++) free[i] = -logsum[i];
}


/**************************************************************
 *  FSBindersearch  - searches for the max of Binder cumulant
 *  and the ising value for it
 */

double
FSBindersearch(double **dat,double *betain, double beta2) 
{
  int i,j,k,iter;
  double bval,loc,beta;
  FILE *out;
  double *op, *jp, *dp, *sig, val[3], p[3], pe[3];

  bnum = betanumber = 3;
  loc = *betain;

  iter = 0;
  do {
    iter++;
    beta = loc;

    betav[0] = beta*(1.0-1e-5);
    betav[1] = beta;
    betav[2] = beta*(1.0+1e-5);
  
    beta2v[0] = beta2v[1] = beta2v[2] = beta2;

    FSinimeas();  /* logsum has at least 10 elements */
    FSBinder(val,dat,(double *)NULL);

    /* fit parabola */
    polyfit(3,betav,val,NULL,2,p,pe);
    /* location of max */
    loc = -p[1]/(2*p[2]);
    bval = p[0] + p[1]*(loc) + p[2]*sqr(loc);

    /* fprintf(stderr,"beta %.11g, loc %.11g val %.11g\n",beta,loc,bval);
     */
  } while (fabs(loc-beta) > fabs(1e-9*beta) && iter < 30);

  /*  fprintf(stderr,"Binder norm %.11g\n",fabs(loc-beta)); */

  if (iter >= 30) 
    halt("Bindersearch: no convergence, loc %g\n",loc);
  
  *betain = loc;
  return(bval);
}

/**********************************************************
 *  calculate the second moment eigenmatrix
 */

void
FSm2mat(double **dv[],int nd,double *val)
{
  int i,j,k,r,p,cutbuf;
  double *ap,*a2p,*wp,*dp,beta,beta2,dval,*ave;

  if (bnum > 1) 
    halt("Only one beta-value!\n",NULL);
  
  for (i=0; i<nd*nd; i++) val[i] = 0;

  ave = dblarr(nd);

  cutbuf = is_cut; is_cut = 0;
  for (i=0; i<nd; i++) FSval(ave+i,dv[i]);
  is_cut = cutbuf;

  for (i=0; i<n; i++) {
    ap = a[i];
    a2p= a2[i];
    wp = wt[i];
    dp = down[i];
    
    beta = betav[0];
    beta2= beta2v[0];
    map_beta(&beta,&beta2);
    dval = logsum[0] - log(g[i]);

    if (!isact2) { 
      for (j=0;j<nv[i];j++) {
	for (r=0; r<nd; r++) for (p=0; p<=r; p++)
	  val[r+nd*p] += 
	    (dv[r][i][j] - ave[r])*(dv[p][i][j] - ave[p])*
	    exp(-beta*ap[j] - dp[j] - dval);
      }
    } else {
      for (j=0;j<nv[i];j++) {
	for (r=0; r<nd; r++) for (p=0; p<=r; p++)
	  val[r+nd*p] += 
	    (dv[r][i][j] - ave[r])*(dv[p][i][j] - ave[p])*
	    exp(-beta*ap[j] - beta2*a2p[j] - dp[j] - dval);
      }      
    }
  }
  free(ave);

  for (r=0; r<nd; r++) for (p=r+1; p<nd; p++) val[r+nd*p] = val[p+nd*r];
}



/**********************************************************
 *  calculate the second moment eigenmatrix
 */

void
FSeigJack(double **dv[],int nd,int jack,int jprint, char *js, char *ss)
{
  int i,j,k,ntop,nbot,iter;
  FILE *out;
  double beta,beta2, *jp, *dp;
  double *resm,*eval,*evals,*evec,*evecs,*evecj,*evalj, ***jdat;

  /* allocate workspace, if needed */

  if (bnum != 1) 
    halt("One and only one beta-value!",NULL);

  eval = dblarr(nd); evals = dblarr(nd); evalj = dblarr(nd);
  evec = dblarr(nd*nd); evecs = dblarr(nd*nd); evecj = dblarr(nd*nd);
  resm = dblarr(nd*nd); 

  jdat = (double ***)calloc(nd,sizeof(double ***));
  for (j=0; j<nd; j++) {
    jdat[j] = (double **)calloc(n,sizeof(double **));
    for (i=0; i<n; i++) 
      jdat[j][i] = dblarr(nn[i] - nn[i]/jack + 1);
  }

  if (jprint) if (jstdout) out = stdout; else out = fopen(js,"w");

  for (iter=0; iter<jack; iter++) {
    FSiniJack(iter,jack,ss);

    for (i=0; i<n; i++) {
      nbot  = iter*nn[i]/jack;
      ntop  = (iter+1)*nn[i]/jack;

      for (k=0; k<nd; k++) {
	jp = jdat[k][i];
	dp = dv[k][i];
#pragma ivdep
	for (j=0; j<nbot; j++) jp[j] = dp[j];
#pragma ivdep
	for (j=ntop; j<nn[i]; j++) jp[j-ntop+nbot] = dp[j];
      }
    }

    FSm2mat(jdat,nd,resm);
    jacobi(nd,resm,evalj,evecj,1);

    for (k=0; k<nd; k++) {
      eval[k]  += evalj[k];
      evals[k] += sqr(evalj[k]);
    }
    for (k=0; k<nd*nd; k++) {
      evec[k] += evecj[k];
      evecs[k] += sqr(evecj[k]);
    }

    /* write stuff */

    if (jprint) {
      for (i=0; i<nd; i++) {
	fprintf(out,"%d %d eigenvalue %.14lg\n",iter,i,evalj[i]);
	fprintf(out,"%d %d vector ",iter,i);
	for (j=0; j<nd; j++) fprintf(out,"%.12lg ",evecj[i*nd+j]);
	printf("\n");
      }
    }
  }
  if (jprint && !jstdout) fclose(out);

  for (k=0; k<nd; k++) {
    eval[k] /= jack;
    for (j=0; j<nd; j++) evec[k+nd*j] /= jack;
  }

  for (k=0; k<nd; k++) {
    evals[k] = sqrt((jack-1)*(evals[k]/jack - sqr(eval[k])));
    for (j=0; j<nd; j++) 
      evecs[k+nd*j] = sqrt((jack-1)*(evecs[k+nd*j]/jack - sqr(evec[k+nd*j])));
  }
  if (!jstdout) for (i=0; i<nd; i++) {
    printf("%d eigenvalue %.14lg  %.14lg\n",i,eval[i],evals[i]);
    printf("%d vector\n",i);
    for (j=0; j<nd; j++)
      printf("%d %.12lg  %.12lg\n",j,evec[i*nd+j],evecs[i*nd+j]);
  }
}
  
/*******************************************************************
 *  Coupling constant mapping functions
 */

/* First, `old' Higgs */

#define G     (2./3.)               /* Gauge coupling */
#define G2    (G*G)
#define MW    (80.6)                /* W mass         */
#define MD    (sqrt(5./6.)*G)       /* Debye mass     */

double
beta2_Higgs(double mh,double betah)
{
  return(sqr(betah)/betag *
	 ((1./8.)*sqr(mh/MW) - 3*G2/(128*pi*MD))/(1. - G2/(24*pi*MD)));
}

/* Then, new u1 x su2 -mapping */

#define Sigma 3.1759115

double
beta_u1su2(double p_y,double p_x)
{
  double g2,betau1,betah,beta2,p_z;


  if (!isact2) p_x = p_x_global;
  p_z = 0.3;         /* fixed now */
  g2  = 4.0/betag;   /* g2 a */

  betau1 = betag/p_z;
  betah  = 2.0 * g2;                /* gives phi_L = phi/sqrt(g2) */
  /*  beta4  = p_x * g2*g2*g2;  */
  beta2  = betah * (8*p_y/sqr(betag) -  /* +3 used to be here */
		    Sigma*(3 + 12*p_x + p_z)/(4*pi*betag) -
		    1.0/(2*sqr(pi*betag)) *
		    ((51.0/16 - 9*p_z/8 - 5*sqr(p_z)/16 + 
		      9*p_x - 12*sqr(p_x) + 3*p_x*p_z)*(log(3*betag/2) + 0.09)
		     + 5 - 0.9*p_z + 0.01*sqr(p_z)
		     + 5.2*p_x + 1.7*p_x*p_z));

  beta2  *= sqr(betag)/(8*betah); /* this, since we operate on `betay',
				     which factorizes out this number */
  /* betay  = betah*8/sqr(betag); */

  return(beta2);
}

/*******************************************************************
 *  mapping function for su3h along the physical curve
 *  converts parameter y to beta2 and beta4
 */


double 
beta_su3h(double p_y,double *b4)
{
  double p_x,betaA,beta2,beta4;

  /* first, get x-y pair */
  p_x = 3.0/( 8*sqr(pi) * ( p_y - 9.0/(16*sqr(pi))) );

  /* unimprove it */
  p_x = p_x + ( 0.328432 - 0.835282 * p_x + 1.167759 * sqr(p_x) ) / betag;

  /* copied directly from su3h-code 6/2000 */
  betaA = 12/betag;
  beta4 = p_x * 1.5 * sqr(betaA)/betag;
  beta2 = 3*betaA * (1 + 6*p_y/sqr(betag) 
		     - (6 + 10*p_x)*Sigma/(4*pi*betag)
		     - 6/(16*sqr(pi*betag)) * 
		     ((60*p_x - 20*sqr(p_x))*(log(betag) + 0.08849)
		      + 34.768*p_x + 36.130));
  /* now beta4, beta2 multiply directly a2, a4 - which have
     been multiplied by volume */

  *b4 = beta4;
  return( beta2 );
}

/******************************************
 * finally, routine which transforms the couplings 
 */

void 
map_beta(double *beta, double *beta2)
{
  if (Higgs) *beta2 = beta2_Higgs(*beta2,*beta);
  else if (is_u1su2) *beta = beta_u1su2(*beta,*beta2);
  else if (is_su3h) *beta = beta_su3h(*beta, beta2);
}


/*********************************************************************/

void
errorcalc(int nconf,double svect[],double *avep,double *sigp,
	  double *fsp, double *fser,double *w) {

#define TINT 6

  /*
   *     this subroutine analyses the data given in svect
   *     nconf is the number of the data
   *     ave   gives the average value
   *     sig   gives the naive error
   *     fsp   integrated autocorrelation time
   *     fser  error of fsp
   *
   *     TRUE SIGMA = sig*sqrt(2*fsp)
   */

  double ave,fs,sig,ax,ay,fi,norm;
  int i,j,nc,it,nis,ns,nm ,k;

  ave = sig = norm = 0.0;

  if (w == NULL) {
    for(i=0; i<nconf; i++) ave += svect[i];
    ave /= nconf;
    for(i=0; i<nconf; i++) sig += sqr(svect[i]-ave);
    sig /= nconf;
  } else {
    double t;
    for(i=0; i<nconf; i++) {
      t = exp(-w[i]);
      ave += svect[i]*t;
      norm += t;
    }
    ave /= norm;
    for(i=0; i<nconf; i++) sig += sqr(svect[i]-ave)*exp(-w[i]);
    sig /= norm;
  }
  if (sig < 1.0e-12*sqr(ave) || sig < 1.0e-12 || noauto) {
    if (!noauto) fprintf(stderr," ** small sigma/ave: %lg/%lg\n",sqrt(sig),ave);
    *avep = ave;
    *sigp = sqrt(sig/nconf);
    *fsp  = *fser = 1.0;
    return;
  }

  /* time correlations */

  it = 0;
  fs = 0.5;
  do {
    it++;
    ax = ay = fi = 0.0;

    nc=nconf-it;
    for (j=0; j<nc; j++) {
      fi += svect[j]*svect[j+it];
      ax += svect[j];
      ay += svect[j+it];
    }

    fi = (fi-ax*ay/nc)/(nc*sig);
    fs = fs+fi*nc/nconf;
  } while (TINT*fs > it && it < nconf/2);

  if (it >= nconf/2) fprintf(stderr," ** correlation > N/2*%d\n",TINT);

  sig=sqrt(sig/(nconf-1));

  *sigp = sig;
  *avep = ave;
  *fsp  = fs;
  *fser = fs*sqrt(2.*(2.*it+1)/nc);

  return;
}



/**************************************************************************
 *
 *     Normal I/O system
 *
 *     Kari Rummukainen July 1994
 *
 *     reads formatted datafiles (with fixed length records!)
 */

static int n_data = 0;
static FILE *in_file;
static char file_name[500];
static int file_counter;


int 
inicfil(char *fn) {

  char dstr[10000], *cp, *tp;
  double dfloat;
  int i;

  if ((in_file = fopen(fn,"r")) == NULL) 
    halt(" *** could not open file \'%s\'\n",fn);

  fgets(dstr,10000,in_file);
  i = 0;
  cp = dstr;
  do {
    dfloat = strtod(cp,&tp);
    if (tp == cp) break;
    i++;
    cp = tp;
  } while (1);

  rewind(in_file);
  n_data = i;
  return(i);
}

void
getcdbldat(double d[]) 
{
  int i;
  char dstr[10000], *cp, *tp;

  fgets(dstr,10000,in_file);

  cp = dstr;
  for (i=0; i<n_data; i++) {
    d[i] = strtod(cp,&tp);
    if (tp == cp) {
      fprintf(stderr,"Element %d in input file \'%s\' ?\n",file_counter,file_name);
      exit(0);
    }
    cp = tp;
  }
}

int 
closefil() {
  fclose(in_file);
  return(1);
}

/**************************************************************************
 *
 *     Unformatted I/O system - use io_unformat.c
 *
 */


int 
inifil(char *fn,double *inv_vol) 
{

  file_counter = 0;
  strcpy(file_name,fn);

  if (is_ascii) n_data = inicfil(fn);
  else {
    if ((in_file = fopen(fn,"r")) == NULL) 
      halt(" *** could not open file \'%s\'\n",fn);
    n_data = readheader(in_file,&h);
  }

  *inv_vol = 1.0/(h.lx * h.ly * h.lz * h.lt);
  return(n_data);
}

void
getdbldat(double d[]) 
{
  long k;
  int i;

  if (is_ascii) getcdbldat(d);
  else {
    if (++file_counter != readdata( in_file,d )) {
      fprintf(stderr,"Element %d in input file \'%s\' ?\n",file_counter,file_name);
      exit(-1);
    }
  }
}    


