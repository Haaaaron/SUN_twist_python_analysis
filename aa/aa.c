/*************************************************************************
 *
 * This c-program analyzes autocorrelations from
 * histogram datafiles
 *
 * aa [opt] 'name'
 * input files:
 *        histogram file
 * output to standard-output
 *
 * COMPILE:
 *  cc -o aa aa.c io_unformat.c read_ascii.c calclist.c dblarr.c jacobi.c halt.c -lm
 *
 *************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
/* #include <float.h> */
#include <memory.h>

#include "stuff.h"
 
/* #include <io.h> */
  
/* #define max(x,y) (((x) > (y)) ? (x) : (y))  */
 
/* prototypes... 
 */
// alloc in units of 60000 measurements
#define ALLOCATE_UNIT 60000

#define MAX_DATA 10000

#define PI 3.14159265358979
#define pi2 (PI*2.0)

double * dblalloc(int num);
void halt(char *message,char *msg);
int ccorr_matrix(double *dat[],int n,int length,int jack,
		 double *wd,int weight,double wnorm);
int history(int nv, double a[],int showlevel,double level_val,
	    double wd[],int weight);
int errorcalc(int nconf,double svect[],double *avep,double *sigp,
	       double *fsp,double *fser);
int autocorr(double *d,int nd,double *res,int nres);
double timecorr(double d[],int nd,double res[],int nres);
int histgr2(double *a,double *b,int ndata,int bins,double *wd,int weight,
	    double low1,double high1,double low2,double high2);
double tunnelcalc(int nd,double *d,double tmin,double tmax);
int printhgram(double *p,int bins,int ndata,double *wd,int weight,int jack,double hbin);
int ftrans(double *d,int nd);
int printtimeh(double *p,int bin,int ndata);
double dissipation_coeff(double d[],int nd,double res[],int nres);


int noauto = 0;
int ascii = 0;
char * label = NULL;

char usage[]  = " usage: aa [-opt] datafile\n\
 Analyses measurement data files, with formats\n\
 plain ascii/labelled ascii/special binary/\n\
 Options:\n\
  +LABEL         : use labelled ascii datafiles (no - -sign!)\n\
  ++             : force plain unlabelled ascii [DEFAULT, if invalid binary]\n\
  q              : quiet\n\
  d num          : use only data num\n\
  D #2+#3,#4     : collect data and do the arithmetic [+-*/]\n\
  I n            : use only measurements with index divisible by n\n\
  j nblocks      : jackknife data\n\
  J nblocks      : print jackknife blocks\n\
  A              : print all measurements as jackknife blocks\n\
  i iters        : use only i iterations\n\
  s iters        : skip iteration\n\
  S n            : use only every n iteration\n\
  R r1:r2        : print data from r1 to r2\n\
  f num=val      : filter data with vector num = value\n\
  l[L] a:b       : data range from a to b [element 2]\n\
  X              : print just the data indicated\n\
  C              : cross-correlation of data\n\
  E              : eigenvalue analysis of cross-correlation matrix\n\
  c num          : correlation number num\n\
  cl length      : force vector length\n\
  cc/cd num      : print 2-dim. array number num\n\
  cos            : cosine transform the data (only with -j, -A)\n\
  a              : sum -D or -c into one data element\n\
  x d,i,num      : datavector element d[index[num]], num optional\n\
    K n          : take modulo of index\n\
  h              : print histogram plot\n\
  v value        : show value in the histogram\n\
  w              : use the weight factor (elem 1)\n\
  0              : subtract the 0-moment (average)\n\
  2              : print second moment\n\
  p power        : print moment power\n\
  F bins         : print histogram in weight format\n\
  r block        : do the 'running blocking'\n\
  b block        : do the 'normal blocking'\n\
  m/M value      : analyze only up to minimum/maximum value\n\
  O              : list measurements up to minimum/maximum value\n\
  P              : lag the list by one unit\n\
  g [bins]       : print histogram to the standard output, with bins bins\n\
  B [binw]       : print histogram, with binwidth [1]\n\
  H              : print history to the standard output\n\
  G bins         : 2-component histogram\n\
  W b_old,b_new,\'act_string\'  : reweight to new beta value (up to 10 times)\n\
  y tmin:tmax    : print tunneling time\n\
  t              : no autocorrelations\n\
  T length       : print autocorrelation function for distance length\n\
  N length       : print autocorrelation (non-normalized)\n\
  Q length       : measure <(x(t+dt)-x(t))^2> to distance dt=length";

// Old, not used any more
//  Y bg,mU_old,mU_new,T_old,T_new : reweight to new T value (susy)\n   \
//  Z mH_old,mH_new,b_old,b_new : reweight to new mH, betaH -values (su2-Higgs)\n \

 
int main(argc,argv)
     int argc;
     char * argv[];
 
{
  double *fp,*dat[MAX_DATA],*tmparr,*wd;
  double *tmpx;
  int block,indexed,iv,dv,index;
  int nblocking,length,jack,sum;
  double temps,minval,maxval,tmin,tmax;
  char * ss,*lists[MAX_DATA];
  int raw,eig;
  int idata,nd;
  int i,j,jj,k,iters,hist,printhist;
  long lk;
  int gram,rblock,bblock,skip,weight,binder=0;
  double hbin,error,naive,aveg,average[MAX_DATA],timerel,level_val;
  int list,limits,lagged,showlevel,stagger,mom2,ccorr;
  int wgram,irange1,irange2,icorr,printjack,index_mod;
  int four,h2gram,reweight,isarray,tunnel;
  int quiet,sqrtV,icorrlen;
  int opt_on = 0;
  int filter;
  double filtervalue;
  double betar,wnorm;
  int autocorrelation = 0,sub_ave = 0;
  int timecorrelation = 0;
  int dissipation = 0;
  double power = 0,beta1[10],beta2[10];
  char rwstring[10][100];
  FILE *ff;
  e_header h;
  double low1,high1,low2,high2;
  char *legend = NULL;  // legend string for ascii
  int data_index = 1;
  // int susy,higgs,
  // double susy_mul,betag;
  // double mU,mUn,mH,mH2,T,Tn;

  hbin = 0.0;
  block = nblocking = iters = hist = printhist = gram = idata =
    rblock = bblock = jack = list = limits = stagger = mom2 = quiet = 0;
  lagged = showlevel = weight = wgram = 0;
  irange1 = irange2 = icorr = icorrlen = 0;
  ccorr = sqrtV = sum = four = h2gram = 0;
  printjack = indexed = index_mod = reweight = 0;
  raw = eig = isarray = 0;
  filter = tunnel = 0;
  nd = 0;
  minval = -1e-60;
  maxval = 1e60;
  
  low1=low2=high1=high2=0;

  skip = 0;
  if (argc <= 1) halt(usage,NULL);

  for (i=0; i<MAX_DATA; i++) dat[i] = NULL;

  while (--argc > 0 && 
	 ( (*++argv)[0] == '-' || (*argv)[0] == '+') ) {

    ss = argv[0] + 1;

    /* check '+' */
    if (*(ss-1) == '+') {
      if (*ss == '+') {
	ss++;
	label = NULL;
	ascii = -1;
      } else {
	getstring(label);
	ascii = 1;
      }

    } else {
      opt_on = 1;
      while (*ss) {
	switch(*ss++) {

	case 'q': quiet   = 1; break;
	case 'O': list    = 1; break;
	case 'P': lagged  = 1; break;
	case 't': noauto  = 1; break;
	case 'T': getnum("%d",&autocorrelation); break;
	case 'N': getnum("%d",&timecorrelation); break;
	case 'Q': getnum("%d",&dissipation); break;
	case 'y': tunnel  = 1; get2num("%lg:%lg",&tmin,&tmax); break;
	case 'w': weight  = 1; break;
	case '0': sub_ave = 1; break;
	case '2': mom2    = 1; break;
	case 'a': sum     = 1; break;
	case 'p': getnum("%lg",&power); break;
          //case '': binder  = 1; break;
	case 'F': wgram = weight = 1; idata = 4; 
	  getnum("%d",&gram); break;
	case 'c':
	  if (strcmp(ss,"os") == 0) { four = 1; ss+=2; }
	  else if (*ss == 'l') { ss++; getnum("%d",&icorrlen); }
	  else {
	    if (*ss == 'c') { isarray=1; ss++; }  
	    if (*ss == 'd') { isarray=2; ss++; }
	    getnum("%d",&icorr);
	  }
	  break;
	case 'm': limits  = 1; getnum("%lg",&minval); break;
	case 'M': limits  = 1; getnum("%lg",&maxval); break;
	case 'i': getnum("%d",&iters); break;
	case 's': getnum("%d",&skip); break;
	case 'S': getnum("%d",&stagger); break;
	case 'r': getnum("%d",&rblock); break;
	case 'b': getnum("%d",&bblock); break;
	case 'h': hist = 1; break;
	case 'H': printhist = 1; break;
	case 'v': showlevel = 1; getnum("%lg",&level_val); break;
	case 'j': getnum("%d",&jack); break;
	case 'J': printjack = 1; getnum("%d",&jack); break;
	case 'A': printjack = 1; jack = -1; break;
	case 'W': 
	  /* take away white space */	
	  { char *p,*q;
	    p = q = ss;
	    while (*p) { 
	      if (*p != ' ' && *p != '\t') *(q++) = *p;
	      p++;
	    }	
	    *q = 0;
	  }
	  get3num("%lg,%lg,%s",&beta1[reweight],&beta2[reweight],rwstring[reweight]);
	  beta2[reweight] -= beta1[reweight]; 
	  reweight++;
	  break;

          //	case 'Y': reweight = susy = 1;
          //	  get5num("%lg,%lg,%lg,%lg,%lg",&betag,&mU,&mUn,&T,&Tn);
          //	  break;
          //	case 'Z': reweight = higgs = 1;
          //	  get4num("%lg,%lg,%lg,%lg",&mH,&mH2,&beta1[0],&beta2[0]);
          //	  break;
	  
	case 'd': getnum("%d",&idata); nd = 1; break;
	case 'D': getlist(lists,nd); break;

        case 'I': getnum("%d",&data_index); break;
	  
	case 'f':
	  get2num("%d=%lg",&filter,&filtervalue);
	  break;
	  
	case 'X': raw = 1; break;
	case 'E': eig = 1; break;
	  
	case 'R':
	  if (!(*ss) && --argc) ss = (++argv)[0];
	  if (sscanf(ss,"%d:%d",&irange1,&irange2) != 2) halt(usage,NULL);
	  ss = strchr(ss,0);
          nd = irange2-irange1+1;
	  break;

	case 'x':
	  indexed = 1;
	  if (!(*ss) && --argc) ss = (++argv)[0];
	  if ((i = sscanf(ss,"%d,%d,%d",&dv,&iv,&index)) == 3) {
	    irange1 = irange2 = index;
	  } else if (i == 2) index = 0;
	  else halt(usage,NULL);
	  ss = strchr(ss,0);
          nd = 1;
	  break;
       
	case 'K': getnum("%d",&index_mod); break;
	  
	case 'C': ccorr = 1; break;

	case 'B': getoptnum("%lg",&hbin,1.0); gram = 1; break;
	case 'g': getoptnum("%d",&gram,100); break;
	case 'G': getnum("%d",&h2gram); break;
	
	case 'l': get2num("%lg:%lg",&low1,&high1); break;
	case 'L': get2num("%lg:%lg",&low2,&high2); break;

	default: halt(usage,NULL);

	} /* switch */
      } /* while */
    } /* else '-' */
  } /* while */

  if (irange1 && sum) nd = 1;
  if (idata) nd = 1;

  if (argc == 0) halt(usage,NULL);

  int ndata;
  length = 0;
  int allocated = 0;

  // over files
  int filep;
  for (filep=0; filep<argc; filep++) {

    ss = argv[filep];
    lk = 0;

    ff = fopen(ss,"r");
    if (ff == NULL) halt("Could not open file %s",ss);
  
    if (!ascii) {
      int bl = readheader(ff,&h);
      if (filep == 0) block = bl;
      // set ascii if not binary
      // fprintf(stderr,"* Defaulting to ascii\n");
      if (bl < 0) ascii = -1;
    }
    // only for the 1st file 
    if (filep == 0) {
      if (ascii) {
        int rows;
        h.n_double = block = read_ascii_info(ff,label,&rows,&legend,data_index);
        h.n_float = h.n_long = h.n_char = 0;
        if (rows > 1) {
          h.d2 = rows;
          h.d1 = block / h.d2;
        }
        if (legend != NULL) 
          fprintf(stderr,"- Legend: %s\n",legend); 
      }
      
      if (nd) ndata = nd;
      else ndata = block;

      if (weight || reweight) {
        wd = NULL;
        ndata++;
      }

      if (icorrlen) h.d2 = icorrlen;
      tmparr = dblarr(block);

      if (indexed && !irange1) icorr = 1;

      if (icorr) {
        if (h.d1 >= icorr || icorrlen) {
          irange1 = 1 + (icorr-1)*h.d2;
          irange2 = irange1 + h.d2 - 1;
        } else if (h.d1 + h.d3 >= icorr) {
          irange1 = 1 + h.d1*h.d2 + (icorr-1-h.d1)*h.d4;
          irange2 = irange1 + h.d4 - 1;
        } else halt("No such thing!",NULL);
      }

      if (indexed) tmpx = dblarr(irange2-irange1+1);

    }  // fp == 0, the first file

    /*
    do {
      length++;
      if (!ascii) lk = readdata(ff,tmparr);
      else {
        if (read_ascii(ff,label,tmparr)) lk = length;
        else lk = 0;
      }
    } while (lk == length && iters != length);
    if (lk != length) length--;
    rewind(ff);
    if (!ascii) skipheader(ff);
    */

    if (icorrlen) h.d2 = icorrlen;
    
    // length -= (skip);

    /*
    
    if (length <= 0) halt("Measurements == 0",NULL);
    */
    if ((ccorr || h2gram) && nd < 2) 
      halt("at least 2 datavectors required",NULL);

    if (idata > block) halt("Illegal data number",NULL); 
    if (irange1 > irange2 || irange2 > block) halt("Illegal range spec",NULL);

    // if (susy) susy_mul = pow(4.0/(betag * 4.0/9.0),3.0) * (h.lx*h.ly*h.lz);

    j = skip;
    while (j-- > 0) {
      if (!ascii) lk = readdata(ff,tmparr);
      else read_ascii(ff,label,tmparr,data_index);
    }

    // MAIN READ LOOP STARTS
    for ( ; ; length++) {
      j = length;

      if (iters && iters == length) break;

      if (!raw) {      
        if (length >= allocated) {
          allocated += ALLOCATE_UNIT;
          for (i=0; i<ndata; i++) {
//            fprintf(stderr,"REALLOC: size %d\n",allocated);
            dat[i] = (double *)realloc( dat[i], allocated * sizeof(double) );
          }
          if (weight || reweight) wd = dat[ndata-1];
        }
      }

      extern int calclist_index;

      if (!ascii) {
        if (!readdata(ff,tmparr)) break;
      } else {
        if (!read_ascii(ff,label,tmparr,data_index)) break;
      }
    
      calclist_index = length;

      k=0;
      if (!filter || tmparr[filter-1] == filtervalue) {
        
        if (reweight || weight) {
          double w;
          static int first=1;
          static double w0;
	
          /* if (susy) */ 
          /**** THis is the `old' susy weight
                w = 0.379*(sqr(100.0/T) - sqr(100.0/Tn)) * tmparr[8] 
                + 0.849*(sqr(mU/T) - sqr(mUn/Tn)) * tmparr[16];
          ****/

          /***** New susy weight 
                 w = -susy_mul*( (1.0/sqr(T) - 1.0/sqr(Tn)) *
                 ( 18384.1*tmparr[3] - 3984.08*tmparr[4] - 2*1191.72*tmparr[7] 
                 + 2*96.6867*tmparr[8] )
                 - (sqr(mU/T) - sqr(mUn/Tn))*tmparr[14] );
          *****/
          /* else if (higgs) w = (beta2[0]/beta1[0]-1)*tmparr[2] + 
           *	  (sqr(beta2[0]*mH2/(beta1[0]*mH))-1)*tmparr[4];
           * else 
           */
          if (reweight) 
            for (w=i=0; i<reweight; i++) w += beta2[i]*calclist(tmparr,block,rwstring[i]);
          /* beta2[i]*tmparr[rwi[i]-1];*/
          else w = 0.0;
          if (weight) w += tmparr[0];
	
          /* rebalance weight with the 1st element of the vector */
          /* This screws up normalisation of histograms, but should
             not matter in practice */
          if (first) {
            w0 = w;
            first = 0;
          }
          wd[j] = exp(-(w-w0));
        }
      
        if (indexed) {
          int id;
          for (i=0,k=irange1-1; k<irange2; k++,i++) {
            id = (int)tmparr[(iv-1)*h.d2 + k];
            if (index_mod) id = id % index_mod;
            tmpx[i] = tmparr[id];
          }
          for (k=irange1-1; k<irange2; k++) tmparr[k] = tmpx[k-irange1+1];
        }
    
        if (raw) {
          if (stagger == 0 || j % stagger == 0) {
            if (idata) printf("%.14lg\n",tmparr[idata-1]);
            else if (irange1)
              for (i=irange1-1; i<irange2; i++) 
                printf("%d %.14lg\n",i-irange1+1,tmparr[i]);
            else if (nd) {
              for (i=0; i<nd; i++) 
                printf("%.14lg ",calclist(tmparr,block,lists[i]));
              printf("\n");
            }
          }
        } else if (sum && irange1) {
          fp = dat[0];
          for (i=irange1-1; i<irange2; i++) {
            fp[j] += tmparr[i];
          }
        } else if (irange1) {
          for (i=irange1-1; i<irange2; i++) {
            fp = dat[i-irange1+1];
            fp[j] = tmparr[i];
          }
        } else if (idata) {
          fp = dat[0];
          fp[j] = tmparr[idata-1];
        } else if (nd) {
          for (i=0; i<nd; i++) dat[i][j] = calclist(tmparr,block,lists[i]);
        } else {
          for (i=0; i<block; i++) {
            fp = dat[i];
            fp[j] = tmparr[i];
          }
        }
      } /* filtervalue */
    } /* for over lines */

    fclose(ff); // close file
  } // read new file

  if (raw) return(1);

  if (jack < 0) jack = length;

  if (!quiet) {
    fprintf(stderr,"- Data: %d Measurements: %d\n",block,length);
    fprintf(stderr,"- double %ld, float %ld, long %ld, char %ld\n",
	    h.n_double,h.n_float,h.n_long,h.n_char);
  }

  if (sum && irange1) idata = irange2 = irange1 = 1;
  if (irange1) nd = irange2 - irange1 + 1;

  if (reweight) weight = 1;

  if (weight) {
    wnorm = 0;
    for (j=0; j<length; j++) wnorm += wd[j];
    wnorm = length/wnorm; 
  }
  
  if (eig) {
    ccorr_matrix(dat,nd,length,jack,wd,weight,wnorm);
    return(1);
  }
  if (ccorr) {
    double a1,a2,w1,w2,*fp,*cp;

    a1 = a2 = 0;
    fp = dat[0];
    cp = dat[1];
    for (i=0; i<length; i++) {
      a1 += fp[i];
      a2 += cp[i];
    }
    a1 /= length;
    a2 /= length;
    w1 = w2 = 0;
    for (i=0; i<length; i++) {
      w1 += sqr(fp[i] - a1);
      w2 += sqr(cp[i] - a2);
    }
    w1 = sqrt(w1/length); 
    w2 = sqrt(w2/length);
    for (i=0; i<length; i++) 
      fp[i] = (fp[i] - a1)*(cp[i] - a2)/(w1*w2);
    
    nd = 1;
  }

  if (!nd) {
    if (!quiet) 
    fprintf(stderr,"-------------- av ----- error ------- autoc -------- ind.meas\n");
    
    for (i=0; i<block; i++) {
      fp = dat[i];
      if (weight) for (j=0; j<length; j++) fp[j] *= wnorm*wd[j];
      errorcalc(length,dat[i],&(average[i]),&naive,&timerel,&temps);
      error = naive*(double)sqrt(2.0 * fabs(timerel));
 
      printf("%3d: %12.6g %11.5g %10.4g(%9.4g) %9d",
	     i+1,average[i],error,timerel,temps,(int)(length/(timerel+.5)));
      putchar('\n');
    }
    return(0);
  }

  if (h2gram) {
    histgr2(dat[0],dat[1],length,h2gram,wd,weight,
	    low1,high1,low2,high2);
    return(1);
  }

  irange1 = 0; irange2 = nd-1;
  idata = irange1;

  for (idata=irange1; idata<=irange2; idata++) {
    if (sub_ave || mom2 || binder || power) {
      double av = 0;

      if (weight) {
	for (i=0; i<length; i++) av += dat[idata][i]*wd[i]*wnorm;
      } else {
	for (i=0; i<length; i++) av += dat[idata][i];
      }
      av /= length;
      for (i=0; i<length; i++) dat[idata][i] -= av;
    }	
    
    if (mom2) {
      for (i=0; i<length; i++) dat[idata][i] = sqr(dat[idata][i]);
    }

    if (power != 0) {
      double p, ave;

      if (fmod(power,1.0) == 0) {
        
        for (i=0; i<length; i++) dat[idata][i] = pow(dat[idata][i],power);
  
      } else {
	fprintf(stderr,"warning: pow non-integer:  %g\n",power);
	for (i=0; i<length; i++) 
	  dat[idata][i] = pow(fabs(dat[idata][i]),power);
      }
    }
  }

  idata = irange1;
  if (list) {
    double val;
    int isaccept,begmeas,lag_on;
    
    isaccept = 0;
    lag_on = 1;
    for (i=j=0; i<length; i++) {
      val = dat[idata][i];
      if (val <= maxval && val >= minval) {
	if (lagged && !lag_on) lag_on = 1;
	else {
	  j++;
	  if (!isaccept) {
	    isaccept = 1;
	    begmeas = i;
	  } 
	}
      } else {
	lag_on = 0;
	if (isaccept) {
	  isaccept = 0;
	  printf("%d %d\n",begmeas,i-1);
	}
      }
    }
    if (isaccept) printf("%d %d\n",begmeas,length-1);

    if (!quiet) fprintf(stderr," - accepted %d of %d measurements\n",j,length);
    exit(-1);
  }
    
  if (limits) {
    double *dd[MAX_DATA],mval,Mval;

    for (idata=irange1; idata<=irange2; idata++) dd[idata] = dblarr(length);
    for (i=j=0; i<length; i++) {
      mval = Mval = dat[irange1][i];
      for (idata=irange1+1; idata<=irange2; idata++) {
	mval = smaller(mval,dat[idata][i]);
	Mval = greater(Mval,dat[idata][i]);
      }
      if (Mval <= maxval && mval >= minval) {
	for (idata=irange1; idata<=irange2; idata++) 
	  dd[idata][j] = dat[idata][i];
	j++;
      }
    }
    for (idata=irange1; idata<=irange2; idata++) {
      free(dat[idata]);
      dat[idata] = dd[idata];
    }
    if (!quiet) fprintf(stderr," - accepted %d of %d measurements\n",j,length);
    length = j;
    if (jack > length) jack = length;
  }

  if (stagger) {
    double * dd,val;
    
    for (idata=irange1; idata<=irange2; idata++) {
      dd = (double *)calloc(length/stagger+1,sizeof(double));
      for (i=0; i<length; i++) 
	if ((i)%stagger == 0) dd[(i)/stagger] = dat[idata][i];
      free(dat[idata]);
      dat[idata] = dd;
    }
    length = (length-1)/stagger + 1;
    if (jack > length) jack = length;
  }

  if (rblock) {
    double * dd,val;

    for (idata=irange1; idata<=irange2; idata++) {
      dd = (double *)calloc(length-rblock+1,sizeof(double));
      for (val=i=0; i<rblock-1; i++) val += dat[idata][i];
      for (i=rblock-1; i<length; i++) {
	val += dat[idata][i];
	dd[i-rblock+1] = val/rblock;
	val -= dat[idata][i-rblock+1];
      }
      free(dat[idata]);
      dat[idata] = dd;
    }
    length -= rblock - 1;
  }

  if (bblock) {
    double * dd,val;
    
    for (idata=irange1; idata<=irange2; idata++) {
      dd = (double *)calloc(length/bblock+1,sizeof(double));
      i = k = 0;
      while (i < length) {
	for (val=j=0; j<bblock && i < length; j++) val += dat[idata][i++];
	dd[k++] = val/j;
      }
      free(dat[idata]);
      dat[idata] = dd;
    }
    length = k;
    if (jack > length) jack = length;
  }

  idata = irange1;

  if (printjack && gram) {
    printhgram(dat[0],gram,length,wd,weight,jack,hbin);
    return(1);
  }

  if (printjack || jack) {
    double *jave,*jnorm;
    int dl,id;
    
    dl = irange2-irange1+1;
    jave = dblarr(dl*jack);
    jnorm = dblarr(jack);

    for (j=0; j<jack; j++) {
      if (weight) {
	for (i=0; i<j*length/jack; i++) jnorm[j] += wd[i];
	for (i=(j+1)*length/jack; i<length; i++) jnorm[j] += wd[i];
      }
      else jnorm[j] = (j*length/jack + (length - (j+1)*length/jack));

      for (idata=irange1; idata<=irange2; idata++) {
	id = idata-irange1;
	fp = dat[idata];
	if (!weight) {
	  for (i=0; i<j*length/jack; i++) jave[id+j*dl] += fp[i];
	  for (i=(j+1)*length/jack; i<length; i++) jave[id+j*dl] += fp[i];
	  jave[id+j*dl] /= jnorm[j];
	} else {
	  for (i=0; i<j*length/jack; i++) jave[id+j*dl] += fp[i]*wd[i];
	  for (i=(j+1)*length/jack; i<length; i++) jave[id+j*dl] += fp[i]*wd[i];
	  jave[id+j*dl] /= jnorm[j];
	}
      }
      if (four && dl > 1) ftrans(jave+j*dl,dl);
      if (printjack) {
	for (idata=irange1; idata<=irange2; idata++) {
	  id = idata-irange1;
	  if (!isarray) printf("%d  %.16lg\n",id,jave[id+j*dl]);
	  else if (isarray == 1)
	    printf("%ld  %ld  %.16lg\n",id/h.d3,id%h.d3,jave[id+j*dl]);
	  else {
	    printf("%.16lg ",jave[id+j*dl]);
	    if ((id+1)%h.d3 == 0) printf("\n");
	  }
	  
	  /* fprintf(stderr,"."); fflush(stderr); */
	}
      }
    }
    if (!printjack) {
      double aveg,error;

      for (idata=irange1; idata<=irange2; idata++) {
	id = idata-irange1;
	aveg = error = 0.0;
	for (j=0; j<jack; j++) aveg += jave[id+j*dl];
	aveg /= jack;
	for (j=0; j<jack; j++) error += sqr(jave[id+j*dl] - aveg);
	error = sqrt((jack-1)*error/jack);

	/* fprintf(stderr,
	   "------------- av ----- error\n");
	   fprintf(stderr,"%d %16.9g %11.5g\n",idata-irange1,aveg,error);
	   */

	if (irange1 == irange2) {
	  if (!quiet) fprintf(stderr,"------------- av ----- error\n");
	  printf("%d  %.16g %.6g\n",irange1,aveg,error);
	} else {
	  if (!isarray) printf("%d  %.16g  %.6g\n",id,aveg,error);
	  else if (isarray == 1)
	    printf("%d  %d  %.16g  %.6g\n",(int)(id/h.d3),(int)(id%h.d3),aveg,error);
	  else {
	    printf("%.16lg ",aveg); if ((id+1)%h.d3 == 0) printf("\n");
	  }
	}
      }
    }
    return(1);
  } 

  if (irange1 == irange2) {
    idata = irange1;

    if (autocorrelation) {
      double *at;
      at = dblarr(autocorrelation);
      autocorr(dat[idata],length,at,autocorrelation);
      for (i=0; i<autocorrelation; i++) printf("%d  %g\n",i,at[i]);
      return(0);
    }

    if (timecorrelation) {
      double *at;
      at = dblarr(timecorrelation);
      fprintf(stderr,"Integrated: %g\n",
	      timecorr(dat[idata],length,at,timecorrelation));
      for (i=0; i<timecorrelation; i++) printf("%d  %g\n",i,at[i]);
      return(0);
    }
    
    if (dissipation) {
      double *at;
      at = dblarr(dissipation);
      dissipation_coeff(dat[idata],length,at,dissipation);
      for (i=0; i<dissipation; i++) printf("%d  %g\n",i,at[i]);
      return(0);
    }

    if (!(printhist || gram || hist)) {
      if (weight) for(j=0; j<length; j++) dat[idata][j] *= wnorm*wd[j];

      if (tunnel) {
	printf("%d %17.10g \n",idata+1,tunnelcalc(length,dat[idata],tmin,tmax)); 
	return(0);
      }

      errorcalc(length,dat[idata],&aveg,&naive,&timerel,&temps);
      error = naive*sqrt(2.0 * fabs(timerel));
      if (!quiet) fprintf(stderr,
       "------------------ av ----- error ------- autoc -------- ind.meas\n");

      if (!(printhist || wgram || gram || hist)) {
	printf("%3d  %20.15g %11.5g %8.4g(%7.2g) %9d\n",
	       idata+1,aveg,error,timerel,temps,(int)(length/(timerel+.5)));
      } else {
	fprintf(stderr,"%3d  %20.15g %11.5g %8.4g(%7.2g) %9d\n",
		idata+1,aveg,error,timerel,temps,(int)(length/(timerel+.5)));
      }

      if (weight) for(j=0; j<length; j++) dat[idata][j] /= wnorm*wd[j];
    }
    if (printhist) printtimeh(dat[idata],printhist,length);
  
    if (wgram) printhgram(dat[idata],gram,length,wd,1,-1,hbin);
    else if (gram) printhgram(dat[idata],gram,length,wd,weight,-1,hbin);

    if (hist) history(length,dat[idata],showlevel,level_val,wd,weight);

  } else {

    if (gram) {
      double *sd, *vd;
      int si;
      sd = dblarr(length*(irange2-irange1+1));
      if (weight) vd = dblarr(length*(irange2-irange1+1));
      si = 0;

      for (idata=irange1; idata<=irange2; idata++) for(j=0; j<length; j++) {
	if (weight) vd[si]   = wnorm*wd[j];
	sd[si++] = dat[idata][j];
      }
      if (gram) printhgram(sd,gram,si,vd,weight,-1,hbin);
      free(sd); if (weight) free(vd);
    } else {

      int id;

      for (idata=irange1; idata<=irange2; idata++) {
	if (weight) for(j=0; j<length; j++) dat[idata][j] *= wnorm*wd[j];

	if (!tunnel) {
	  errorcalc(length,dat[idata],&aveg,&naive,&timerel,&temps);
	  error = naive*sqrt(2.0 * fabs(timerel));
	  id = idata-irange1;
	  if (!isarray) 
	    printf("%3d  %16.9g %11.5g %10.4g(%9.4g) %9d\n",
		   idata,aveg,error,timerel,temps,(int)(length/(timerel+.5)));
	  else if (isarray == 1)
	    printf("%d  %d  %16.9g  %11.5g\n",(int)(id/h.d3),(int)(id%h.d3),aveg,error);
	  else {
	    printf("%.16lg ",aveg);
	    if (!((id+1)%h.d3)) printf("\n");
	  }
	} else {
	  printf("%d %17.10g \n",id+1,tunnelcalc(length,dat[idata],tmin,tmax)); 
	}
      }
    }
  }
  return(1);
}


/*****************************************************/
 
int 
printtimeh(double *p,int bin,int ndata)
{
  int i,j,k;
  double dat;

  if (bin <= 0) bin = 1;
  k = 0;
  for (i=1; i<=ndata/bin; i++) {
    for (j=dat=0; j<bin && k<ndata; j++) dat += p[k++];
    printf("%d %.15lg\n",i,dat/bin);
  }
  return(0);
}


int 
printhgram(double *p,int bins,int ndata,double *wd,int weight,int jack,double hbin)
{
  int i,j,k,bn;
  double *hg;
  double maxv,minv,d,xadd;
  
  maxv = -1e60; minv = -maxv;
  for (i=0; i<ndata; i++) {
    maxv = (maxv < p[i]) ? p[i] : maxv;
    minv = (minv > p[i]) ? p[i] : minv;
  }

  if (hbin <= 0) {
    hg = (double *)calloc(bins,sizeof(double));
    d = (maxv-minv)/(bins-1);
    bn = bins;
    maxv += d/2; minv -= d/2;
    xadd = 0.5;
  } else {
    maxv = hbin * ceil(maxv/hbin);
    minv = hbin * floor(minv/hbin);
    //bins = maxv/hbin - minv/hbin + 1;
    bins = maxv/hbin - minv/hbin+1;
    // bn = bins-1;
    bn = bins;
    hg = (double *)calloc(bins,sizeof(double));
    d = hbin;
    xadd = 0;
  }

  if (jack < 2) {

    if (!weight) 
      for (i=0; i<ndata; i++) hg[(int)((p[i]-minv)*bn/(maxv-minv))] += 1.0;
    else
      for (i=0; i<ndata; i++) hg[(int)((p[i]-minv)*bn/(maxv-minv))] += wd[i];

    for (i=0; i<bins; i++) printf("%g %g\n",minv + d*(i+xadd),hg[i]);
    
  } else {
    for (j=0; j<jack; j++) {
      for (i=0; i<bins; i++) hg[i] = 0;

      if (!weight) {
	for (i=0; i<j*ndata/jack; i++) 
	  hg[(int)((p[i]-minv)*bn/(maxv-minv))] += 1.0;
	for (i=(j+1)*ndata/jack; i<ndata; i++)
	  hg[(int)((p[i]-minv)*bn/(maxv-minv))] += 1.0;
      } else {
	for (i=0; i<j*ndata/jack; i++) 
	  hg[(int)((p[i]-minv)*bn/(maxv-minv))] += wd[i];
	for (i=(j+1)*ndata/jack; i<ndata; i++)
	  hg[(int)((p[i]-minv)*bn/(maxv-minv))] += wd[i];
      }

      for (i=0; i<bins; i++) printf("%g %g\n",minv + d*(i+xadd),hg[i]);
    }
  }

  free(hg);
  return(0);
}


int 
histgr2(double *a,double *b,int ndata,int bins,double *wd,int weight,
	double mina, double maxa, double minb, double maxb)
{
  int i,j,k;
  double *hg,d;
  
  fprintf(stderr,"%g %g %g %g\n",mina,maxa,minb,maxb);

  hg = dblarr(bins*bins);

  if (maxa <= mina) {
    maxa = -1e60; mina = -maxa;
    for (i=0; i<ndata; i++) {
      maxa = (maxa < a[i]) ? a[i] : maxa;
      mina = (mina > a[i]) ? a[i] : mina;
    }
  }
  d = (maxa-mina)/(bins-1);
  maxa += d/2; mina -= d/2;
  
  if (maxb <= minb) {
    maxb = -1e60; minb = -maxb;
    for (i=0; i<ndata; i++) {
      maxb = (maxb < b[i]) ? b[i] : maxb;
      minb = (minb > b[i]) ? b[i] : minb;
    }
  }
  d = (maxb-minb)/(bins-1);
  maxb += d/2; minb -= d/2;

  for (i=0; i<ndata; i++) {
    j = (a[i]-mina)*bins/(maxa-mina);
    k = (b[i]-minb)*bins/(maxb-minb);
    
    if (!weight) hg[j + bins*k] += 1;
    else hg[j + bins*k] += wd[i];
  }

  for (j=0; j<bins; j++) {
    for (i=0; i<bins; i++) {
      printf("%g %g %g\n",
	     mina + (maxa-mina)/(bins-1)*(j+0.5),
	     minb + (maxb-minb)/(bins-1)*(i+0.5),
	     hg[j + bins*i]);
    }
  }    
  return(0);
}



/*********************************************/
 
int 
history(nv,ac,showlevel,level_val,wd,weight)
     int nv,showlevel,weight;
     double ac[],level_val,wd[]; 

#define NLEV 30
#define BINS 78
 
{
  int j,i,k,bins;
  double min, max, an, dave,dmin,dmax,dsig;
  char * sp, level[NLEV+1][BINS+2];
  double *his,hmax;

  bins = (BINS < nv) ? BINS : nv;
  his = (double *)calloc(bins+2,sizeof(double));
  printf("\n- Histogram - number of bins %d\n",bins);
 
  for (i=0; i <= NLEV; i++) {
    sp = level[i];
    for (j=0; j< bins; j++) sp[j] = ' ';
    sp[j] = 0;
  }
 
  min = 1e300;
  max = -min;
  for (j=0; j<nv; j++) {
    min = (min < ac[j]) ? min : ac[j];
    max = (max > ac[j]) ? max : ac[j];
  }
 
  for (j=0; j<bins; j++) his[j] = 0;
  for (j=0; j<nv; j++) {
    k = ((ac[j]-min)/(max-min))*(bins-1);
    if (!weight) his[k] += 1.0; else his[k] += wd[j];
  }
 
  for (j=hmax=0; j<bins; j++) if (hmax < his[j]) hmax = his[ k=j ];
 
  printf("min %g  max %g;  maximum count %g in bin %d, x = %g\n",
	 min,max,hmax,k+1,min+k*(max-min)/(bins-1));
 
  /*  if (histnorm) an = 1.0/nv; else an = 1.0; */
  an = ((double)NLEV)/hmax;
 
  if (showlevel && level_val >= min && level_val <= max) {
    for (j=0; j<NLEV; j++) {
      level[j][ (int) (((level_val-min)/(max-min))*(bins-1)) ] = '|';
    }
  }

  for (j=0; j<bins; j++) {
    sp = level[ (int)( his[j]*an) ];
    sp[ j ] = '*';
  }

  for (j=NLEV; j>=0; j--) printf("%s\n",level[j]);
  
  free(his);
 
  for (i=0; i <= NLEV; i++) {
    sp = level[i];
    for (j=0; j< bins; j++) sp[j] = ' ';
    sp[j] = 0;
  }

  if (showlevel && level_val >= min && level_val <= max) {
    for (j=0; j<bins; j++) {
      level[ (int)(NLEV*(level_val-min)/(max-min)) ][j] = '=';
    }
  }
 
  for (j=0; j<bins; j++) {
    dave = dsig = 0.; dmin = 1e20; dmax = -1e20;
    k = 0;
    for (i=j*nv/bins; i<(j+1)*nv/bins; i++) {
      dave += ac[i];
      dsig += sqr(ac[i]);
      dmin = (dmin < ac[i]) ? dmin : ac[i];
      dmax = (dmax > ac[i]) ? dmax : ac[i];
      k++;
    }
    dave /= k;
    dsig = sqrt(dsig/k - sqr(dave));
    level[ (int)(NLEV*(dmin-min)/(max-min)) ][j] = '.';
    level[ (int)(NLEV*(dmax-min)/(max-min)) ][j] = '.';
    if (dave+dsig <= max) 
      level[ (int)(NLEV*(dave+dsig-min)/(max-min)) ][j] = '-';
    if (dave-dsig >= min)
      level[ (int)(NLEV*(dave-dsig-min)/(max-min)) ][j] = '-';
    level[ (int)(NLEV*(dave-min)/(max-min)) ][j] = '*';
  }
 
  printf("\nblocked history with %d blocks, size %d ...\n",bins,nv/bins);
 
  for (j=NLEV; j>=0; j--) printf("%s\n",level[j]);
 
  return(1);
}

/********************************************************************
 * Assume here periodicity so that d[2*nd - 1 - i] = d[i]
 */

int
ccorr_matrix(double *dat[],int n,int length,int jack,
	     double *wd,int weight,double wnorm)
{
  int i,j,k;
  double *a,*w,*eval,*evec;

  a = dblarr(n);
  w = dblarr(n*n);
  eval = dblarr(n);
  evec = dblarr(n*n);

  if (!weight) wnorm = 1.0/length; else wnorm /= length;

  for (j=0; j<n; j++) {
    if (weight) for (a[j]=i=0; i<length; i++) a[j] += dat[j][i]*wd[i];
    else for (a[j]=i=0; i<length; i++) a[j] += dat[j][i];
    a[j] *= wnorm;
  }

  for (j=0; j<n; j++) for (k=0; k<=j; k++) {
    w[j+k*n] = 0;
    if (!weight) for (i=0; i<length; i++) {
      w[j+k*n] += (dat[j][i] - a[j])*(dat[k][i] - a[k]);
    } else for (i=0; i<length; i++) {
      w[j+k*n] += (dat[j][i] - a[j])*(dat[k][i] - a[k])*wd[i];
    }
    w[k+j*n] = (w[j+k*n] *= wnorm);
  }

  i = jacobi(n,w,eval,evec,1);

  printf("%d jacobi rotations\n",i);
  
  for (i=0; i<n; i++) {
    printf("%d: eigenvalue %.10lg, with vector\n",i,eval[i]);
    for (j=0; j<n; j++) printf("%.10lg ",evec[i*n+j]);
    printf("\n");
  }
  return(0);
}


/********************************************************************
 * We want here \sum_x f(x) cos(x n 2pi/L) , where L = 2*(nd-1)
 * fill in the function periodically to 2*(nd-1)
 */

int 
ftrans(double *d,int nd)
{
  int k,i,l;
  double *t,a;

  /* if (nd % 2 == 0) fprintf(stderr,"Even vector to ft!\n");
   */

  l = 2*(nd-1);
  t = dblarr(nd);
  for (k=0; k<nd; k++) {
    t[k] = d[0];
    for (i=1; i<nd; i++) t[k] += d[i]*cos(pi2*k*i/l);
    for (i=nd; i<l; i++) t[k] += d[l-i]*cos(pi2*k*i/l);
  }
  for (k=0; k<nd; k++) d[k] = t[k];
  free(t);
}

 
 
/*********************************************************************/
 
int 
errorcalc(nconf,svect,avep,sigp, fsp, fser)
     int nconf;
     double svect[], *avep, *sigp, *fsp, *fser;

#define TINT 6
 
/*
 *
 *     this subroutine analyses the data given in
 *                svect
 *     nconf is the number of the data
 *     ave   gives the average value
 *     sig   gives the naive error
 *     fsp   integrated autocorrelation time
 *     fser  error of fsp
 *
 *     TRUE SIGMA = sig*sqrt(2*fsp)
 */
 
{
  double ave,fs,sig,ax,ay,fi;
  int i,j,nc,it,nis,ns,nm ,k;
  
  ave = sig = 0.0;
 
  for(i=0; i<nconf; i++) ave += svect[i];
  ave /= nconf;
  for(i=0; i<nconf; i++) sig += sqr(svect[i]-ave);
  sig /= nconf;
  if (nconf > 1 && (sig == 0.0 || sig < 1.0e-12*sqr(ave) || noauto)) {
    /* if (!noauto) fprintf(stderr," ** small sigma/ave: %g/%g\n",sqrt(sig),ave); */
    *avep = ave;
    *sigp = sqrt(sig/(nconf-1));
    *fsp  = *fser = 0.5;
    return(1);
  }

  if (nconf <= 1) {
    *avep = ave;
    *sigp = 0.0;
    *fsp = *fser = 0.5;
    return(1);
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
  
  return(0);
} 


/********************************************************************
 *   this routine calculates autocorrelation function
 */

int 
autocorr(d,nd,res,nres)
     int nd,nres;
     double d[], res[];
{
  double ave,fs,sig,ax,ay,fi;
  int i,j,nc,it,k;
  
  ave = sig = 0.0;
 
  for(i=0; i<nd; i++) ave += d[i];
  ave /= nd;
  for(i=0; i<nd; i++) sig += sqr(d[i]-ave);
  sig /= nd;
  
  /* time correlations */
  
  fs = 0.5;

  for (it=0; it<nres; it++) {
    ax = ay = fi = 0.0;
    
    nc=nd-it;
    for (j=0; j<nc; j++) {
      fi += d[j]*d[j+it];
      ax += d[j];
      ay += d[j+it];
    }
    
    fi = (fi-ax*ay/nc)/(nc*sig);
    fs = fs+fi*nc/nd;

    res[it] = fi;
  }

  return(0);
} 

/* calculate autocorrelation, but 
 * a) do not normalise
 * b) assume <d>=0
 */


double timecorr(double d[],int nd,double res[],int nres)
{
  double ave,fs,fi;
  int nc,it,j;
  
  ave = 0.0;
 
  // for(i=0; i<nd; i++) ave += d[i];
  // ave /= nd;
  
  /* time correlations */
  
  for (it=0; it<nres; it++) {
    nc=nd-it;
    fi = 0.0;
    for (j=0; j<nc; j++) {
      fi += d[j]*d[j+it];
    }
    fi /= nc;

    if (it > 0)
      fs += fi;
    else
      fs = 0.5*fi;

    res[it] = fi;
  }

  return(fs);
} 


double dissipation_coeff(double d[],int nd,double res[],int nres)
{
  double ave,fs,fi;
  int nc,it,j;
  
  
  /* time correlations */
  
  for (it=0; it<nres; it++) {
    nc=nd-it;
    fi = 0.0;
    for (j=0; j<nc; j++) {
      double diff = d[j] - d[j+it];
      fi += diff*diff;
    }
    fi /= nc;

    res[it] = fi;
  }

  return(fs);
}





/********************************************************************
 *   this routine calculates tunneling time
 */

double
tunnelcalc(int nd,double *d,double tmin,double tmax)
{
  double min,max;
  int i,iold,ismin,ismax,num;

  ismin = ismax = num = 0;

  for (i=0; i<nd; i++) {
    if (d[i] >= tmax) {
      if (!ismax) num++;
      ismax = 1; ismin = 0;
    } else if (d[i] <= tmin) {
      if (!ismin) num++;
      ismax = 0; ismin = 1;
    }
  }

  return(1.0*nd/num);
}



