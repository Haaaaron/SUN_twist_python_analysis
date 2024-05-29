#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
/* #include <float.h> */
#include <memory.h>

#include "stuff.h"
#ifdef cray
double asinh(double x) { return(log(sqrt(x*x+1.0) + x)); }
double acosh(double x) { return(log(sqrt(x*x-1.0) + x)); }
double atanh(double x) { return(0.5*log((1.0+x)/(1.0-x))); }
#endif

#define pi 3.1415926535897929

char * getnumber(char *s,double *d);
double evallist(double d[],int dn, int prec);
void eval2(double *val1, double *val2, double d[],int dn);


static char *cmd, *cmd0;
/* This is modified by the calling program, to produce correct
 * index
 */
int calclist_index = 0;

/**********************************************************
 * this evaluates the data arithmetic string
 */

double 
calclist(double d[],int dn, char *incmd)
{
  double val;

  /* save current number of calclist calls - return with #i */

  cmd0 = cmd = incmd;
  val = evallist(d,dn,0);
  if (*cmd == 0) return(val);
  if (*cmd == ')') {
    fprintf(stderr,"Extra \')\' : %s\n",cmd0);
    exit(0);
  }
  if (*cmd == ',') {
    fprintf(stderr,"Extra \',\' : %s\n",cmd0);
    exit(0);
  }
  fprintf(stderr,"Parser error: %s\n",cmd0);
}


double
evallist(double d[],int dn,int prec)
{
  int id;
  char *sp;
  double val,val2;
  
  /* get argument first */
  while (*cmd == ' ') cmd++;
  if (*cmd == '(') {
    cmd++; 
    val = evallist(d,dn,0);
    if (*cmd != ')') {
      fprintf(stderr,"Expecting \')\': \'%s\'\n",cmd0);
      exit(0);
    }
    cmd++;
  } else if (*cmd == '#') {
    cmd++;
    if (*cmd == 'i') {
      /* print out index */
      val = calclist_index;
      cmd++;
    } else {
      if ((cmd=getnumber(cmd,&val)) == NULL) {
	fprintf(stderr,"Expecting 'i'/number after #: \'%s\'\n",cmd0);
	exit(0);
      }
      id = val;
      if (id > dn) {
	fprintf(stderr,"#%d too large: max %d\n",id,dn);
	exit(0);
      }
      val = d[id-1];
    }
  }
  else if ((sp=getnumber(cmd,&val)) != NULL) cmd = sp;
  else {
    /* now something else as ( or #; check for function name */
    if      (strncmp(cmd,"sqrt(",5) == 0) { cmd+=5; val = sqrt(evallist(d,dn,0)); }
    else if (strncmp(cmd,"abs(",4) == 0)  { cmd+=4; val = fabs(evallist(d,dn,0)); }
    else if (strncmp(cmd,"sin(",4) == 0)  { cmd+=4; val = sin(evallist(d,dn,0)); }
    else if (strncmp(cmd,"cos(",4) == 0)  { cmd+=4; val = cos(evallist(d,dn,0)); }
    else if (strncmp(cmd,"tan(",4) == 0)  { cmd+=4; val = tan(evallist(d,dn,0)); }
    else if (strncmp(cmd,"exp(",4) == 0)  { cmd+=4; val = exp(evallist(d,dn,0)); }
    else if (strncmp(cmd,"log(",4) == 0)  { cmd+=4; val = log(evallist(d,dn,0)); }
    else if (strncmp(cmd,"asin(",5) == 0) { cmd+=5; val = asin(evallist(d,dn,0)); }
    else if (strncmp(cmd,"acos(",5) == 0) { cmd+=5; val = acos(evallist(d,dn,0)); }
    else if (strncmp(cmd,"atan(",5) == 0) { cmd+=5; val = atan(evallist(d,dn,0)); }
    else if (strncmp(cmd,"sinh(",5) == 0) { cmd+=5; val = sinh(evallist(d,dn,0)); }
    else if (strncmp(cmd,"cosh(",5) == 0) { cmd+=5; val = cosh(evallist(d,dn,0)); }
    else if (strncmp(cmd,"tanh(",5) == 0) { cmd+=5; val = tanh(evallist(d,dn,0)); }
    else if (strncmp(cmd,"asinh(",6) == 0){ cmd+=6; val = asinh(evallist(d,dn,0)); }
    else if (strncmp(cmd,"acosh(",6) == 0){ cmd+=6; val = acosh(evallist(d,dn,0)); }
    else if (strncmp(cmd,"atanh(",6) == 0){ cmd+=6; val = atanh(evallist(d,dn,0)); }
    else if (strncmp(cmd,"ceil(",5) == 0) { cmd+=5; val = ceil(evallist(d,dn,0)); }
    else if (strncmp(cmd,"floor(",6) == 0){ cmd+=6; val = floor(evallist(d,dn,0)); }
    else if (strncmp(cmd,"min(",4) == 0)  { 
      cmd+=4; eval2(&val,&val2,d,dn); val = smaller(val,val2); }
    else if (strncmp(cmd,"max(",4) == 0)  {
      cmd+=4; eval2(&val,&val2,d,dn); val = greater(val,val2); }
    else if (strncmp(cmd,"atan2(",6) == 0)  {
      cmd+=6; eval2(&val,&val2,d,dn); val = atan2(val,val2); }	
    else if (strncmp(cmd,"pi",2) == 0) {
      cmd+=2; val = pi; }
    else { fprintf(stderr,"Unknown stuff: %s\n",cmd0); exit(0); }
    
    if (*cmd != ')') {
      fprintf(stderr,"Expecting \')\' after the function name: \'%s\'\n",cmd0);
      exit(0);
    }
    cmd++;
  }

  /* now the operator */
  while (*cmd == ' ') cmd++;
    
  while (*cmd) {
    switch (*cmd) {
    case ',': return(val);
    case ')': return(val);
    case '^': 
      if (prec >= 5) return(val); 
      cmd++; 
      val = pow(val,evallist(d,dn,5));
      break;
    case '*': 
      if (prec >= 4) return(val); 
      cmd++; 
      val *= evallist(d,dn,4);
      break;
    case '/': 
      if (prec >= 4) return(val); 
      cmd++; 
      val /= evallist(d,dn,4);
      break;
    case '%': 
      if (prec >= 4) return(val); 
      cmd++; 
      val = fmod(val,evallist(d,dn,4));
      break;
    case '+':
      if (prec >= 3) return(val); 
      cmd++;
      val += evallist(d,dn,3);
      break;
    case '-':
      if (prec >= 3) return(val); 
      cmd++; 
      val -= evallist(d,dn,3);
      break;

    case '=':
      if (cmd[1] != '=') {  fprintf(stderr,"Expecting '==': %s\n",cmd0); exit(0); }
      if (prec >= 2) return(val);
      cmd += 2;
      val = (val == evallist(d,dn,2));
      break;
    case '<':
      if (prec >= 2) return(val);
      cmd++;
      if (cmd[0] == '=') {
	val = (val <= evallist(d,dn,2)); cmd++;
      } else val = (val < evallist(d,dn,2));
      break;
    case '>':
      if (prec >= 2) return(val);
      cmd++;
      if (cmd[0] == '=') {
	val = (val >= evallist(d,dn,2)); cmd++;
      } else val = (val > evallist(d,dn,2));
      break;
    case '!':
      if (cmd[1] != '=') {  fprintf(stderr,"Expecting '!=': %s\n",cmd0); exit(0); }
      if (prec >= 2) return(val);
      cmd += 2;
      val = (val != evallist(d,dn,2));
      break;

    default:
      fprintf(stderr,"Unknown operator: %s\n",cmd0);
      exit(0);
    }
  }
  return(val);
}

/**************************************************************
 *    eval 2 arguments -- must be comma! 
 */

void
eval2(double *val1, double *val2, double d[],int dn)
{
  *val1 = evallist(d,dn,0);
  if (*cmd != ',') { fprintf(stderr,"Expecting ',' : %s\n",cmd0); exit(0);}
  cmd++;
  *val2 = evallist(d,dn,0);
}


/**************************************************************
 *   hop number
 */
char *
getnumber(char *s,double *d)
{
  int dot = 0;

  if (sscanf(s,"%lg",d) != 1) return(NULL);

  while (*s == ' ' || *s == '\t') s++;
  if (*s == '+' || *s == '-') s++;
  if (*s == '.') { s++; dot = 1; }
  while (*s <= '9' && *s >= '0') s++;
  if ((!dot) && *s == '.') s++;
  while (*s <= '9' && *s >= '0') s++;
  if (*s == 'e' || *s == 'E') {
    s++;
    if (*s == '+' || *s == '-') s++;
    while (*s <= '9' && *s >= '0') s++;
  }
  return(s);
}

