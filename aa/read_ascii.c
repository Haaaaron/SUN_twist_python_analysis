/* Ascii input package for labelled (or not) ascii data,
 *mainly for aa-prg.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
/* #include <float.h> */
#include <memory.h>

#include "stuff.h"

/* This routine inputs ascii data files, with two formats: 
 * 1: if label = NULL, input file is std ascii file with 
 * one record/line, and fields are separated by spaces/tabs.
 * 2: if label = "STRING", we inspect labelled measure files.  
 * Now one record consists of lines beginning with STRING.
 */

#define LINESIZE 100000
static char line[LINESIZE];
static int ndata = 0, nrows = 1;
static int nd_row0=0;

/* is this line labelled line? */
char * islabel(char *cp, char *label)
{
  while (*cp == ' ' || *cp == '\t') cp++;
  while (*cp == *label && *label != 0) { cp++; label++; }
  if (*label != 0) return( NULL );
  while (*cp != 0 && *cp != ' ' && *cp != '\t') cp++;
  while (*cp == ' ' || *cp == '\t') cp++;
  return (cp);
}

/* is this line labelled line? */
char * islegend(char *cp, char *label)
{
  static char legendstr[1000];
  char *p;

  while (*cp == ' ' || *cp == '\t') cp++;
  if (strncmp(cp,"Legend",6) == 0 || strncmp(cp,"legend",6) == 0) {
    cp += 6;
    while (*cp == ' ' || *cp == '\t') cp++;
    if ((cp = islabel(cp,label)) != NULL) {
      // found label, copy and return
      p = legendstr;
      while (*cp != 0 && *cp != '\n') {
        *p = *cp;
        p++;
        cp++;
      }
      return (legendstr);
    }
  }
  return (NULL);
}



void checkwarn(char *line)
{
  static int done = 0;

  if (done) return;
  if (strstr(line,"WARNING") != NULL) {
    fprintf(stderr,"* File contains warning(s), first:\n");
    fprintf(stderr,"%s\n",line);
    done = 1;
  }
}

int is_measure_start(char *line,int iskip)
{
  char *p;
  int it;
  p = strstr(line,"Measure_start");
  if (p != NULL) {
    p += 13;
    if (sscanf(p,"%d",&it) == 1) {
      if (it % iskip == 0) return(1);
    }
  }
  return(0);
}


int read_ascii_info( FILE *f,      // input file
		     char *label,   // label string - if NULL, plain line-ascii
		     int *rows,     // how many rows for label
                     char **legend,   // legend string
                     int iskip     // skip over (take mod of index)
		     )
{
  double d;
  char *tp,*cp;
  int i;

  *legend = NULL;

  ndata = 0;
  *rows = 1;
  if (label == NULL) {
    if (fgets(line,LINESIZE,f) == NULL) {
      fprintf(stderr,"Error - no data in file?\n");
      exit(0);
    }
    cp = line;
    do {
      d = strtod(cp,&tp);
      if (tp == cp) break;
      ndata++;
      if (*tp == ',') tp++; /* allow ,-separated numbers */
      cp = tp;
    } while (1);
    rewind (f);
    return ndata;
  }

  /* now it is labelled file */
  do {
    tp = fgets(line,LINESIZE,f);
    if (tp != NULL) checkwarn(line);
  } while (tp != NULL && !is_measure_start(line,iskip));
  if (tp == NULL) {
    fprintf(stderr,"Error - no 'Measure_start' in file?\n");
    if (iskip>1) fprintf(stderr,"with index %% %d == 0\n",iskip);
    exit(0);
  }

  nrows = ndata = 0;
  do {
    tp=fgets(line,LINESIZE,f);
    if (tp == NULL) {
      fprintf(stderr,"Error - no 'Measure_end' in file?\n");
      exit(0);
    }
    checkwarn(line);

    if (is_measure_start(line,1)) {
      // fprintf(stderr,"Warning - loose Measure_start? ignoring\n");
      nrows = ndata = 0;
      if (!is_measure_start(line,iskip)) {
        fprintf(stderr,"Confused - Measure_start index mismatch at the beginning\n");
        exit(0);
      }
    } else if (strstr(line,"Measure_end") != NULL && ndata != 0) {
      rewind(f);
      *rows = nrows;
      return ndata;
    } else {
    
      if (*legend == NULL) 
        *legend = islegend(line,label);

      cp = islabel(line,label);

      // fprintf(stderr,"Label %s  compare %s \n",label,line);

      if (cp != NULL) {
        nrows++;
        do {
          d = strtod(cp,&tp);
          if (tp == cp) break;
          ndata++;
          cp = tp;
        } while (1);
        if (nrows == 1) nd_row0 = ndata;
        else if (nrows == 2) {
          fprintf(stderr,"More than 1 row - assuming correlation length\n");
        } 
        if (ndata != nd_row0 * nrows) {
          fprintf(stderr,"Error - label %s lines not of equal length\n",label);
          exit(0);
        }
      }
    }
  } while(1);
}


int read_ascii( FILE *f,       // input file
		char *label,   // label string - if NULL, plain line-ascii
		double *dat,   // data array
                int iskip
		)
{
  int n = 0, i, c, row;
  static int res = 0;
  double d;
  char *tp,*cp;
  
  
  if (label == NULL) {
    // Read lines until it's not comment
    do {
      if (fgets(line,LINESIZE,f) == NULL) return(0);
      for (tp=line; *tp==' ' || *tp == '\t'; tp++) ;
    } while (*tp == '#');

    cp = line;
    for (i=0; i<ndata; i++) {
      dat[i] = strtod(cp,&tp);
      if (tp == cp) return(0);
      if (*tp == ',') tp++;
      cp = tp;
    }
    return(1);
  }

  /* now it is labelled file - scan it */
  res = 0;

 reallyreset:   // if we find measure_start with wrong index go here
  do {
    tp = fgets(line,LINESIZE,f);
    if (tp != NULL) checkwarn(line);
  } while (tp != NULL && !is_measure_start(line,iskip));
  if (tp == NULL) return(0);
  
 reset:   // reset everything if we find Measure_start
  if (res) {
    fprintf(stderr,"Warning: loose Measure_start? ignoring\n");
    if (!is_measure_start(line,iskip)) goto reallyreset;
  }


again:
  res = 1;

  i = 0;
  for (row=0; row<nrows; row++) {
    do {
      tp=fgets(line,LINESIZE,f);
      if (tp == NULL) return(0);  /* Only valid exit through "Measure_end"*/
      checkwarn(line);
      if (strstr(line,"Measure_start") != NULL) goto again;
      cp = islabel(line,label);
    } while (cp == NULL);

    /* handle found line */
    for (c=0; c<nd_row0; c++) {
      dat[row + nrows*c] = strtod(cp,&tp);
      if (tp == cp) {
	fprintf(stderr,"Error: expecting %d items on each line\n",nd_row0);
	exit(0);
      }
      cp = tp;
    }
  }
  do {
    tp = fgets(line,LINESIZE,f);
    if (tp != NULL && strstr(line,"Measure_start") != NULL) goto reset;
    if (tp != NULL) checkwarn(line);
  } while (tp != NULL && strstr(line,"Measure_end") == NULL);
  if (tp == NULL) return(0);
  return(1);
}



