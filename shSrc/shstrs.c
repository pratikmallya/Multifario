/* 
    @(#)shstrs.c	1.2
    02/04/19 16:41:14
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>

#ifdef __cplusplus
 extern "C" {
#endif

extern int SHnstr;
extern char **SHstr;
extern int *SHxstr;
extern int *SHystr;
extern double *SHzstr;

void shstrs(int xx,int yy,char *str)
 {
  SHstr=(char**)realloc(SHstr,(SHnstr+1)*sizeof(char*));
  SHxstr=(int*)realloc(SHxstr,(SHnstr+1)*sizeof(int));
  SHystr=(int*)realloc(SHystr,(SHnstr+1)*sizeof(int));
  SHzstr=(double*)realloc(SHzstr,(SHnstr+1)*sizeof(double));

  SHxstr[SHnstr]=xx;
  SHystr[SHnstr]=yy;
  SHzstr[SHnstr]=-1.;
  SHstr[SHnstr]=(char*)malloc((strlen(str)+1)*sizeof(char));
  strcpy(SHstr[SHnstr],str);
  SHnstr++;

  return;
 }

#ifdef __cplusplus
}
#endif
