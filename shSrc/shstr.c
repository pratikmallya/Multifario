/* 
    @(#)shstr.c	1.5
    02/04/19 16:41:09
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#ifdef __cplusplus
 extern "C" {
#endif

int SHnstr=0;
char **SHstr=(char**)NULL;
int *SHxstr=(int*)NULL;
int *SHystr=(int*)NULL;
double *SHzstr=(double*)NULL;

void shstr(float x,float y,float z,char *str)
 {
  float xx=0.;
  float yy=0.;
  float zz=0.;

  SHstr=(char**)realloc(SHstr,(SHnstr+1)*sizeof(char*));
  SHxstr=(int*)realloc(SHxstr,(SHnstr+1)*sizeof(int));
  SHystr=(int*)realloc(SHystr,(SHnstr+1)*sizeof(int));
  SHzstr=(double*)realloc(SHzstr,(SHnstr+1)*sizeof(double));

  shpers(&x,&y,&z,&xx,&yy,&zz);

  SHxstr[SHnstr]=round(xx*shMax);
  SHystr[SHnstr]=round(yy*shMax);
  SHzstr[SHnstr]=zz;
  SHstr[SHnstr]=(char*)malloc((strlen(str)+1)*sizeof(char));
  strcpy(SHstr[SHnstr],str);
  SHnstr++;

  return;
 }

void shdoStrings(FILE *fid,float x,float y)
 {
  int i;
  int ix,iy;

  fprintf(fid,"initclip\n");
  fprintf(fid,"/Helvetica findfont %f scalefont setfont\n",20.*shIMax/1024);
  fprintf(fid,"0. setgray\n");
  for(i=0;i<SHnstr;i++)
   {
    ix=round(SHxstr[i]/x);
    iy=round(SHystr[i]/y);
    fflush(stdout);
    if((SHstr[i])[0]!='/')
      fprintf(fid,"%d %d moveto (%s) show\n",ix,iy,SHstr[i]);
     else
      fprintf(fid,"%s\n",SHstr[i]);
   }
  return;
 }

void shfreeStrings()
 {
  int i;

  for(i=0;i<SHnstr;i++)
    if(SHstr[i]!=(char*)NULL)free(SHstr[i]);
  if(SHxstr!=(int*)NULL)free(SHxstr);
  if(SHystr!=(int*)NULL)free(SHystr);
  if(SHzstr!=(double*)NULL)free(SHzstr);
  
  return;
 }

#ifdef __cplusplus
}
#endif
