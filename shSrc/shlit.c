/* 
    @(#)shlit.c	1.3
    02/04/19 16:38:47
   
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

int sh_nlit=0;
int sh_mlit=0;
int *sh_rs=(int*)NULL;
int *sh_gs=(int*)NULL;
int *sh_bs=(int*)NULL;
int *sh_type=(int*)NULL;
double *sh_lit=(double*)NULL;

void shlit(int *m,float *x,float *y,float *z,int *r,int *g,int *b)
 {
  float xt=0.;
  float yt=0.;
  float zt=0.;

/*     This Routine Sets the parameters associated with light source m. */

/*              x,y,z   (real) position.                                */
/*              r,g,b   (integer 0-255) color.                          */
/*              type    (integer 0-1) 0=point,1=direction.              */

  if((*m)>=sh_mlit)
   {
    if(sh_mlit==0)
     {
      sh_mlit=(*m)+10;
      sh_rs=(int*)malloc(sh_mlit*sizeof(int));
      sh_gs=(int*)malloc(sh_mlit*sizeof(int));
      sh_bs=(int*)malloc(sh_mlit*sizeof(int));
      sh_type=(int*)malloc(sh_mlit*sizeof(int));
      sh_lit=(double*)malloc(3*sh_mlit*sizeof(double));
     }else{
      sh_mlit=(*m)+10;
      sh_rs=(int*)realloc((void*)sh_rs,sh_mlit*sizeof(int));
      sh_gs=(int*)realloc((void*)sh_gs,sh_mlit*sizeof(int));
      sh_bs=(int*)realloc((void*)sh_bs,sh_mlit*sizeof(int));
      sh_type=(int*)realloc((void*)sh_type,sh_mlit*sizeof(int));
      sh_lit=(double*)realloc((void*)sh_lit,sh_mlit*sizeof(double));
     }
   }
  sh_rs[(*m)-1]=*r;
  sh_gs[(*m)-1]=*g;
  sh_bs[(*m)-1]=*b;

  shpers(x,y,z,&xt,&yt,&zt);
  sh_lit[0+3*((*m)-1)]=xt;
  sh_lit[1+3*((*m)-1)]=yt;
  sh_lit[2+3*((*m)-1)]=zt;
  sh_type[(*m)-1]=1;

  return;
 }

#ifdef __cplusplus
}
#endif
