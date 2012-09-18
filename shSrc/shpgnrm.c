/* 
    @(#)shpgnrm.c	1.4
    02/04/19 16:39:47
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#include <math.h>
#include <stdio.h>

#ifdef __cplusplus
 extern "C" {
#endif

void shpgnrm(int *nv,float *xt,float *yt,float *zt,float *nrmt)
 {
  int r=0;
  int g=0;
  int b=0;
  int r0=0;
  int g0=0;
  int b0=0;
  int zero=0;
  int one=1;
  int full=255;
  int i;

/*  shade and z-buffer a polygon defined by the vertices xt,yt,zt. */

  float *x,*y,*z;
  float nx,ny,nz;
  int  *ict,*ic;
  float *nrm;
  int nn=0;
  int np=0;

  x=(float*)NULL;
  y=(float*)NULL;
  z=(float*)NULL;
  nrm=(float*)NULL;
  ic=(int*)NULL;

  ict=(int*)malloc((*nv)*sizeof(int));

  for(i=0;i<*nv;i++)ict[i]=1;
  shclippg(nv,xt,yt,zt,nrmt,ict,&nn,&x,&y,&z,&nrm,&ic);
  shdopg(&nn,x,y,z,nrm);

#ifdef FGHJKL
  x=(float*)realloc((void*)x,(nn+1)*sizeof(float));
  y=(float*)realloc((void*)y,(nn+1)*sizeof(float));
  z=(float*)realloc((void*)z,(nn+1)*sizeof(float));
  ic=(int*)realloc((void*)ic,(nn+1)*sizeof(int));

  x[nn]=x[0];
  y[nn]=y[0];
  z[nn]=z[0];
  ic[nn]=ic[0];

  shqnpln(&np);
  shnpln(&zero);
  for(i=0;i<nn;i++)
   {
    if(0||ic[i+1]!=1)
     {
/*    shqlinc(&r,&g,&b);
      shlinc(&full,&full,&full);*/
      if(fabs(x[i]-x[i+1])>1.e-6||fabs(y[i]-y[i+1])>1.e-6||fabs(z[i]-z[i+1])>1.e-6)
       {
        shlnonrm(&(x[i]),&(y[i]),&(z[i]),&(nrm[3*i]),&(x[i+1]),&(y[i+1]),&(z[i+1]),&(nrm[3*i+3]),&zero,&zero,&one);
/*      shlinc(&r,&g,&b);*/
       }
     }else{
      if(fabs(x[i]-x[i+1])>1.e-6||fabs(y[i]-y[i+1])>1.e-6||fabs(z[i]-z[i+1])>1.e-6)
       {
        shqlinc(&r0,&g0,&b0);
        shqpec(&r,&g,&b);
        if(r0!=r&&g0!=g&&b0!=b)
         {
          shlinc(&r,&g,&b);
          shlnonrm(&(x[i]),&(y[i]),&(z[i]),&(nrm[3*i]),&(x[i+1]),&(y[i+1]),&(z[i+1]),&(nrm[3*i+3]),&zero,&zero,&one);
          shlinc(&r0,&g0,&b0);
         }
       }
     }
   }
  shnpln(&np);
#endif

  free(x);
  free(y);
  free(z);
  free(nrm);
  free(ict);
  free(ic);

  return;
 }

#ifdef __cplusplus
}
#endif
