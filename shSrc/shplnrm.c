/* 
    @(#)shplnrm.c	1.2
    02/04/19 16:40:09

    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#include <math.h>

#ifdef __cplusplus
 extern "C" {
#endif

void shplnrm(int *nv,float *xt,float *yt,float *zt)
 {
  float *n;
  float an;
  int m,m1,m2,m3;
  int i;
  int zero=0;

/*   Shade and render a polyline with normals */

  n=(float*)malloc(3*(*nv)*sizeof(float));

  for(m=0;m<*nv;m++)
   {
    m1=m;
    m3=m-1;
    if(m3<1)m3=*nv;
    m2=m+1;
    if(m2>*nv)m2=1;
    m1=1;
    m2=2;
    m3=3;
    n[3*m]=(yt[m2]-yt[m1])*(zt[m3]-zt[m1])-(yt[m3]-yt[m1])*(zt[m2]-zt[m1]);
    n[1+3*m]=(xt[m3]-xt[m1])*(zt[m2]-zt[m1])-(xt[m2]-xt[m1])*(zt[m3]-zt[m1]);
    n[2+3*m]=(xt[m2]-xt[m1])*(yt[m3]-yt[m1])-(xt[m3]-xt[m1])*(yt[m2]-yt[m1]);
    an=sqrt(n[3*m]*n[3*m]+n[1+3*m]*n[1+3*m]+n[2+3*m]*n[2+3*m]);
    if(an==0.)return;
    an=1./an;
    n[3*m]=n[3*m]*an;
    n[1+3*m]=n[1+3*m]*an;
    n[2+3*m]=n[2+3*m]*an;
   }

  for(i=0;i<*nv-1;i++)
    shlnonrm(xt+i,yt+i,zt+i,&(n[3*i]),xt+i+1,yt+i+1,zt+i+1,&(n[3*i+3]),&zero,&zero,&zero);
  return;
 }

#ifdef __cplusplus
}
#endif
