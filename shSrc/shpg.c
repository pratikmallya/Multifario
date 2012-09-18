/* 
    @(#)shpg.c	1.2
    02/04/19 16:39:30
   
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

void shpg(int *nv,float *xt,float *yt,float *zt)
 {
  int r,g,b;
  int r0,g0,b0;
  int i;
  float an;

/*  shade and z-buffer a polygon defined by the vertices xt,yt,zt. */

  float nx,ny,nz;
  float *nrm;

  nrm=(float*)malloc(3*(*nv)*sizeof(float));

/* find normal */

  nx=(yt[1]-yt[0])*(zt[2]-zt[0])-(yt[2]-yt[0])*(zt[1]-zt[0]);
  ny=(xt[2]-xt[0])*(zt[1]-zt[0])-(xt[1]-xt[0])*(zt[2]-zt[0]);
  nz=(xt[1]-xt[0])*(yt[2]-yt[0])-(xt[2]-xt[0])*(yt[1]-yt[0]);
  an=sqrt(nx*nx+ny*ny+nz*nz);
  if(an==0.)return;
  an=1./an;
  nx=nx*an;
  ny=ny*an;
  nz=nz*an;

  for(i=0;i<*nv;i++)
   {
    nrm[  3*i]=nx;
    nrm[1+3*i]=ny;
    nrm[2+3*i]=nz;
   }

  shpgnrm(nv,xt,yt,zt,nrm);

  free(nrm);

  return;
 }

#ifdef __cplusplus
}
#endif
