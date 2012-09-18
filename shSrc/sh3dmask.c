/* 
    @(#)sh3dmask.c	1.2
    02/04/19 16:37:26
   
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

void sh3dmask(float *x,float *y,float *z,float *nrm,int *imask)
 {
  int cx,cy,cz;
  float tx,ty,tz;
  int n;

  n=5;

  tx=(*x)*n-floor((*x)*n);
  ty=(*y)*n-floor((*y)*n);
  tz=(*z)*n-floor((*z)*n);
  if(tx<0)tx=tx+1.;
  if(ty<0)ty=ty+1.;
  if(tz<0)tz=tz+1.;

  cx=fabs(tx-.5)<.4;
  cy=fabs(ty-.5)<.4;
  cz=fabs(tz-.5)<.4;

  *imask=0;
  if(!cx&&fabs(nrm[0])<.9)*imask=1;
  if(!cy&&fabs(nrm[1])<.9)*imask=1;
  if(!cz&&fabs(nrm[2])<.9)*imask=1;

  return;
 }

#ifdef __cplusplus
}
#endif
