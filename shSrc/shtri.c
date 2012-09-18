/* 
    @(#)shtri.c	1.3
    02/04/19 16:41:22
   
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

void shtri(float *xt1,float *yt1,float *zt1,float *xt2,float *yt2,float *zt2,float *xt3,float *yt3,float *zt3)
 {

/*    Shade a triangle. Note this just calls shpg.*/

  float x[3]={0.,0.,0.};
  float y[3]={0.,0.,0.};
  float z[3]={0.,0.,0.};
  int n;

  x[0]=*xt1;
  x[1]=*xt2;
  x[2]=*xt3;
  y[0]=*yt1;
  y[1]=*yt2;
  y[2]=*yt3;
  z[0]=*zt1;
  z[1]=*zt2;
  z[2]=*zt3;
  n=3;
  shpg(&n,x,y,z);

  return;
 }

#ifdef __cplusplus
}
#endif
