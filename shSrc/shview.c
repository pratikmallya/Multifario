/* 
    @(#)shview.c	1.2
    02/04/19 16:41:33
   
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

void shview(float *dist,float *alpha,float *beta,float *xmin,float *xmax,float *ymin,float *ymax,float *zmin,float *zmax)
 {
  float zero=0.;
  float one=0.;

/*    Set view point and scale.  */

  shsize(&one,&zero,&zero,&one,&zero);
  shvw(dist,alpha,beta,xmin,xmax,ymin,ymax,zmin,zmax);
  shscal(xmin,xmax,ymin,ymax,zmin,zmax);
  return;
 }

#ifdef __cplusplus
}
#endif
