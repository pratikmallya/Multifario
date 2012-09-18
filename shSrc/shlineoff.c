/* 
    @(#)shlineoff.c	1.2
    02/04/19 16:38:42
   
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

void shlineoff(float *x0,float *y0,float *z0,float *x1,float *y1,float *z1,int *ixoff,int *iyoff,int *idoff)
 {

/* Line routine with offset */

  shline(x0,y0,z0,x1,y1,z1);
  return;
 }

#ifdef __cplusplus
}
#endif
