/* 
    @(#)shpl.c	1.2
    02/04/19 16:39:53
   
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

void shpl(int *n,float *x,float *y,float *z)
 {
  int i;

/*  Shade and render a polyline */

  for(i=0;i<*n-1;i++)
    shline(x+i,y+i,z+i,x+i+1,y+i+1,z+i+1);

  return;
 }

#ifdef __cplusplus
}
#endif
