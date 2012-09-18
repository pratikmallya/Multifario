/* 
    @(#)shdraw.c	1.2
    02/04/19 16:38:09
   
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

void shdraw(float *x,float *y,float *z,int *ip)
 {

/*   Calcomp type interface for line drawing routine. */

  if(*ip==2)
     shline(&sh_x0,&sh_y0,&sh_z0,x,y,z);
  sh_x0=*x;
  sh_y0=*y;
  sh_z0=*z;
  return;
 }

#ifdef __cplusplus
}
#endif
