/* 
    @(#)shsrfp.c	1.2
    02/04/19 16:41:05
   
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

void shsrfp(float *a,float *d,float *s,int *n)
 {

/*    This Routine Sets the surface properties used for the Shading */

/*          a ambient intensity.                                    */
/*          d difuse intensity.                                     */
/*          s specular intensity.                                   */
/*          n specular exponent.                                    */

  sh_am=*a;
  sh_ad=*d;
  sh_as=*s;
  sh_nd=*n;

  return;
 }

#ifdef __cplusplus
}
#endif
