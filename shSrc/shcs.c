/* 
    @(#)shcs.c	1.2
    02/04/19 16:37:55
   
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

void shcs(int *ifnt,int *isx,int *isy)
 {
  sh_cFont=*ifnt;
  shadow_idx=*isx;
  shadow_idy=*isy;

  return;
 }

#ifdef __cplusplus
}
#endif
