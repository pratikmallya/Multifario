/* 
    @(#)shmask.c	1.2
    02/04/19 16:39:02
   
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

int shMask=0;

void shmask(int *swtch)
 {
  shMask=*swtch;
  return;
 }

#ifdef __cplusplus
}
#endif
