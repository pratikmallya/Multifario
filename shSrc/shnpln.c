/* 
    @(#)shnpln.c	1.2
    02/04/19 16:39:15
   
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

void shnpln(int *n)
 {

/* Set the number of clipping planes. */

  sh_nplns=*n;
  return;
 }

#ifdef __cplusplus
}
#endif
