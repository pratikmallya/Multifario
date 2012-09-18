/* 
    @(#)shqnpln.c	1.2
    02/04/19 16:40:33
   
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

void shqnpln(int *n)
 {

/* Return the number of clipping planes. */

  *n=sh_nplns;
  return;
 }

#ifdef __cplusplus
}
#endif
