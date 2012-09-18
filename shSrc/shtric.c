/* 
    @(#)shtric.c	1.2
    02/04/19 16:41:28
   
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

void shtric(int *r,int *g,int *b)
 {
  shpgc(r,g,b);
  shpec(r,g,b);
  shbgc(r,g,b);
  shbec(r,g,b);
  return;
 }

#ifdef __cplusplus
}
#endif
