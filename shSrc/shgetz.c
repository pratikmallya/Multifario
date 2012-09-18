/* 
    @(#)shgetz.c	1.2
    02/04/19 16:38:20
   
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

void shgetz(int *ipix,int *jpix, float *zz)
 {
  int index;
 
/*   Return the z-buffer. */

  *zz=0.;

  if((*ipix<0) || (*jpix<0) )return;
  if((*ipix>shIMax-1) || (*jpix>shJMax-1) )return;

  index=(*ipix)+shIMax*(*jpix);

  *zz=shZBuffer[index];

  return;
 }

#ifdef __cplusplus
}
#endif
