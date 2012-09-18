/* 
    @(#)shsetp.c	1.2
    02/04/19 16:40:56
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#include <stdio.h>

#ifdef __cplusplus
 extern "C" {
#endif

void shsetp(int *ipix,int *jpix,int *ired,int *igrn,int *iblu, float *zz)
 {
/*
    Set the pixel color and z-buffer.
*/
  int index;

/*printf("shsetp %d %d %d %d %d %f\n",*ipix,*jpix,*ired,*igrn,*iblu,*zz);fflush(stdout);
  printf("shsetp, shZBuffer=0x%8.8x\n",shZBuffer);fflush(stdout);*/
  if(*ipix<0||*jpix<0)return;
  if(*ipix>=shIMax||*jpix>=shJMax)return;

  index=*ipix+shIMax*(*jpix);
  shRedBuffer[index]=*ired;
  shGreenBuffer[index]=*igrn;
  shBlueBuffer[index]=*iblu;
  shZBuffer[index]=*zz;

  return;
 }

#ifdef __cplusplus
}
#endif
