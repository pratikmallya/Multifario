/* 
    @(#)shpntc.c	1.2
    02/04/19 16:40:19
   
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

int sh_rp=0;
int sh_gp=0;
int sh_bp=0;

void shpntc(int *r,int *g,int *b)
 {
  sh_rp=*r;
  sh_gp=*g;
  sh_bp=*b;
  return;
 }

void shqpntc(int *r,int *g,int *b)
 {
  *r=sh_rp;
  *g=sh_gp;
  *b=sh_bp;
  return;
 }

#ifdef __cplusplus
}
#endif
