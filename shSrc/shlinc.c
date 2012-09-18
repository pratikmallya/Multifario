/* 
    @(#)shlinc.c	1.2
    02/04/19 16:38:29
   
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

int sh_rl=0;
int sh_gl=0;
int sh_bl=0;

void shlinc(int *r,int *g,int *b)
 {
  sh_rl=*r;
  sh_gl=*g;
  sh_bl=*b;
 }

void shqlinc(int *r,int *g,int *b)
 {
  *r=sh_rl;
  *g=sh_gl;
  *b=sh_bl;
 }

#ifdef __cplusplus
}
#endif
