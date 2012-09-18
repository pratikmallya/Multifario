/* 
    @(#)shres.c	1.2
    02/04/19 16:40:42
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

extern int shIMax;
extern int shJMax;

#include <sh.h>

#ifdef __cplusplus
 extern "C" {
#endif

void shSetOutputResolution(int x,int y)
 {
  shIMax=x;
  shJMax=y;
  return;
 }

#ifdef __cplusplus
}
#endif
