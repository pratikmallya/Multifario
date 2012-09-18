/* 
    @(#)shclr.c	1.2
    02/04/19 16:37:47
   
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

void shclr()
 {
  int i,j;
  float two=2.;

  for(i=0;i<shIMax;i++)
   {
    for(j=0;j<shJMax;j++)
     {
      shsetp(&i,&j,&shRedBackground,&shGreenBackground,&shBlueBackground,&two);
     }
   }

  return;
 }

#ifdef __cplusplus
}
#endif
