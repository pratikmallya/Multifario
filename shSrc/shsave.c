/* 
    @(#)shsave.c	1.7
    03/02/20 12:21:03
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#include <multifarioConfig.h>

#ifdef __cplusplus
 extern "C" {
#endif

char *shOutputFormat=(char*)NULL;

void shsave(char *name,int *ln,int nln)
 {
/*  Write image to disk.*/

  if(!strcmp(shOutputFormat,"postscript") || !strcmp(shOutputFormat,"ps") )
    shputps(name,ln,&shIMax,&shJMax,shRedBuffer,shGreenBuffer,shBlueBuffer,nln);
#ifdef HAVE_LIBTIFF
   else if(!strcmp(shOutputFormat,"tiff"))
    shputtiff(name, shIMax,shJMax, shRedBuffer,shGreenBuffer,shBlueBuffer);
#endif
   else
    shputps(name,ln,&shIMax,&shJMax,shRedBuffer,shGreenBuffer,shBlueBuffer,nln);
  return;
 }

#ifdef __cplusplus
}
#endif
