/* 
    @(#)shname.c	1.5
    02/04/19 16:39:07
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

extern char *shOutputName;

#include <sh.h>
#include <shInternal.h>
#include <string.h>
#include <stdlib.h>

#ifdef __cplusplus
 extern "C" {
#endif

void shSetOutputFilename(char *name)
 {
  shOutputName=(char*)realloc(shOutputName,(strlen(name)+10)*sizeof(char));
  strcpy(shOutputName,name);
  if(!strcmp(shOutputFormat,"tiff"))
    strcat(shOutputName,".tif");
   else
    strcat(shOutputName,".ps");

  return;
 }

#ifdef __cplusplus
}
#endif
