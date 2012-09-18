/* 
    @(#)shpause.c	1.3
    02/04/19 16:39:19
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

char *shOutputName=(char*)NULL;

void shpause()
 {
  char string[256]="";
  int n;
  char *basename;

  if(shOutputName==(char*)NULL)
   {
    sprintf(string,"%s%3.3i.ps",getenv("_"),sh_frame);
    sh_frame++;
  
    if((basename=strrchr(string,'/'))==(char*)NULL)basename=string;
     else basename++;
  
    n=strlen(basename);
    shsave(basename,&n,n);
   }else{
    n=strlen(shOutputName);
    shsave(shOutputName,&n,n);
   }
 }

#ifdef __cplusplus
}
#endif
