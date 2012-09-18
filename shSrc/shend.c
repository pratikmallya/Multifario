/* 
    @(#)shend.c	1.4
    02/04/19 16:38:15
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
extern char *shOutputName;
void shfreeStrings(void);

#ifdef __cplusplus
 extern "C" {
#endif

void shend()
 {
  int i;

  if(shOutputName!=(char*)NULL)free(shOutputName);
  shfreeStrings();

  if(shRedBuffer!=(unsigned char*)NULL)free(shRedBuffer);
  if(shGreenBuffer!=(unsigned char*)NULL)free(shGreenBuffer);
  if(shBlueBuffer!=(unsigned char*)NULL)free(shBlueBuffer);
  if(shZBuffer!=(float*)NULL)free(shZBuffer);

  if(sh_kfont!=(int*)NULL)free(sh_kfont);
  if(sh_dir!=(int*)NULL)free(sh_dir);
  if(sh_fmag!=(int*)NULL)free(sh_fmag);
  if(sh_image!=(char*)NULL)free(sh_image);
  for(i=0;i<sh_nFonts;i++)
    if(sh_fontnm[i]!=(char*)NULL)free(sh_fontnm[i]);
  if(sh_fontnm!=(char**)NULL)free(sh_fontnm);

  if(sh_rs!=(int*)NULL)free(sh_rs);
  if(sh_gs!=(int*)NULL)free(sh_gs);
  if(sh_bs!=(int*)NULL)free(sh_bs);
  if(sh_type!=(int*)NULL)free(sh_type);
  if(sh_lit!=(double*)NULL)free(sh_lit);

  if(sh_plno!=(float*)NULL)free(sh_plno);
  if(sh_plnn!=(float*)NULL)free(sh_plnn);
  if(sh_ipln!=(int*)NULL)free(sh_ipln);
  if(sh_oper!=(int*)NULL)free(sh_oper);

  return;
 }

#ifdef __cplusplus
}
#endif
