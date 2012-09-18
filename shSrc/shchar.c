/* 
    @(#)shchar.c	1.3
    02/04/19 16:37:36
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#ifdef __cplusplus
 extern "C" {
#endif

void shchar(int *ix,int *iy,float *z,char *chr,int *ired,int *igrn,int *iblu,int *iw)
 {

/*        places a character string in the image */

/*        Make a list, with ix,iy,iz and the string */
/*        At the end of the postsript file, loop  */
/*          over list, if iz is visible, place postscript commands to
            draw the string in the output file */

/*        places a character in the sh_image */

  int xref=0;
  int yref=0;
  int ipwid=0;
  int iphgt=0;
  int width=0;
  int ibyte;
  int ichr;
  int ifnt;
  int idir;
  int xx,yy;
  float zz;
  int zero=0;
  int i,j;
  int jx,jy;
  int ib;
  int level;

  float conv;

/*     conv converts from sp to pixels */

  conv=240./72.27/65536.;
  ifnt=sh_kfont[sh_cFont];

/*  Blank is width of lower case m (?) */

  ichr=chr[0];
  if(chr[0]==' ')ichr='m';

  idir=sh_dir[ifnt]+4+ichr*16;
  sscanf(sh_image+idir,"%4.4d",&width);
  *iw=round(width*(sh_fmag[ifnt])*0.001*1.2*conv);

  if(chr[0]==' ')return;

  sscanf(sh_image+idir+ 4,"%4.4d",&xref);
  sscanf(sh_image+idir+ 6,"%4.4d",&yref);
  sscanf(sh_image+idir+ 8,"%4.4d",&ipwid);
  sscanf(sh_image+idir+10,"%4.4d",&iphgt);
  sscanf(sh_image+idir+12,"%4.4d",&ibyte);
  ibyte=ibyte+sh_dir[ifnt]-1;

  jx=round((*ix)-xref);
  jy=round((*iy)+yref);

/*    put shadow in frame */

  if(shadow_idx!=0||shadow_idy!=0)
   {
    ib=ibyte;
    for(j=0;j<iphgt;j++)
     {
      for(i=0;i<ipwid;i++)
       {
        level=sh_image[ib];
        if(level!=0)
         {
          zz=*z-.01;
          xx=jx+i+shadow_idx;
          yy=jy-j+shadow_idy;
          shsetp(&xx,&yy,&zero,&zero,&zero,&zz);
         }
        ib++;
       }
     }
   }

/*    put character in frame */

  ib=ibyte;
  for(j=0;j<iphgt;j++)
   {
    for(i=0;i<ipwid;i++)
     {
      level=sh_image[ib];
      if(level!=0)
       {
        zz=*z-.01;
        xx=jx+i;
        yy=jy-j;
        shsetp(&xx,&yy,ired,igrn,iblu,&zz);
       }
      ib++;
     }
   }

  return;
 }

#ifdef __cplusplus
}
#endif
