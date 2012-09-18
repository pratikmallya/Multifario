/* 
    @(#)shline2s.c	1.2
    02/04/19 16:38:38
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#include <math.h>
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#ifdef __cplusplus
 extern "C" {
#endif

void shline2s(int *i0,int *j0,int *i1,int *j1,int *ixoff,int *iyoff)
 {
  float zz;
  float xc,yc;
  float dx,dy;
  int ipix,jpix;

/*   draw a line segment */

  if((*i0)==(*i1)&&(*j0)==(*j1))
   {
    if((*i0)>0&&(*i0)<shIMax&&(*j0)>0&&(*j0)<shJMax)
     {
      zz=0.;
      shsetp(i0,j0,&sh_rl,&sh_gl,&sh_bl,&zz);
     }
    return;
   }

  xc=(*i0);
  yc=(*j0);
  dx=(*i1)-xc;
  dy=(*j1)-yc;
  zz=0.;
  while(fabs(dx)>.5||fabs(dy)>.5)
   {
    ipix=round(xc+.5+(*ixoff));
    jpix=round(yc+.5+(*iyoff));
    if(ipix>0&&ipix<shIMax&&jpix>0&&jpix<shJMax)shsetp(&ipix,&jpix,&sh_rl,&sh_gl,&sh_bl,&zz);

    dx=(*i1)-xc;
    dy=(*j1)-yc;
    if(fabs(dx)>.5||fabs(dy)>.5)
     {
      if(fabs(dx)>fabs(dy))
       {
        if(dx>0)
         {
          xc=xc+1;
          yc=yc+dy/dx;
         }else{
          xc=xc-1;
          yc=yc-dy/dx;
         }
       }else{
        if(dy>0)
         {
          yc=yc+1;
          xc=xc+dx/dy;
         }else{
          yc=yc-1;
          xc=xc-dx/dy;
         }
       }
     }
    dx=(*i1)-xc;
    dy=(*j1)-yc;
   }

  return;
 }

#ifdef __cplusplus
}
#endif
