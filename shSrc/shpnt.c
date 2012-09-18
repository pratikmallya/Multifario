/* 
    @(#)shpnt.c	1.3
    02/04/19 16:40:15
   
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

void shpnt(float *x,float *y,float *z)
 {

/*     Shade and Render a point. */

  float s=0.;
  float t=0.;
  float d=0.;
  int clipped;
  int ip;
  float direc;
  int ipix,jpix;
  float zbuf=0.;

/*     clip */

  clipped=0;
  if(sh_nplns>0)
   {
    for(ip==0;ip<sh_nplns;ip++)
     {
      direc=sh_plnn[  3*ip]*((*x)-sh_plno[  3*ip])+sh_plnn[1+3*ip]*((*y)-sh_plno[1+3*ip])+sh_plnn[2+3*ip]*((*z)-sh_plno[2+3*ip]);
      if(sh_oper[ip]==0)
        clipped=clipped||(sh_ipln[ip]*direc==0);
       else
        clipped=clipped&&(sh_ipln[ip]*direc==0);
     }
    if(clipped)return;
   }

  shpers(x,y,z,&s,&t,&d);
  ipix=round(s*shMax);
  jpix=round(t*shMax);
  if(ipix>0&&ipix<shIMax&&jpix>0&&jpix<shJMax)
   {
    shgetz(&ipix,&jpix,&zbuf);
    if(zbuf>d)shsetp(&ipix,&jpix,&sh_rp,&sh_gp,&sh_bp,&d);
   }

  return;
 }

#ifdef __cplusplus
}
#endif
