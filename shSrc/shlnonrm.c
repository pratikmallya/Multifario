/* 
    @(#)shlnonrm.c	1.2
    99/10/13 11:07:16
   
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

static int shlw=0;

void shSetLineWidth(int w){shlw=w;return;}
void shtmpsetp(int,int,int,int,int,float,int);

void shlnonrm(float *x0,float *y0,float *z0,float *n0,float *x1,float *y1,float *z1,float *n1,int *ixoff,int *iyoff,int *idoff)
 {
  int rr=0;
  int gg=0;
  int bb=0;
  int i0,j0,i1,j1;
  int imask=0;

  float nn[3]={0.,0.,0.};
  float zbuf=0.;
  float dx,dy;
  int ipix,jpix;

  float s0=0.;
  float t0=0.;
  float d0;
  float s1=0.;
  float t1=0.;
  float d1;
  float xx,yy,zz;
  float xxx=0.;
  float yyy=0.;
  float zzz=0.;
  float xc,yc;
  int clipped;
  int ip;
  float direc;

/*   Project into perspective coordinates. */

  shpers(x0,y0,z0,&s0,&t0,&d0);
  shpers(x1,y1,z1,&s1,&t1,&d1);

  if((*idoff)<0)
   {
    d0=1.005*d0;
    d1=1.005*d1;
   }
  if((*idoff)>0)
   {
    d0=.995*d0;
    d1=.995*d1;
   }

/*   Transform to pixel coordinates. */

  i0=round(s0*shMax);
  j0=round(t0*shMax);

  i1=round(s1*shMax);
  j1=round(t1*shMax);

  if(i0==i1&&j0==j1)
   {
    if(i0>0&&i0<shIMax&&j0>0&&j0<shJMax)
     {
      shgetz(&i0,&j0,&zbuf);
      if(zbuf>d0)
       {
        if(shMask)
         {
          sh3dmask(x0,y0,z0,n0,&imask);
          if(imask==1)
           {
            shcolor(x0,y0,z0,n0,&sh_rl,&sh_gl,&sh_bl,&sh_rl,&sh_gl,&sh_bl,&rr,&gg,&bb);
            if(shlw==0)shsetp(&i0,&j0,&rr,&gg,&bb,&d0);
             else shtmpsetp(i0,j0,rr,gg,bb,d0,shlw);
           }
         }else{
          shcolor(x0,y0,z0,n0,&sh_rl,&sh_gl,&sh_bl,&sh_rl,&sh_gl,&sh_bl,&rr,&gg,&bb);
          if(shlw==0)shsetp(&i0,&j0,&rr,&gg,&bb,&d0);
           else shtmpsetp(i0,j0,rr,gg,bb,d0,shlw);
         }
       }
     }
   return;
  }

/*  Trace and shade the line. */

  xc=i0;
  yc=j0;
  dx=i1-xc;
  dy=j1-yc;
  while(fabs(dx)>.5||fabs(dy)>.5)
   {
    if(fabs(i1-i0)>fabs(j1-j0))
     {
      xx=s0+(float)(xc-i0)*(s1-s0)/(float)(i1-i0);
      yy=t0+(float)(xc-i0)*(t1-t0)/(float)(i1-i0);
      zz=d0+(float)(xc-i0)*(d1-d0)/(float)(i1-i0);
      nn[0]=n0[0]+(float)(xc-i0)*(n1[0]-n0[0])/(float)(i1-i0);
      nn[1]=n0[1]+(float)(xc-i0)*(n1[1]-n0[1])/(float)(i1-i0);
      nn[2]=n0[2]+(float)(xc-i0)*(n1[2]-n0[2])/(float)(i1-i0);
     }else{
      xx=s0+(float)(yc-j0)*(s1-s0)/(float)(j1-j0);
      yy=t0+(float)(yc-j0)*(t1-t0)/(float)(j1-j0);
      zz=d0+(float)(yc-j0)*(d1-d0)/(float)(j1-j0);
      nn[0]=n0[0]+(float)(yc-j0)*(n1[0]-n0[0])/(float)(j1-j0);
      nn[1]=n0[1]+(float)(yc-j0)*(n1[1]-n0[1])/(float)(j1-j0);
      nn[2]=n0[2]+(float)(yc-j0)*(n1[2]-n0[2])/(float)(j1-j0);
     }

/*    clip */

    clipped=0;
    shunpers(&xx,&yy,&zz,&xxx,&yyy,&zzz);
    if(sh_nplns>0)
     {
      for(ip=0;ip<sh_nplns;ip++)
       {
        direc=sh_plnn[  3*ip]*(xxx-sh_plno[  3*ip])+sh_plnn[1+3*ip]*(yyy-sh_plno[1+3*ip])+sh_plnn[2+3*ip]*(zzz-sh_plno[2+3*ip]);
        if(sh_oper[ip]==0)
          clipped=clipped||(sh_ipln[ip]*direc<0);
         else
          clipped=clipped&&(sh_ipln[ip]*direc<0);
       }
     }
    if(!clipped)
     {
      ipix=round(xc);
      jpix=round(yc);
      if(ipix>0&&ipix<shIMax&&jpix>0&&jpix<shJMax)
       {
        shgetz(&ipix,&jpix,&zbuf);
        if(zbuf>zz-.005)
         {
          if(shMask)
           {
            sh3dmask(&xxx,&yyy,&zzz,nn,&imask);
            if(imask==1)
             {
              shcolor(&xxx,&yyy,&zzz,nn,&sh_rl,&sh_gl,&sh_bl,&sh_rl,&sh_gl,&sh_bl,&rr,&gg,&bb);
              if(shlw==0)shsetp(&ipix,&jpix,&rr,&gg,&bb,&zz);
               else shtmpsetp(ipix,jpix,rr,gg,bb,zz,shlw);
             }
           }else{
            shcolor(&xxx,&yyy,&zzz,nn,&sh_rl,&sh_gl,&sh_bl,&sh_rl,&sh_gl,&sh_bl,&rr,&gg,&bb);
            if(shlw==0)shsetp(&ipix,&jpix,&rr,&gg,&bb,&zz);
             else shtmpsetp(ipix,jpix,rr,gg,bb,zz,shlw);
           }
         }
       }
     }

    dx=i1-xc;
    dy=j1-yc;
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
   }

  return;
 }

void shtmpsetp(int i0,int j0,int r,int g,int b,float z,int w)
 {
  int iii,jjj,ii,jj;
  float zbuf;

  for(iii=-w;iii<w+1;iii++)
   {
    for(jjj=-w;jjj<w+1;jjj++)
     {
      ii=i0+iii;jj=j0+jjj;
      shgetz(&ii,&jj,&zbuf);
      if(sqrt(1.*(iii*iii+jjj*jjj)/w/w)<=1)
       {
        if(i0>0&&i0<shIMax&&j0>0&&j0<shJMax)
         {
          if(zbuf>z-.005)shsetp(&ii,&jj,&r,&g,&b,&z);
         }
       }
     }
   }

  return;
 }

#ifdef __cplusplus
}
#endif
