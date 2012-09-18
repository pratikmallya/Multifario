/* 
    @(#)shsphere.c	1.3
    02/04/19 16:41:00
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#include <math.h>
#include <stdio.h>

#ifdef __cplusplus
 extern "C" {
#endif

void shsphere(float *xo,float *yo,float *zo,float *r)
 {
  float *xt,*yt,*zt;
  float x[3]={0.,0.,0.};
  float y[3]={0.,0.,0.};
  float z[3]={0.,0.,0.};
  float nrm[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float xa,ya,za,xb,yb,zb,xc,yc,zc;
  float x0,y0,z0,x1,y1,z1,x2,y2,z2;
  float d;
  int maxtri;
  int i;
  int itri,ntri,mtri;
  int n,ismooth;
  int nlevel;

/* Draws a tesselated sphere, 8*4**n triangles, if ismooth^=0 smoothed. */

  ismooth=1;
  nlevel=4;

  ismooth=0;
  nlevel=3;

  maxtri=8;
  for(i=0;i<nlevel;i++)maxtri=maxtri*4;

  xt=(float*)malloc(3*maxtri*sizeof(float));
  yt=(float*)malloc(3*maxtri*sizeof(float));
  zt=(float*)malloc(3*maxtri*sizeof(float));

/* 1 */
  ntri=0;
  xt[  3*ntri]= 0.;
  yt[  3*ntri]= 1.;
  zt[  3*ntri]= 0.;
  xt[1+3*ntri]= 0.;
  yt[1+3*ntri]= 0.;
  zt[1+3*ntri]= 1.;
  xt[2+3*ntri]= 1.;
  yt[2+3*ntri]= 0.;
  zt[2+3*ntri]= 0.;
  ntri++;

/* 2 */
  xt[  3*ntri]= 0.;
  yt[  3*ntri]= 1.;
  zt[  3*ntri]= 0.;
  xt[1+3*ntri]=-1.;
  yt[1+3*ntri]= 0.;
  zt[1+3*ntri]= 0.;
  xt[2+3*ntri]= 0.;
  yt[2+3*ntri]= 0.;
  zt[2+3*ntri]= 1.;
  ntri++;

/* 3 */
  xt[  3*ntri]= 0.;
  yt[  3*ntri]= 1.;
  zt[  3*ntri]= 0.;
  xt[1+3*ntri]= 0.;
  yt[1+3*ntri]= 0.;
  zt[1+3*ntri]=-1.;
  xt[2+3*ntri]=-1.;
  yt[2+3*ntri]= 0.;
  zt[2+3*ntri]= 0.;
  ntri++;

/* 4 */
  xt[  3*ntri]= 0.;
  yt[  3*ntri]= 1.;
  zt[  3*ntri]= 0.;
  xt[1+3*ntri]= 1.;
  yt[1+3*ntri]= 0.;
  zt[1+3*ntri]= 0.;
  xt[2+3*ntri]= 0.;
  yt[2+3*ntri]= 0.;
  zt[2+3*ntri]=-1.;
  ntri++;

/* 5 */
  xt[  3*ntri]= 0.;
  yt[  3*ntri]=-1.;
  zt[  3*ntri]= 0.;
  xt[2+3*ntri]= 0.;
  yt[2+3*ntri]= 0.;
  zt[2+3*ntri]= 1.;
  xt[1+3*ntri]= 1.;
  yt[1+3*ntri]= 0.;
  zt[1+3*ntri]= 0.;
  ntri++;

/*  6 */
  xt[  3*ntri]= 0.;
  yt[  3*ntri]=-1.;
  zt[  3*ntri]= 0.;
  xt[2+3*ntri]=-1.;
  yt[2+3*ntri]= 0.;
  zt[2+3*ntri]= 0.;
  xt[1+3*ntri]= 0.;
  yt[1+3*ntri]= 0.;
  zt[1+3*ntri]= 1.;
  ntri++;

/*  7 */
  xt[  3*ntri]= 0.;
  yt[  3*ntri]=-1.;
  zt[  3*ntri]= 0.;
  xt[2+3*ntri]= 0.;
  yt[2+3*ntri]= 0.;
  zt[2+3*ntri]=-1.;
  xt[1+3*ntri]=-1.;
  yt[1+3*ntri]= 0.;
  zt[1+3*ntri]= 0.;
  ntri++;

/*  8 */
  xt[  3*ntri]= 0.;
  yt[  3*ntri]=-1.;
  zt[  3*ntri]= 0.;
  xt[2+3*ntri]= 1.;
  yt[2+3*ntri]= 0.;
  zt[2+3*ntri]= 0.;
  xt[1+3*ntri]= 0.;
  yt[1+3*ntri]= 0.;
  zt[1+3*ntri]=-1.;
  ntri++;

  if(nlevel>0)
   {
    for(i=0;i<nlevel;i++)
     {
      mtri=ntri;
      for(itri=0;itri<ntri;itri++)
       {
        x0=xt[3*itri];
        y0=yt[3*itri];
        z0=zt[3*itri];
        x1=xt[1+3*itri];
        y1=yt[1+3*itri];
        z1=zt[1+3*itri];
        x2=xt[2+3*itri];
        y2=yt[2+3*itri];
        z2=zt[2+3*itri];

        xa=.5*(x0+x1);
        ya=.5*(y0+y1);
        za=.5*(z0+z1);
        d=1./sqrt(xa*xa+ya*ya+za*za);
        xa=xa*d;
        ya=ya*d;
        za=za*d;

        xb=.5*(x1+x2);
        yb=.5*(y1+y2);
        zb=.5*(z1+z2);
        d=1./sqrt(xb*xb+yb*yb+zb*zb);
        xb=xb*d;
        yb=yb*d;
        zb=zb*d;

        xc=.5*(x2+x0);
        yc=.5*(y2+y0);
        zc=.5*(z2+z0);
        d=1./sqrt(xc*xc+yc*yc+zc*zc);
        xc=xc*d;
        yc=yc*d;
        zc=zc*d;

        xt[  3*itri]=x0;
        yt[  3*itri]=y0;
        zt[  3*itri]=z0;
        xt[1+3*itri]=xa;
        yt[1+3*itri]=ya;
        zt[1+3*itri]=za;
        xt[2+3*itri]=xc;
        yt[2+3*itri]=yc;
        zt[2+3*itri]=zc;

        xt[  3*mtri]=xa;
        yt[  3*mtri]=ya;
        zt[  3*mtri]=za;
        xt[1+3*mtri]=x1;
        yt[1+3*mtri]=y1;
        zt[1+3*mtri]=z1;
        xt[2+3*mtri]=xb;
        yt[2+3*mtri]=yb;
        zt[2+3*mtri]=zb;
        mtri++;

        xt[  3*mtri]=xb;
        yt[  3*mtri]=yb;
        zt[  3*mtri]=zb;
        xt[1+3*mtri]=x2;
        yt[1+3*mtri]=y2;
        zt[1+3*mtri]=z2;
        xt[2+3*mtri]=xc;
        yt[2+3*mtri]=yc;
        zt[2+3*mtri]=zc;
        mtri++;

        xt[  3*mtri]=xa;
        yt[  3*mtri]=ya;
        zt[  3*mtri]=za;
        xt[1+3*mtri]=xb;
        yt[1+3*mtri]=yb;
        zt[1+3*mtri]=zb;
        xt[2+3*mtri]=xc;
        yt[2+3*mtri]=yc;
        zt[2+3*mtri]=zc;
        mtri++;
       }
      ntri=mtri;
     }
   }

  for(itri=0;itri<ntri;itri++)
   {
    for(i=0;i<3;i++)
     {
      nrm[0+3*i]=     xt[i+3*itri];
      nrm[1+3*i]=     yt[i+3*itri];
      nrm[2+3*i]=     zt[i+3*itri];
      x[i]=(*xo)+(*r)*xt[i+3*itri];
      y[i]=(*yo)+(*r)*yt[i+3*itri];
      z[i]=(*zo)+(*r)*zt[i+3*itri];
     }
    n=3;
    if(ismooth!=0)shpgnrm(&n,x,y,z,nrm);
     else if(ismooth<0)shpl(&n,x,y,z);
     else shpg(&n,x,y,z);
   }

  free(xt);
  free(yt);
  free(zt);

  return;
 }

#ifdef __cplusplus
}
#endif
