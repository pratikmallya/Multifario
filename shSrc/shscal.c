/* 
    @(#)shscal.c	1.3
    02/04/19 16:40:51
   
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

void shscal(float *xmin,float *xmax,float *ymin,float *ymax,float *zmin,float *zmax)
 {
  int i,j;

  float one=1.;
  float zero=0.;
  float smax,smin;
  float tmax,tmin;
  float dmax,dmin;
  float s=0.;
  float t=0.;
  float d=0.;
  float x,y,z;
  float as,at,pcale;
  float soff;
  float toff;
  float dcale;
  float doff;

  float I[8*6]={0.,0.,0.,0.,1.,1.,1.,1.,
                1.,1.,1.,1.,0.,0.,0.,0.,
                0.,0.,1.,1.,0.,0.,1.,1.,
                1.,1.,0.,0.,1.,1.,0.,0.,
                0.,1.,0.,1.,0.,1.,0.,1.,
                1.,0.,1.,0.,1.,0.,1.,0.};

/* Perform Scaling for the perspective transform. */

  shsize(&one,&zero,&zero,&one,&zero);

  smax=-1.e5;
  smin=1.e5;
  tmax=-1.e5;
  tmin=1.e5;
  dmax=-1.e5;
  dmin=1.e5;

  for(j=0;j<8;j++)
   {
    x=I[j+8*0]*(*xmax)+I[j+8*1]*(*xmin);
    y=I[j+8*2]*(*ymax)+I[j+8*3]*(*ymin);
    z=I[j+8*4]*(*zmax)+I[j+8*5]*(*zmin);
    shpers(&x,&y,&z,&s,&t,&d);
    if(s>smax)smax=s;
    if(s<smin)smin=s;
    if(t>tmax)tmax=t;
    if(t<tmin)tmin=t;
    if(d>dmax)dmax=d;
    if(d<dmin)dmin=d;
   }

  as=1.*shIMax/shMax;
  at=1.*shJMax/shMax;
  pcale=.9*as/(smax-smin);
  if(pcale*(tmax-tmin)>.9*at)pcale=.9*at/(tmax-tmin);
  soff=.5*(as-pcale*(smax+smin));
  toff=.5*(at-pcale*(tmax+tmin));
  dcale=.9/(dmax-dmin);
  doff=.5*(1.-dcale*(dmax+dmin));
  shsize(&pcale,&soff,&toff,&dcale,&doff);
  return;
 }

#ifdef __cplusplus
}
#endif
