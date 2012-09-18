/* 
    @(#)shpers.c	1.2
    02/04/19 16:39:24
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#include <math.h>

#ifdef __cplusplus
 extern "C" {
#endif

void shpers(float *x,float *y,float *z,float *s,float *t,float *depth)
 {
  float c1,c2,c3;
  float dtemp;

/*    Project x,y,z to perspective coordinates. */

  c1= (*x-shv_x0)*shv_n11+(*y-shv_y0)*shv_n12+(*z-shv_z0)*shv_n13;
  c2= (*x-shv_x0)*shv_n21+(*y-shv_y0)*shv_n22+(*z-shv_z0)*shv_n23;
  c3= (*x-shv_x0)*shv_n31+(*y-shv_y0)*shv_n32+(*z-shv_z0)*shv_n33;

  dtemp=(shv_l-shv_d)/(c1-shv_d);
  *s=dtemp*c2*shv_pcale+shv_soff;
  *t=dtemp*c3*shv_pcale+shv_toff;
  *depth=-dtemp*c1*shv_dscl+shv_doff;

  return;
 }

void shunpers(float *ss,float *tt,float *dd,float *xx,float *yy,float *zz)
 {
  float c1,c2,c3;

/*    Entry to invert perspective transform. */

  c1=shv_d*(*dd-shv_doff)/(*dd-shv_doff+shv_dscl*(shv_l-shv_d));
  c2=(*ss-shv_soff)*(c1-shv_d)/(shv_l-shv_d)/shv_pcale;
  c3=(*tt-shv_toff)*(c1-shv_d)/(shv_l-shv_d)/shv_pcale;

  *xx=shv_x0+c1*shv_n11+c2*shv_n21+c3*shv_n31;
  *yy=shv_y0+c1*shv_n12+c2*shv_n22+c3*shv_n32;
  *zz=shv_z0+c1*shv_n13+c2*shv_n23+c3*shv_n33;

  return;
 }

void shvw(float *dist,float *alpha,float *beta,float *xmin,float *xmax,float *ymin,float *ymax,float *zmin,float *zmax)
 {
  float a,b,pi;

/*    Entry to set view point. */

  shv_x0=((*xmin)+(*xmax))/2.;
  shv_y0=((*ymin)+(*ymax))/2.;
  shv_z0=((*zmin)+(*zmax))/2.;

  pi=3.1415926;
  a=(*alpha-90.)*pi/180.;
  b=*beta*pi/180.;

  shv_n11=cos(b)*cos(a);
  shv_n12=cos(b)*sin(a);
  shv_n13=sin(b);

  shv_eye[0]=*dist*shv_n11;
  shv_eye[1]=*dist*shv_n12;
  shv_eye[2]=*dist*shv_n13;

  shv_n21=-sin(a);
  shv_n22=cos(a);
  shv_n23=0.0;

  shv_n31=-sin(b)*cos(a);
  shv_n32=-sin(b)*sin(a);
  shv_n33=cos(b);

  shv_d=*dist;
  shv_l=0.75*shv_d;
  return;
 }

void shqeye(float *xe,float *ye,float *ze)
 {

/*    Entry to return viewpoint. */

  *xe=shv_d*shv_n11;
  *ye=shv_d*shv_n12;
  *ze=shv_d*shv_n13;
  return;
 }

void shsize(float *scl,float *s0,float *t0,float *dsc,float *d0)
 {

/*    Entry to set scale factors. */

  shv_pcale=*scl;
  shv_soff=*s0;
  shv_toff=*t0;
  shv_doff=*d0;
  shv_dscl=*dsc;

  return;
 }

#ifdef __cplusplus
}
#endif
