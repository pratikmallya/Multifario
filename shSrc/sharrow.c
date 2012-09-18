/* 
    @(#)sharrow.c	1.3
    02/04/19 16:37:31
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <shInternal.h>
#include <math.h>

/*   maximg */

#define  f0(s) 1+s*s*(2*s-3)
#define  f1(s) -s*s*(2*s-3)
#define  f2(s) s*(s-1)*(s-1)
#define  f3(s) s*s*(s-1)
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#ifdef __cplusplus
 extern "C" {
#endif

void sharrow(float *x0,float *y0,float *z0,int *idx0,int *idy0,float *x1,float *y1,float *z1,int *idx1,int *idy1)
 {
  int rl,gl,bl;
  int orl=0;
  int ogl=0;
  int obl=0;
  float d,dx,dy;
  float s0=0.;
  float t0=0.;
  float d0=0.;
  float s1=0.;
  float t1=0.;
  float d1=0.;
  int ip;
  int i0,j0;
  int i1,j1;
  float dx0,dy0;
  float dx1,dy1;
  float s;
  int i,j,io,jo;
  int zero=0;
  int one=1;

  shpers(x0,y0,z0,&s0,&t0,&d0);
  i0=round(s0*shMax);
  j0=round(t0*shMax);

  shpers(x1,y1,z1,&s1,&t1,&d1);
  i1=round(s1*shMax);
  j1=round(t1*shMax);

  shqlinc(&orl,&ogl,&obl);

  d=1./sqrt((*idx0)*(*idx0)+(*idy0)*(*idy0));
  dx0=(*idx0)*d >= 0 ? (int)((*idx0)*d+0.5) : (int)((*idx0)*d-0.5);
  dy0=(*idy0)*d >= 0 ? (int)((*idy0)*d+0.5) : (int)((*idy0)*d-0.5);
  d=1./sqrt((*idx1)*(*idx1)+(*idy1)*(*idy1));
  dx1=(*idx1)*d;
  dy1=(*idy1)*d;

  ip=3;
  for(s=0.;s<1.;s+=.01)
   {
    i=f0(s)*i0+f1(s)*i1+f2(s)*dx0+f3(s)*dx1 >= 0 ? (int)( f0(s)*i0+f1(s)*i1+f2(s)*dx0+f3(s)*dx1 +0.5) : (int)( f0(s)*i0+f1(s)*i1+f2(s)*dx0+f3(s)*dx1 -0.5);
    j=f0(s)*j0+f1(s)*j1+f2(s)*dy0+f3(s)*dy1 >= 0 ? (int)( f0(s)*j0+f1(s)*j1+f2(s)*dy0+f3(s)*dy1 +0.5) : (int)( f0(s)*j0+f1(s)*j1+f2(s)*dy0+f3(s)*dy1 -0.5);
    if(ip==2)
     {
      shlinc(&rl,&gl,&bl);
      shline2s(&io,&jo,&i,&j,&one,&one);
     }else{
      rl=orl;
      gl=ogl;
      bl=obl;
      shline2s(&io,&jo,&i,&j,&zero,&zero);
     }
    ip=2;
    io=i;
    jo=j;
   }

  rl=0;
  gl=0;
  bl=0;
  io=i1;
  jo=j1;
  i=round(i1-(*idx1)+.3*(*idy0));
  j=round(j1-(*idy0)-.3*(*idx0));
  shline2s(&io,&jo,&i,&j,&one,&one);
  rl=orl;
  gl=ogl;
  bl=obl;
  shline2s(&io,&jo,&i,&j,&zero,&zero);

  rl=0;
  gl=0;
  bl=0;
  i=round(i1-(*idx1)-.15*(*idy0));
  j=round(j1-(*idy0)+.15*(*idx0));
  shline2s(&io,&jo,&i,&j,&one,&one);
  rl=orl;
  gl=ogl;
  bl=obl;
  shline2s(&io,&jo,&i,&j,&zero,&zero);

  rl=0;
  gl=0;
  bl=0;
  io=round(i1-(*idx1)-.15*(*idy0));
  jo=round(j1-(*idy0)+.15*(*idx0));
  i=round(i1-(*idx1)+.15*(*idy0));
  j=round(j1-(*idy0)-.15*(*idx0));
  shline2s(&io,&jo,&i,&j,&one,&one);
  rl=orl;
  gl=ogl;
  bl=obl;
  shline2s(&io,&jo,&i,&j,&zero,&zero);

  return;
 }

#ifdef __cplusplus
}
#endif
