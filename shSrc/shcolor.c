/* 
    @(#)shcolor.c	1.3
    02/04/19 16:37:51
   
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

void shcolor(float *x,float *y,float *z,float *n,int *rf,int *gf,int *bf,int *rb,int *gb,int *bb,int *r,int *g,int *b)
 {

/*  Shade and render a point. */

   int rd[2]={0,0};
   int gd[2]={0,0};
   int bd[2]={0,0};

/* Geometric Information */

/*    EYE[2]    position of the viewer                            */
/*    NLIT      number of light sources                           */
/*    LIT(3,M)  position of light source m                        */
/*    RS(M)     color of light source m                           */
/*    GS(M)                                                       */
/*    BS(M)                                                       */
/*    AM        intensity of ambient illumination                 */
/*    AD        intensity of diffuse illumination                 */
/*    AS        intensity of specular illumination                */
/*    ND        exponent  of specular illumination                */

/* Clipping Planes                            */

/*    NPLNS     Number of clipping planes */
/*    PLNO(3,M) point on plane M */
/*    PLNN(3,M) normal to plane M */
/*    IPLN(M)   side of plane M to clip (IPLN*PLNN<0 ==> clipped) */
/*    oper(M)   whether clipping is ored or anded. */

  float e[3]={0.,0.,0.};
  float l[3]={0.,0.,0.};
  float h[3]={0.,0.,0.};

  float an;
  float ah;
  float ae;
  float sd;
  float ssr;
  float ssg;
  float ssb;
  int   k;
  float xxx=0.;
  float yyy=0.;
  float zzz=0.;
  float al;
  int icol;

  an=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  if(an<1.e-7)
   {
    *r=*rf;
    *g=*gf;
    *b=*bf;
    return;
   }

/*             e is the vector from the pixel to the eye */

  e[0]=shv_eye[0]-*x;
  e[1]=shv_eye[1]-*y;
  e[2]=shv_eye[2]-*z;
  ae=1./sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
  e[0]=e[0]*ae;
  e[1]=e[1]*ae;
  e[2]=e[2]*ae;

  sd =0.;
  ssr=0.;
  ssg=0.;
  ssb=0.;
  for(k=0;k<sh_nlit;k++)
   {

/*             l is the vector from the pixel to light source k */

    if(sh_type[k]==0)
     {
      l[0]=sh_lit[  3*k]-xxx;
      l[1]=sh_lit[1+3*k]-yyy;
      l[2]=sh_lit[2+3*k]-zzz;
     }else{
      l[0]=sh_lit[  3*k];
      l[1]=sh_lit[1+3*k];
      l[2]=sh_lit[2+3*k];
     }
    al=1./sqrt(l[0]*l[0]+l[1]*l[1]+l[2]*l[2]);
    l[0]=l[0]*al;
    l[1]=l[1]*al;
    l[2]=l[2]*al;
/*  printf("light %d, l=(%f,%f,%f) n=(%f,%f,%f\n",k,l[0],l[1],l[2],n[0],n[1],n[2]);*/

/*             h will control the specular reflection */

    h[0]=e[0]+l[0];
    h[1]=e[1]+l[1];
    h[2]=e[2]+l[2];
    ah=1./sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]);
    h[0]=h[0]*ah;
    h[1]=h[1]*ah;
    h[2]=h[2]*ah;
/*  printf("h=(%f,%f,%f)\n",h[0],h[1],h[2]);*/

/*  printf("sd=%f sscol=(%f,%f,%f)\n",fabs(n[0]*l[0]+n[1]*l[1]+n[2]*l[2]),sh_rs[k]*pow(fabs(n[0]*h[0]+n[1]*h[1]+n[2]*h[2]),sh_nd),sh_gs[k]*pow(fabs(n[0]*h[0]+n[1]*h[1]+n[2]*h[2]),sh_nd),sh_bs[k]*pow(fabs(n[0]*h[0]+n[1]*h[1]+n[2]*h[2]),sh_nd));*/

    sd =sd +fabs(n[0]*l[0]+n[1]*l[1]+n[2]*l[2]);
    ssr=ssr+sh_rs[k]*pow(fabs(n[0]*h[0]+n[1]*h[1]+n[2]*h[2]),sh_nd);
    ssg=ssg+sh_gs[k]*pow(fabs(n[0]*h[0]+n[1]*h[1]+n[2]*h[2]),sh_nd);
    ssb=ssb+sh_bs[k]*pow(fabs(n[0]*h[0]+n[1]*h[1]+n[2]*h[2]),sh_nd);
   }
/*printf("total sd=%f sscol=(%f,%f,%f)\n",sd,ssr,ssg,ssb);*/

  icol=0;
  if(n[0]*(shv_eye[0]-*x)+n[1]*(shv_eye[1]-*y)+n[2]*(shv_eye[2]-*z)<0.)icol=1;
/*printf("side %d\n",icol);*/

  rd[0]=*rf;
  gd[0]=*gf;
  bd[0]=*bf;
  rd[1]=*rb;
  gd[1]=*gb;
  bd[1]=*bb;

  *r=round(sh_am*rd[icol]+sh_ad*rd[icol]*sd/sh_nlit+sh_as*ssr/sh_nlit);
  *g=round(sh_am*gd[icol]+sh_ad*gd[icol]*sd/sh_nlit+sh_as*ssg/sh_nlit);
  *b=round(sh_am*bd[icol]+sh_ad*bd[icol]*sd/sh_nlit+sh_as*ssb/sh_nlit);

  if(*r>255)*r=255;
  if(*g>255)*g=255;
  if(*b>255)*b=255;
  if(*r<0)*r=0;
  if(*g<0)*g=0;
  if(*b<0)*b=0;

  return;
 }

#ifdef __cplusplus
}
#endif
