/* 
    @(#)shinit.c	1.5
    02/04/19 16:38:24
   
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

void shinit(int *r,int *g,int *b)
 {
  int i,j;
  int full=255;
  int zero=0;
  float am;
  float ad;
  float as;
  int ae;
  float dist,alpha,beta,xmin,xmax;
  float ftwo=2.;
  int one=1;
  int two=2;
  float xp,yp,zp;
/*
      Initialization Routine.
            r,g,b  Background Color, integer 0-&full.
 
         This routine fills the image buffer iax with the
            supplied color.
         Then sets the z-buffer to be infinitly far. (z is scaled from
            zero to 1., so it is set to 2.)
 
         Two light sources are defined.
         No clipping planes.
         Default color is white (&full,&full,&full)
         Default view is looking down the x-axis, distance=10.
                                                  box=(0.,1.)**3
*/

  dist=10.;
  alpha=0.;
  beta=0.;
  xmin=0.;
  xmax=1.;
  shview(&dist,&alpha,&beta,&xmin,&xmax,&xmin,&xmax,&xmin,&xmax);

  shnlit(&two);
  xp=.3;
  yp=.3;
  zp=.5;
  shlit(&one,&xp,&yp,&zp,&full,&full,&full);
  xp=.7;
  yp=.3;
  zp=-1.;
  shlit(&two,&xp,&yp,&zp,&full,&full,&full);

  shnpln(&zero);
  shmask(&zero);

  shRedBackground=*r;
  shGreenBackground=*g;
  shBlueBackground=*b;

  shtric(&full,&full,&full);
  shlinc(&full,&full,&full);
  shpntc(&full,&full,&full);

/*    shsrfp(.3,.7,.5,5)*/
  am=.3;
  ad=.9;
  as=.3;
  ae=5;
  shsrfp(&am,&ad,&as,&ae);

  shRedBuffer=(unsigned char*)malloc(shIMax*shJMax*sizeof(unsigned char));
  if(shRedBuffer==(unsigned char*)NULL)
   {
    printf("Out of memory in shinit, allocating %d, line %d in file %s\n",shIMax*shJMax*sizeof(unsigned char),__LINE__,__FILE__);fflush(stdout);
    abort();
   }
  shGreenBuffer=(unsigned char*)malloc(shIMax*shJMax*sizeof(unsigned char));
  if(shGreenBuffer==(unsigned char*)NULL)
   {
    printf("Out of memory in shinit, allocating %d, line %d in file %s\n",shIMax*shJMax*sizeof(unsigned char),__LINE__,__FILE__);fflush(stdout);
    abort();
   }
  shBlueBuffer=(unsigned char*)malloc(shIMax*shJMax*sizeof(unsigned char));
  if(shBlueBuffer==(unsigned char*)NULL)
   {
    printf("Out of memory in shinit, allocating %d, line %d in file %s\n",shIMax*shJMax*sizeof(unsigned char),__LINE__,__FILE__);fflush(stdout);
    abort();
   }
  shZBuffer=(float*)malloc(shIMax*shJMax*sizeof(float));
  if(shZBuffer==(float*)NULL)
   {
    printf("Out of memory in shinit, allocating %d, line %d in file %s\n",shIMax*shJMax*sizeof(float),__LINE__,__FILE__);fflush(stdout);
    abort();
   }
/*printf("shinit, shZBuffer=0x%8.8x\n",shZBuffer);fflush(stdout);*/

  for(j=0;j<shJMax;j++)
   for(i=0;i<shIMax;i++)shsetp(&i,&j,r,g,b,&ftwo);

/*    shlss('CMR10',5,1000,0,ierr) */
/*    shcs(0,-3,3) */

  sh_frame=0;
  return;
 }

/* Lights */

#ifdef SHINITIALIZEHERE
int sh_nlit=0;
int sh_mlit=0;
int *sh_rs=(int*)NULL;
int *sh_gs=(int*)NULL;
int *sh_bs=(int*)NULL;
int *sh_type=(int*)NULL;
double *sh_lit;
#endif

/* The Frame Buffer */

unsigned char *shRedBuffer=(unsigned char*)NULL;
unsigned char *shGreenBuffer=(unsigned char*)NULL;
unsigned char *shBlueBuffer=(unsigned char*)NULL;
float *shZBuffer=(float*)NULL;
int shMax=480;
int shIMax=480;
int shJMax=480;

/* Color for Background */

int shRedBackground=0;
int shGreenBackground=0;
int shBlueBackground=0;

/* Color for Polylines */

#ifdef SHINITIALIZEHERE
int sh_rl=255;
int sh_gl=255;
int sh_bl=255;
#endif

/* Mask for Polygons */

#ifdef SHINITIALIZEHERE
int shMask=0;
#endif

/* Frame count for Saving Images */

int sh_frame=0;

/* Colors for Polygons */

#ifdef SHINITIALIZEHERE
int sh_rd[4]={255,255,255,255};
int sh_gd[4]={255,255,255,255};
int sh_bd[4]={255,255,255,255};
int sh_ram[4]={255,255,255,255};
int sh_gam[4]={255,255,255,255};
int sh_bam[4]={255,255,255,255};

/* Colors for Points */

int sh_rp=255;
int sh_gp=255;
int sh_bp=255;
#endif

/* Current point for shdraw */

float sh_x0=0.;
float sh_y0=0.;
float sh_z0=0.;

/* Surface Properties for Shading Polygons */

float sh_am=.3;
float sh_ad=.9;
float sh_as=.3;
int sh_nd=5;

/* Clipping Planes */

int sh_nplns=0;
int sh_mplns=0;
float *sh_plno;
float *sh_plnn;
int *sh_ipln=(int*)NULL;
int *sh_oper=(int*)NULL;

/* Fonts */

int sh_cFont=0;
int shadow_idx=0;
int shadow_idy=0;

int sh_kFonts=0;
int *sh_kfont=(int*)NULL;

int sh_nFonts=0;
int sh_mFonts=0;
int *sh_dir=(int*)NULL;

int *sh_fmag=(int*)NULL;
char **sh_fontnm=(char**)NULL;

int sh_lenimg=0;
int sh_nimg=0;
char *sh_image=(char*)NULL;

/* View */

float shv_eye[3]={1.,1.,1.};
float shv_d=10.;
float shv_l=100.;
float shv_n11=1.;
float shv_n12=0.;
float shv_n13=0.;
float shv_n21=0.;
float shv_n22=1.;
float shv_n23=0.;
float shv_n31=0.;
float shv_n32=0.;
float shv_n33=1.;
float shv_x0=0.;
float shv_y0=0.;
float shv_z0=0.;
float shv_pcale=1.;
float shv_soff=0.;
float shv_toff=0.;
float shv_dscl=1.;
float shv_doff=0.;

#ifdef __cplusplus
}
#endif
