/* 
    @(#)shInternal.h	1.5
    02/04/19 14:43:17
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#ifndef __SHINTERNAL_H__
#define __SHINTERNAL_H__
#include <sh.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

/* Lights */

extern int sh_nlit;
extern int sh_mlit;
extern int *sh_rs;
extern int *sh_gs;
extern int *sh_bs;
extern int *sh_type;
extern double *sh_lit;

/* The Frame Buffer */

extern unsigned char *shRedBuffer;
extern unsigned char *shGreenBuffer;
extern unsigned char *shBlueBuffer;
extern float *shZBuffer;
extern int shMax;
extern int shIMax;
extern int shJMax;

/* Color for Background */

extern int shRedBackground;
extern int shGreenBackground;
extern int shBlueBackground;

/* Color for Polylines */

extern int sh_rl;
extern int sh_gl;
extern int sh_bl;

/* Mask for Polygons */

extern int shMask;

/* Frame count for Saving Images */

extern int sh_frame;

/* Colors for Polygons */

extern int sh_rd[];
extern int sh_gd[];
extern int sh_bd[];
extern int sh_ram[];
extern int sh_gam[];
extern int sh_bam[];

/* Colors for Points */

extern int sh_rp;
extern int sh_gp;
extern int sh_bp;

/* Current point for shdraw */

extern float sh_x0;
extern float sh_y0;
extern float sh_z0;

/* Surface Properties for Shading Polygons */

extern float sh_am;
extern float sh_ad;
extern float sh_as;
extern int sh_nd;

/* Clipping Planes */

extern int sh_nplns;
extern int sh_mplns;
extern float *sh_plno;
extern float *sh_plnn;
extern int *sh_ipln;
extern int *sh_oper;

/* Fonts */

extern int sh_cFont;
extern int shadow_idx;
extern int shadow_idy;

extern int sh_kFonts;
extern int *sh_kfont;

extern int sh_nFonts;
extern int sh_mFonts;
extern int *sh_dir;

extern int *sh_fmag;
extern char **sh_fontnm;

extern int sh_lenimg;
extern int sh_nimg;
extern char *sh_image;

/* View */

extern float shv_eye[];
extern float shv_d;
extern float shv_l;
extern float shv_n11;
extern float shv_n12;
extern float shv_n13;
extern float shv_n21;
extern float shv_n22;
extern float shv_n23;
extern float shv_n31;
extern float shv_n32;
extern float shv_n33;
extern float shv_x0;
extern float shv_y0;
extern float shv_z0;
extern float shv_pcale;
extern float shv_soff;
extern float shv_toff;
extern float shv_dscl;
extern float shv_doff;

extern char *shOutputFormat;

void shputps(char*,int*,int*,int*,unsigned char*,unsigned char*,unsigned char*,int);
void shputtiff(char*,int,int,unsigned char*,unsigned char*,unsigned char*);
void sh3dmask(float*,float*,float*,float*,int*);
void shdopg(int*,float*,float*,float*,float*);
void shgetz(int*,int*,float*);
void shlss(char*,int*,int*,int*,int*,int);
void shmask(int*);
void shcolor(float*,float*,float*,float*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
void shclippg(int*,float*,float*,float*,float*,int*,int*,float**,float**,float**,float**,int**);
void shsetp(int*,int*,int*,int*,int*,float*);
void shdoStrings(FILE*,float,float);
void shfreeStrings();
void shlineoff(float*,float*,float*,float*,float*,float*,int*,int*,int*);
void shpers(float*,float*,float*,float*,float*,float*);
void shunpers(float*,float*,float*,float*,float*,float*);
void shsize(float*,float*,float*,float*,float*);

#endif

#ifdef __cplusplus
}
#endif
