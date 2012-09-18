/* 
    @(#)shpgc.c	1.3
    02/04/19 16:39:39
   
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

int sh_rd[4]={0,0,0,0};
int sh_gd[4]={0,0,0,0};
int sh_bd[4]={0,0,0,0};
int sh_ram[4]={0,0,0,0};
int sh_gam[4]={0,0,0,0};
int sh_bam[4]={0,0,0,0};

void shpgc(int *r,int *g,int *b)
 {
  sh_rd[0] =*r;
  sh_gd[0] =*g;
  sh_bd[0] =*b;
  sh_ram[0]=*r;
  sh_gam[0]=*g;
  sh_bam[0]=*b;
  return;
 }

void shpec(int *r,int *g,int *b)
 {
  sh_rd[1] =*r;
  sh_gd[1] =*g;
  sh_bd[1] =*b;
  sh_ram[1]=*r;
  sh_gam[1]=*g;
  sh_bam[1]=*b;
  return;
 }

void shbgc(int *r,int *g,int *b)
 {
  sh_rd[2] =*r;
  sh_gd[2] =*g;
  sh_bd[2] =*b;
  sh_ram[2]=*r;
  sh_gam[2]=*g;
  sh_bam[2]=*b;
  return;
 }

void shbec(int *r,int *g,int *b)
 {
  sh_rd[3] =*r;
  sh_gd[3] =*g;
  sh_bd[3] =*b;
  sh_ram[3]=*r;
  sh_gam[3]=*g;
  sh_bam[3]=*b;
  return;
 }

void shqpgc(int *r,int *g,int *b)
 {
  *r=sh_rd[0];
  *g=sh_gd[0];
  *b=sh_bd[0];
  return;
 }

void shqpec(int *r,int *g,int *b)
 {
  *r=sh_rd[1];
  *g=sh_gd[1];
  *b=sh_bd[1];
  return;
 }

void shqbgc(int *r,int *g,int *b)
 {
  *r=sh_rd[2];
  *g=sh_gd[2];
  *b=sh_bd[2];
  return;
 }

void shqbec(int *r,int *g,int *b)
 {
  *r=sh_rd[3];
  *g=sh_gd[3];
  *b=sh_bd[3];
  return;
 }

#ifdef __cplusplus
}
#endif
