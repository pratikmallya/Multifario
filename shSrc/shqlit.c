/* 
    @(#)shqlit.c	1.2
    02/04/19 16:40:28
   
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

void shqlit(int *m, float *x, float *y, float *z,int *r,int *g,int *b)
 {

/*    This Routine Returns the parameters associated with light source */

/*             x,y,z   (real) position.                                */
/*             r,g,b   (integer 0-255) color.                          */
/*             type    (integer 0-1) 0=point,1=direction.              */

  *r=sh_rs[*m];
  *g=sh_gs[*m];
  *b=sh_bs[*m];

  *x=sh_lit[ +3*(*m)];
  *y=sh_lit[1+3*(*m)];
  *z=sh_lit[2+3*(*m)];

  return;
 }

#ifdef __cplusplus
}
#endif
