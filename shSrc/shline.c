/* 
    @(#)shline.c	1.2
    02/04/19 16:38:33
   
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

void shline(float *x0,float *y0,float *z0,float *x1,float *y1,float *z1)
 {

/*  Shade and render a line segment. */

  float nn[3]={0.,0.,0.};
  int zero=0;

  shlnonrm(x0,y0,z0,nn,x1,y1,z1,nn,&zero,&zero,&zero);

/*fprintf(stderr,"%e %e %e\n",*x0,*y0,*z0);
  fprintf(stderr,"%e %e %e\n",*x1,*y1,*z1);*/

  return;
 }

#ifdef __cplusplus
}
#endif
