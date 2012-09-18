/* 
    @(#)shcube.c	1.2
    02/04/19 16:37:59
   
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

void shcube(float *xmin,float *xmax,float *ymin,float *ymax,float *zmin,float *zmax)
 {

/*  Draw a wireframe cube. */

  shline(xmin,ymin,zmin,xmax,ymin,zmin);
  shline(xmax,ymin,zmin,xmax,ymin,zmax);
  shline(xmax,ymin,zmax,xmin,ymin,zmax);
  shline(xmin,ymin,zmax,xmin,ymin,zmin);

  shline(xmin,ymax,zmin,xmax,ymax,zmin);
  shline(xmax,ymax,zmin,xmax,ymax,zmax);
  shline(xmax,ymax,zmax,xmin,ymax,zmax);

  shline(xmin,ymax,zmax,xmin,ymax,zmin);

  shline(xmin,ymin,zmin,xmin,ymax,zmin);
  shline(xmin,ymin,zmax,xmin,ymax,zmax);

  shline(xmax,ymin,zmin,xmax,ymax,zmin);

  shline(xmax,ymin,zmax,xmax,ymax,zmax);

  return;
 }

#ifdef __cplusplus
}
#endif
