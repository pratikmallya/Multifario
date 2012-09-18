/* 
    @(#)shpln.c	1.2
    02/04/19 16:39:59
   
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

void shpln(int *i,float *ox,float *oy,float *oz,float *nx,float *ny,float *nz,int *iside)
 {
  int ip;

  ip=(*i)-1;

/* Set the clipping plane info. */

  if(*i>=sh_mplns)
   {
    if(sh_mplns==0)
     {
      sh_mplns=10;
      sh_plno=(float*)malloc(3*sh_mplns*sizeof(float));
      sh_plnn=(float*)malloc(3*sh_mplns*sizeof(float));
      sh_ipln=(int*)malloc(sh_mplns*sizeof(int));
      sh_oper=(int*)malloc(sh_mplns*sizeof(int));
     }else{
      sh_mplns+=10;
      sh_plno=(float*)realloc((void*)sh_plno,3*sh_mplns*sizeof(float));
      sh_plnn=(float*)realloc((void*)sh_plnn,3*sh_mplns*sizeof(float));
      sh_ipln=(int*)realloc((void*)sh_ipln,sh_mplns*sizeof(int));
      sh_oper=(int*)realloc((void*)sh_oper,sh_mplns*sizeof(int));
     }
   }

  sh_plno[  3*ip]=*ox;
  sh_plno[1+3*ip]=*oy;
  sh_plno[2+3*ip]=*oz;

  sh_plnn[  3*ip]=*nx;
  sh_plnn[1+3*ip]=*ny;
  sh_plnn[2+3*ip]=*nz;

  sh_ipln[ip]=*iside;
  sh_oper[ip]=0;

/*printf("plane %d is ((x,y,z)-(%f,%f,%f)).(%f,%f,%f)*%d>0\n",ip,sh_plno[  3*ip],sh_plno[1+3*ip],sh_plno[2+3*ip],sh_plnn[  3*ip],sh_plnn[1+3*ip],sh_plnn[2+3*ip],sh_ipln[ip]);fflush(stdout);*/


  return;
 }

#ifdef __cplusplus
}
#endif
