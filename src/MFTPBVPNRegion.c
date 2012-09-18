/* 
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   November 11, 1997
 *              February 2, 1999   Ported to C
 *              October 6, 2000    Added Polyhedron
 */

static char *id="@(#) $Id: MFTPBVPNRegion.c,v 1.3 2007/02/13 01:22:34 mhender Exp $";

static char MFNRegionErrorMsg[256]="";

#include <MFNVector.h>
#include <MFNRegion.h>
#include <MFErrorHandler.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
 extern "C" {
#endif

static int TPBVPTest(MFNVector,void*,MFErrorHandler);
static void TPBVPFree(void*,MFErrorHandler);

MFNRegion MFNRegionCreateTPBVP(int nx,int nu,int np,double *p0,double *p1,double u0, double u1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreateTPBVP"};
  MFNRegion tpbvp;
  int i;
  double *data;

  tpbvp=MFNRegionCreateBaseClass("TPBVP",e);
  MFNRegionSetTest(tpbvp,TPBVPTest,e);
  MFNRegionSetFreeData(tpbvp,TPBVPFree,e);

  data=(double*)malloc((5+2*np)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",(5+2*np)*sizeof(double));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(tpbvp);
    return NULL;
   }
#endif

  ((int*)(data))[0]=nx;
  ((int*)(data))[1]=nu;
  ((int*)(data))[2]=np;
  for(i=0;i<np;i++)
   {
    ((double*)(data))[3+2*i]=p0[i];
    ((double*)(data))[3+2*i+1]=p1[i];
   }
  ((double*)(data))[3+2*np]=u0;
  ((double*)(data))[3+2*np+1]=u1;

  MFNRegionSetData(tpbvp,data,e);

  return(tpbvp);
 }

int TPBVPTest(MFNVector v,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"TPBVPTest"};
  double *d;
  double u;
  int rc;
  int n;
  double *x;
  int nx,nu,np;
  int i,j;
  static int verbose=0;

  n=MFNV_NC(v,e);
  x=MFNV_CStar(v,e);

  d=(double*)data;
  nx=((int*)d)[0];
  nu=((int*)d)[1];
  np=((int*)d)[2];
  rc=1;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s. nx=%d, nu=%d, np=%d\n",RoutineName,nx,nu,np);}
#endif

  for(i=0;i<np;i++)
   {
    rc=rc&&(x[nx*nu+i]>=d[3+2*i])&&(x[nx*nu+i]<=d[3+2*i+1]);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("%d, %lf<=%lf<=%lf\n",i,d[3+2*i],x[nx*nu+i],d[3+2*i+1]);fflush(stdout);}
#endif

   }
  u=0;
  i=0;
  for(j=0;j<nu;j++)u+=(x[i*nu+j]+x[(i+1)*nu+j])*x[nx*nu+np+i]/2;
  for(i=0;i<nx-1;i++)
   {
    for(j=0;j<nu;j++)
     {
      u+=(x[i*nu+j]+x[(i+1)*nu+j])*(x[nx*nu+np+i+1]-x[nx*nu+np+i])/2;
     }
   }
  i=nx-1;
  for(j=0;j<nu;j++)u+=(x[i*nu+j]+x[(i+1)*nu+j])*x[nx*nu+np+i]/2;
  u=u/nu;
  rc=rc&&(u>=d[3+2*np])&&(u<=d[3+2*np+1]);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  %lf<=%lf<=%lf\n",d[3+2*np],u,d[3+2*np+1]);fflush(stdout);}
#endif

  return rc;
 }

void TPBVPFree(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"TPBVPFree"};

  free(data);
  return;
 }

#ifdef __cplusplus
}
#endif
