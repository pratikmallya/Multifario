/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 */

static char *id="@(#) $Id: MFComputePointsOn3dEdge.c,v 1.4 2007/02/13 01:22:33 mhender Exp $";

static char MFCompute3dEdgeErrorMsg[256]="";

#include <MFAtlas.h>
#include <MFAtlasFriends.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <MFErrorHandler.h>
#include <sh.h>
#include <math.h>

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#ifdef __cplusplus
 extern "C" {
#endif

int MFComputePointsOn3dEdge(double r,double **x,double xa,double ya,double za,double xb,double yb,double zb, MFErrorHandler e)
 {
  static char RoutineName[]={"MFComputePointsOn3dEdge"};
  double R,t;
  int i,n;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("MFComputePointsOn3dEdge(%lf,%8.8x,%lf,%lf,%lf,%lf,%lf,%lf)\n",r,*x,xa,ya,za,xb,yb,zb);fflush(stdout);}
#endif

  R=sqrt(pow(xa-xb,2)+pow(ya-yb,2)+pow(za-zb,2));
  n=round(R/r);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  allocating space for %d points\n",n);fflush(stdout);}
#endif

  *x=(double*)realloc((void*)(*x),3*n*sizeof(double));
#ifndef MFNOSAFETYNET
  if(*x==NULL)
   {
    sprintf(MFCompute3dEdgeErrorMsg,"Out of memory, trying to allocate %d bytes",3*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFCompute3dEdgeErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose) 
   {
    printf("  space at %8.8x\n",*x);fflush(stdout);
    printf("  Distance between points is %lf, r=%lf\n",R/n,r);fflush(stdout);
   }
#endif

  for(i=0;i<n;i++)
   {
    t=(i+1)*1./(n+1);
    (*x)[3*i  ]=xa+t*(xb-xa);
    (*x)[3*i+1]=ya+t*(yb-ya);
    (*x)[3*i+2]=za+t*(zb-za);
   }

  return n;
 }

#ifdef __cplusplus
}
#endif
