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
 *              October 6, 2000   Added Polyhedron
 *              September 22, 2005   Ported from MFCSGBallsNRegion
 */

static char *id="@(#) $Id: MFCSGBallsNRegion.c,v 1.3 2007/02/13 01:22:33 mhender Exp $";

#include <stdio.h>
#include <MFNVector.h>
#include <MFNRegion.h>
#include <stdlib.h>
#include <math.h>

static char MFNRegionErrorMsg[256]="";

#ifdef __cplusplus
 extern "C" {
#endif

int CSGBallsTest(MFNVector,void*,MFErrorHandler);
void CSGBallsFree(void*,MFErrorHandler);

struct MFCSGBallData
 {
  int n;
  int nc;
  double *x0;
  double *R;
  int *dir;
 };

MFNRegion MFNRegionCreateCSGBalls(int n,int nc,double *x0,double *R,int *dir, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreateCSGBalls"};
  MFNRegion csg;
  int i;
  struct MFCSGBallData *data;

  csg=MFNRegionCreateBaseClass("CSGBalls",e);
  MFNRegionSetTest(csg,CSGBallsTest,e);
  MFNRegionSetFreeData(csg,CSGBallsFree,e);

  data=(struct MFCSGBallData*)malloc(sizeof(struct MFCSGBallData));
#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFCSGBallData));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->n=n;
  data->nc=nc;

  data->x0=(double*)malloc(n*nc*sizeof(double));
#ifndef MFNOSAFETYNET
  if(data->x0==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",n*nc*sizeof(double));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<n*nc;i++)(data->x0)[i]=x0[i];
  data->R=(double*)malloc(nc*sizeof(double));
#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",nc*sizeof(double));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<nc;i++)(data->R)[i]=R[i];

  data->dir=(int*)malloc(nc*sizeof(int));
#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",nc*sizeof(int));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<nc;i++)(data->dir)[i]=dir[i];

  MFNRegionSetData(csg,data,e);

  return(csg);
 }

int CSGBallsTest(MFNVector v,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"CSGBallsTest"};
  struct MFCSGBallData *d;
  double r;
  double *x;
  int rc;
  int i,j;
  static int verbose=0;

  d=(struct MFCSGBallData*)data;
  x=MFNV_CStar(v,e);

  rc=1;
  for(i=0;i<d->nc;i++)
   {
    r=0.;
    for(j=0;j<d->n;j++)r+=(x[j]-d->x0[j+d->n*i])*(x[j]-d->x0[j+d->n*i]);
    rc=rc&&d->dir[i]*(d->R[i]*d->R[i]-r*r)<0.;
   }

  return rc;
 }

void CSGBallsFree(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"CSGBallsFree"};
  struct MFCSGBallData *d;

  d=(struct MFCSGBallData*)data;

  free(d->x0);
  free(d->R);
  free(d->dir);
  free(d);

  return;
 }

#ifdef __cplusplus
}
#endif
