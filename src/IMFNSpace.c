/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   December 9, 2002   modified InvMF
 */

static char *id="@(#) $Id: IMFNSpace.c,v 1.6 2011/07/21 17:42:46 mhender Exp $";

static char IMFNSpaceErrorMsg[256]="";

#include <MFNSpace.h>
#include <MFNVector.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

static double MFInvMFNSpaceDistance(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFInvMFNSpaceDirection(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFInvMFNSpaceAdd(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFInvMFNSpaceScale(MFNSpace,double,MFNVector,MFNVector,void*,MFErrorHandler);
static double MFInvMFNSpaceInner(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
static void MFInvMFNSpaceFreeData(void*,MFErrorHandler);

static void MFWriteInvMFNSpace(FILE*,MFNSpace,void*,MFErrorHandler);
static MFNSpace MFReadInvMFNSpace(FILE*,MFErrorHandler);

extern double scaled_dp(double*,double*,void*);

struct InvMFNSpaceData
 {
  int n;
 };

/*! \fn MFNSpace InvMFCreateNSpace(int n, MFErrorHandler e);
 *  \brief Creates an inner product space for the intersection of a sphere and expansion.
 *
 *  \param n The dimension of the space.
 *  \param e An error handler.
 *  \returns A new MFNSpace
 */
MFNSpace InvMFCreateNSpace(int n, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateInvMFNSpace"};
  MFNSpace thisSpace;
  struct InvMFNSpaceData *data;


  thisSpace=MFCreateNSpaceBaseClass("InvMFNSpace",e);

  data=(struct InvMFNSpaceData*)malloc(sizeof(struct InvMFNSpaceData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(IMFNSpaceErrorMsg,"Out of memory trying to allocate %d bytes\n",sizeof(struct InvMFNSpaceData));
    MFSetError(e,12,RoutineName,IMFNSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  data->n=n;
  MFNSpaceSetData(thisSpace,data,e);
  MFNSpaceSetFreeData(thisSpace,MFInvMFNSpaceFreeData,e);

  MFNSpaceSetDistance(thisSpace,MFInvMFNSpaceDistance,e);
  MFNSpaceSetInnerProduct(thisSpace,MFInvMFNSpaceInner,e);
  MFNSpaceSetDirection(thisSpace,MFInvMFNSpaceDirection,e);
  MFNSpaceSetAdd(thisSpace,MFInvMFNSpaceAdd,e);
  MFNSpaceSetScale(thisSpace,MFInvMFNSpaceScale,e);

  return thisSpace;
 }

double MFInvMFNSpaceDistance(MFNSpace thisSpace,MFNVector v0,MFNVector v1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFInvMFNSpaceDistance"};
  double result;
  MFNVector dv;

  dv=MFCloneNVector(v0,e);
  MFInvMFNSpaceDirection(thisSpace,v0,v1,dv,d,e);
  result=MFInvMFNSpaceInner(thisSpace,dv,dv,d,e);
  MFFreeNVector(dv,e);

  return sqrt(result);
 }

void MFInvMFNSpaceDirection(MFNSpace thisSpace,MFNVector v0,MFNVector v1,MFNVector diff,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFInvMFNSpaceDirection"};
  int i;


  MFNVDiff(v1,v0,diff,e);

  return;
 }

double MFInvMFNSpaceInner(MFNSpace thisSpace,MFNVector v0,MFNVector v1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFInvMFNSpaceInner"};
  double result;
  double *x0;
  double *x1;
  int i,n;
  struct InvMFNSpaceData *data;

  data=(struct InvMFNSpaceData *)d;

  x0=MFNV_CStar(v0,e);
  x1=MFNV_CStar(v1,e);
  n=data->n;
  result=0.;for(i=0;i<n;i++)result+=x0[i]*x1[i];

  return result;
 }

void MFWriteInvMFNSpace(FILE *fid,MFNSpace thisSpace,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteInvMF"};
  int i;

  return;
 }

void MFInvMFNSpaceAdd(MFNSpace thisSpace,MFNVector v0,MFNVector v1,MFNVector sum,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFInvMFNSpaceAdd"};
  int i;

  MFNVAdd(v0,v1,sum,e);

  return;
 }

void MFInvMFNSpaceScale(MFNSpace thisSpace,double s, MFNVector v,MFNVector w,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFInvMFNSpaceScale"};
  double *pv;
  double *pw;
  int i,n;

  pv=MFNV_CStar(v,e);
  pw=MFNV_CStar(w,e);
  n=MFNV_NC(v,e);

  for(i=0;i<n;i++)pw[i]=s*pv[i];

  return;
 }

void MFInvMFNSpaceFreeData(void *d, MFErrorHandler e)
 {
  struct InvMFNSpaceData *data;

  data=(struct InvMFNSpaceData *)d;
  if(data!=NULL)free(data);
 }

#ifdef __cplusplus
}
#endif
