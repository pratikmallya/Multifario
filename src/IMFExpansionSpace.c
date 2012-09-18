/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   March 1, 1999
 */

static char *id="@(#) $Id: IMFExpansionSpace.c,v 1.6 2011/07/21 17:42:46 mhender Exp $";

static char IMFExpansionSpaceErrorMsg[256]="";

#include <MFNSpace.h>
#include <MFNVector.h>
#include <MFErrorHandler.h>
#include <MFPrint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct IMFExpansionSpaceData
 {
  int n;
 };

static void MFFreeExpansionSpaceData(void*,MFErrorHandler);
static double IMFExpansionSpaceDistance(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);
static void IMFExpansionSpaceDirection(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
static void IMFExpansionSpaceAdd(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler);
static void IMFExpansionSpaceScale(MFNSpace,double,MFNVector,MFNVector,void*,MFErrorHandler);
static double IMFExpansionSpaceInner(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler);

/*! \fn MFNSpace IMFCreateExpansionSpace(int n, MFErrorHandler e);
 *  \brief An MFNSpace defining operations on a point on a fat trajectory.
 *
 *  \param n The dimension of the embedding space of the expansion.
 *  \param e An error handler.
 *  \returns A new IMFExpansionNSpace.
 */
MFNSpace IMFCreateExpansionSpace(int n, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFCreateExpansionSpace"};
  MFNSpace thisSpace;
  struct IMFExpansionSpaceData *data;

  thisSpace=MFCreateNSpaceBaseClass("ExpansionSpace",e);
  MFNSpaceSetDistance(thisSpace,IMFExpansionSpaceDistance,e);
  MFNSpaceSetInnerProduct(thisSpace,IMFExpansionSpaceInner,e);
  MFNSpaceSetDirection(thisSpace,IMFExpansionSpaceDirection,e);
  MFNSpaceSetAdd(thisSpace,IMFExpansionSpaceAdd,e);
  MFNSpaceSetScale(thisSpace,IMFExpansionSpaceScale,e);
  MFNSpaceSetFreeData(thisSpace,MFFreeExpansionSpaceData,e);
  data=(struct IMFExpansionSpaceData*)malloc(sizeof(struct IMFExpansionSpaceData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(IMFExpansionSpaceErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct IMFExpansionSpaceData));
    MFSetError(e,12,RoutineName,IMFExpansionSpaceErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return thisSpace;
   }
#endif

  data->n=n;

  MFNSpaceSetData(thisSpace,(void*)data,e);

  return thisSpace;
 }

void MFFreeExpansionSpaceData(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeExpansionSpaceData"};
  struct IMFExpansionSpaceData *d;

  d=(struct IMFExpansionSpaceData*)data;
  if(d!=NULL)free(d);

  return;
 }

double IMFExpansionSpaceDistance(MFNSpace thisSpace,MFNVector v0,MFNVector v1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionSpaceDistance"};
  double result;
  int i;
  struct IMFExpansionSpaceData *data;

  data=(struct IMFExpansionSpaceData*)d;

  result=0.;
  for(i=0;i<data->n;i++)result+=pow(MFNV_C(v0,i,e)-MFNV_C(v1,i,e),2);
  result=sqrt(result);

  return result;
 }

void IMFExpansionSpaceDirection(MFNSpace thisSpace,MFNVector v0,MFNVector v1,MFNVector diff,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionSpaceDirection"};
  int i;
  struct IMFExpansionSpaceData *data;

  data=(struct IMFExpansionSpaceData*)d;

  for(i=0;i<data->n;i++)MFNVSetC(diff,i,MFNV_C(v1,i,e)-MFNV_C(v0,i,e),e);

  return;
 }

double IMFExpansionSpaceInner(MFNSpace thisSpace,MFNVector v0,MFNVector v1,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionSpaceInner"};
  double result;
  int i;
  int verbose=0;
  struct IMFExpansionSpaceData *data;

  data=(struct IMFExpansionSpaceData*)d;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s 0x%8.8x\n",RoutineName,thisSpace);
    printf(" v0 ");MFPrintNVector(stdout,v0,e);printf("\n");
    printf(" v1 ");MFPrintNVector(stdout,v1,e);printf("\n");
    printf(" n=%d\n",data->n);fflush(stdout);
   }
#endif

  result=0.;
  for(i=0;i<data->n;i++)
    result+=MFNV_C(v0,i,e)*MFNV_C(v1,i,e);

  return result;
 }

void IMFExpansionSpaceAdd(MFNSpace thisSpace,MFNVector v0,MFNVector v1,MFNVector sum,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionSpaceAdd"};
  int i;
  struct IMFExpansionSpaceData *data;

  data=(struct IMFExpansionSpaceData*)d;

  for(i=0;i<data->n;i++)MFNVSetC(sum,i,MFNV_C(v1,i,e)+MFNV_C(v0,i,e),e);

  return;
 }

void IMFExpansionSpaceScale(MFNSpace thisSpace,double s, MFNVector v,MFNVector prod,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFExpansionSpaceScale"};
  int i;
  struct IMFExpansionSpaceData *data;

  data=(struct IMFExpansionSpaceData*)d;

  for(i=0;i<data->n;i++)MFNVSetC(prod,i,s*MFNV_C(v,i,e),e);

  return;
 }


#ifdef __cplusplus
}
#endif
