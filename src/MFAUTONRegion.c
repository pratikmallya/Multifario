/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   July 15, 2002 Copied from MFLOCANRegion.c
 */

static char *id="@(#) $Id: MFAUTONRegion.c,v 1.6 2008/05/02 12:44:17 mhender Exp $";

#include <multifarioConfig.h>

#ifdef HAVE_AUTO

static char MFAUTONRegionErrorMsg[256]="";

#include <stdio.h>
#include <MFNRegion.h>
#include <MFNSpace.h>
#include <MFErrorHandler.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
 extern "C" {
#endif

static int MFAUTONSpaceTest(MFNVector,void*,MFErrorHandler);
static void MFAUTONRegionFreeData(void*,MFErrorHandler);

double *MFAUTOBVNVGetPar(MFNVector,MFErrorHandler);

struct MFAUTONRegionData
 {
  MFNSpace space;
  int k;
  long *icp;
  double *RL0;
  double *RL1;
  double A0;
  double A1;
 };

MFNRegion MFNRegionCreateAUTO(MFNSpace space, int k, long *icp, double *RL0, double *RL1, double A0,double A1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreateAUTO"};
  MFNRegion Auto;
  struct MFAUTONRegionData *data;
  int i;

  Auto=MFNRegionCreateBaseClass("AUTO",e);
  MFNRegionSetTest(Auto,MFAUTONSpaceTest,e);
  data=malloc(sizeof(struct MFAUTONRegionData));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFAUTONRegionErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFAUTONRegionData));
    MFSetError(e,12,RoutineName,MFAUTONRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  MFNRegionSetData(Auto,data,e);
  MFNRegionSetFreeData(Auto,MFAUTONRegionFreeData,e);
  MFRefNSpace(space,e);

  data->space=space;
  data->k=k;
  data->icp=icp;
  data->RL0=malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->RL0==NULL)
   {
    sprintf(MFAUTONRegionErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTONRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<k;i++)(data->RL0)[i]=RL0[i];

  data->RL1=malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data->RL1==NULL)
   {
    sprintf(MFAUTONRegionErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,MFAUTONRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<k;i++)(data->RL1)[i]=RL1[i];

  data->A0=A0;
  data->A1=A1;

  return(Auto);
 }

int MFAUTONSpaceTest(MFNVector u,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFAUTONSpaceTest"};
  double *par;
  struct MFAUTONRegionData *data;
  double norm;
  int i;
  int result;
  int verbose=0;

  data=(struct MFAUTONRegionData*)d;

  result=1;
  par=MFAUTOBVNVGetPar(u,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  for(i=0;i<data->k;i++)
   {
#ifdef MFALLOWVERBOSE
    if(verbose){printf(" l%2d %lf <= %lf <= %lf\n",i,data->RL0[i],par[i],data->RL1[i]);fflush(stdout);}
#endif

    if( par[i]<data->RL0[i] || par[i]>data->RL1[i])result=0;
   }
  norm=sqrt(MFNSpaceInner(data->space,u,u,e));

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" |u|  %lf <= %lf <= %lf\n",data->A0,norm,data->A1);fflush(stdout);}
#endif

  if( norm<data->A0 || norm>data->A1)result=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" result = %d\n",result);fflush(stdout);}
#endif

  return result;
 }

void MFAUTONSpaceFreeData(void *d, MFErrorHandler e)
 {
  double *par;
  struct MFAUTONRegionData *data;

  data=(struct MFAUTONRegionData*)d;

  MFFreeNSpace(data->space,e);
  free(data->RL0);
  free(data->RL1);
  free(data);
  return;
 }

void MFAUTONRegionFreeData(void *d, MFErrorHandler e)
 {
  struct MFAUTONRegionData *data;

  data=(struct MFAUTONRegionData*)d;
 
  MFFreeNSpace(data->space,e);
  free(data->RL0);
  free(data->RL1);
  free(data);
 }

#else

int MFThereIsNoAUTO_MFAUTONRegion()
 {
  return 0;
 }

#endif

#ifdef __cplusplus
}
#endif
