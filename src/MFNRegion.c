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
 */

static char *id="@(#) $Id: MFNRegion.c,v 1.3 2007/02/13 01:22:33 mhender Exp $";

static char MFNRegionErrorMsg[256]="";

#include <MFErrorHandler.h>
#include <stdio.h>
#include <MFNVector.h>
#include <MFNRegion.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct MFNRegionSt
 {
  char *id;
  int (*test)(MFNVector,void*,MFErrorHandler);
  void *data;
  void (*freedata)(void*,MFErrorHandler);
  void (*writedata)(FILE*,void*,MFErrorHandler);
 };

int MFNRegionInterior(MFNRegion Omega,MFNVector v, MFErrorHandler e)
 {
  static char RoutineName[]={"MFRegionInterior"};
  int result;

  result=Omega->test(v,Omega->data,e);

  return result;
 }

void MFFreeNRegion(MFNRegion Omega, MFErrorHandler e)
 {
  char RoutineName[]={"MFFreeNRegion"};

  if(Omega->id!=NULL)free(Omega->id);
  if(Omega->freedata!=NULL)Omega->freedata(Omega->data,e);
  free(Omega);

  return;
 }

MFNRegion MFNRegionCreateBaseClass(char *id, MFErrorHandler e)
 {
  MFNRegion result;
  static char RoutineName[]={"MFNRegionCreateBaseClass"};

  result=(MFNRegion)malloc(sizeof(struct MFNRegionSt));

#ifndef MFNOSAFETYNET
  if(result==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFNRegionSt));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  result->id=(char*)malloc((strlen(id)+1)*sizeof(char));

#ifndef MFNOSAFETYNET
  if(result->id==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",(strlen(id)+1)*sizeof(char));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  strcpy(result->id,id);

  return result;
 }

void MFNRegionSetTest(MFNRegion Omega,int (*test)(MFNVector,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionSetTest"};

#ifdef MFNOCONFIDENCE
  if(Omega==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Region (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  Omega->test=test;

  return;
 }

void MFNRegionSetData(MFNRegion Omega,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionSetData"};

#ifdef MFNOCONFIDENCE
  if(Omega==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Region (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  Omega->data=data;

  return;
 }

void *MFNRegionGetData(MFNRegion Omega, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionGetData"};

#ifdef MFNOCONFIDENCE
  if(Omega==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Region (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    return NULL;
   }
#endif

  return Omega->data;
 }

void MFNRegionSetFreeData(MFNRegion Omega,void (*freedata)(void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionSetFreeData"};

#ifdef MFNOCONFIDENCE
  if(Omega==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Region (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  Omega->freedata=freedata;
  return;
 }

void MFNRegionSetWriteData(MFNRegion Omega,void (*writedata)(FILE*,void*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionSetWriteData"};

#ifdef MFNOCONFIDENCE
  if(Omega==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Region (argument 1) is NULL");
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  Omega->writedata=writedata;
  return;
 }

#ifdef __cplusplus
}
#endif
