/* 
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   February 22, 1999
 */

static char *id="@(#) $Id: MFForcedNRegion.c,v 1.4 2007/04/04 17:51:20 mhender Exp $";

/*      author: Mike Henderson mhender@watson.ibm.com
 *      date:   November 11, 1997
 *              February 2, 1999   Ported to C
 *              October 6, 2000   Added Polyhedron
 */

#include <stdio.h>
#include <MFNVector.h>
#include <MFNRegion.h>
#include <stdlib.h>
#include <math.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

static char MFForcedNRegionMFErrorHandlerMsg[256]="";

int ForcedOscillatorTest(MFNVector,void*,MFErrorHandler);
void ForcedOscillatorFree(void*,MFErrorHandler);

double *MFNV_CStar(MFNVector,MFErrorHandler);

MFNRegion MFNRegionCreateForcedOscillator(int n, int WN,double eMin,double eMax, double BMin,double BMax, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreateForcedOscillator"};
  MFNRegion Omega;
  double *data;

  Omega=MFNRegionCreateBaseClass("ForcedOscillator",e);
  MFNRegionSetTest(Omega,ForcedOscillatorTest,e);
  MFNRegionSetFreeData(Omega,ForcedOscillatorFree,e);

  data=(double*)malloc(4*sizeof(double));
#ifndef MFNPOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFForcedNRegionMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",4*sizeof(double));
    MFSetError(e,12,RoutineName,MFForcedNRegionMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(Omega);
    return NULL;
   }
#endif

  data[0]=eMin;
  data[1]=eMax;
  data[2]=BMax;
  data[3]=BMin;

  MFNRegionSetData(Omega,data,e);

  return(Omega);
 }

int ForcedOscillatorTest(MFNVector v,void *data, MFErrorHandler eh)
 {
  static char RoutineName[]={"ForcedOscillatorTest"};
  double e,B;
  double eMin,eMax,BMin,BMax;
  int i;
  int rc;
  int n;
  double *x;

  n=MFNV_NC(v,eh);
  x=MFNV_CStar(v,eh);

  eMin=((double*)data)[0];
  eMax=((double*)data)[1];
  BMax=((double*)data)[2];
  BMin=((double*)data)[3];

  e=x[n-2];
  B=x[n-1];

  rc=1;
  if(e<eMin)rc=0;
  if(e>eMax)rc=0;
  if(B<BMin)rc=0;
  if(B>BMax)rc=0;

  return rc;
 }

void ForcedOscillatorFree(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"ForcedOscillatorFree"};

  free(data);
  return;
 }

#ifdef __cplusplus
}
#endif
