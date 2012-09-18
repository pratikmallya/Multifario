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

static char *id="@(#) $Id: MFPendulaNRegion.c,v 1.3 2007/02/13 01:22:34 mhender Exp $";

static char MFNRegionErrorMsg[256]="";

#include <MFNVector.h>
#include <MFNRegion.h>
#include <MFErrorHandler.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
 extern "C" {
#endif

static int PendulaTest(MFNVector,void*,MFErrorHandler);
static void PendulaFree(void*,MFErrorHandler);
static int PendulaFourTest(MFNVector,void*,MFErrorHandler);
static void PendulaFourFree(void*,MFErrorHandler);

MFNRegion MFNRegionCreatePendula(int n, double kappa,double gamma,int WN,double IcMin,double IcMax, double InMin,double InMax,double TMax, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreatePendula"};
  MFNRegion pendula;
  double *data;

  pendula=MFNRegionCreateBaseClass("Pendula",e);
  MFNRegionSetTest(pendula,PendulaTest,e);
  MFNRegionSetFreeData(pendula,PendulaFree,e);

  data=(double*)malloc(6*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(double));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(pendula);
    return NULL;
   }
#endif

  data[0]=kappa;
  data[1]=IcMin;
  data[2]=IcMax;
  data[3]=InMin;
  data[4]=InMax;
  data[5]=TMax;

  MFNRegionSetData(pendula,data,e);

  return(pendula);
 }

int PendulaTest(MFNVector v,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"PendulaTest"};
  double kappa,T,In,Ic;
  double IcMin,IcMax,InMin,InMax,TMax;
  int i;
  int rc;
  int n;
  double *x;

  n=MFNV_NC(v,e);
  x=MFNV_CStar(v,e);

  kappa=((double*)data)[0];
  IcMin=((double*)data)[1];   /* -1. */
  IcMax=((double*)data)[2];   /*  1. */
  InMin=((double*)data)[3];   /*  0. */
  InMax=((double*)data)[4];   /* 10. */
  TMax=((double*)data)[5];    /* 50. */

  T=10*x[n-3];
  In=x[n-2];
  Ic=x[n-1];  /* Ic=3.1415926*kappa*x[n-1] */

  rc=1;
  if(Ic<IcMin)rc=0;
  if(Ic>IcMax)rc=0;
  if(In<InMin)rc=0;
  if(In>InMax)rc=0;
  if(T>TMax)rc=0;

  if(0){printf(" Pendula Test (Ic,In,T)=(%lf,%lf,%lf)  Ic in [%lf,%lf], In in [%lf,%lf], T<%lf, rc=%d\n",Ic,In,T,IcMin,IcMax,InMin,InMax,TMax,rc);fflush(stdout);}

  return rc;
 }

void PendulaFree(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"PendulaFree"};

  free(data);
  return;
 }

MFNRegion MFNRegionCreatePendulaFour(int n, int WN, double kappaMin,double kappaMax, double gammaMin,double gammaMax,double IcMin,double IcMax, double InMin,double InMax,double TMax, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreatePendulaFour"};
  MFNRegion pendula;
  double *data;

  pendula=MFNRegionCreateBaseClass("PendulaFour",e);
  MFNRegionSetTest(pendula,PendulaFourTest,e);
  MFNRegionSetFreeData(pendula,PendulaFourFree,e);

  data=(double*)malloc(8*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(double));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(pendula);
    return NULL;
   }
#endif

  data[0]=kappaMin;
  data[1]=IcMin;
  data[2]=IcMax;
  data[3]=InMin;
  data[4]=InMax;
  data[5]=TMax;
  data[6]=kappaMax;
  data[7]=gammaMin;
  data[8]=gammaMax;

  MFNRegionSetData(pendula,data,e);

/*printf("MFNRegion - PendulaFour, n=%d Winding Number=%d\n",n,WN);
  printf("                         kappa in %lf,%lf\n",kappaMin,kappaMax);
  printf("                         gamma in %lf,%lf\n",gammaMin,gammaMax);
  printf("                         Ic    in %lf,%lf\n",IcMin,IcMax);
  printf("                         In    in %lf,%lf\n",InMin,InMax);
  printf("                         T     < %lf\n",TMax);*/
  MFNRegionSetFreeData(pendula,PendulaFourFree,e);

  return(pendula);
 }

int PendulaFourTest(MFNVector v,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"PendulaFourTest"};
  double kappa,gamma,T,In,Ic;
  double IcMin,IcMax,InMin,InMax,TMax;
  double kappaMin,kappaMax,gammaMin,gammaMax;
  int i;
  int rc;
  int n;
  double *x;

  n=MFNV_NC(v,e);
  x=MFNV_CStar(v,e);

  kappaMin=((double*)data)[0];
  IcMin=((double*)data)[1];   /* -1. */
  IcMax=((double*)data)[2];   /*  1. */
  InMin=((double*)data)[3];   /*  0. */
  InMax=((double*)data)[4];   /* 10. */
  TMax=((double*)data)[5];    /* .1  */
  kappaMax=((double*)data)[6];
  gammaMin=((double*)data)[7];
  gammaMax=((double*)data)[8];

  T=x[n-5];
  In=x[n-4];
  Ic=x[n-3];  /* Ic=3.1415926*kappa*x[n-1] */
  kappa=x[n-2];
  gamma=x[n-1];

/*printf("%s\n",RoutineName);
  printf(" %lf < %lf (Ic) < %lf \n",IcMin,Ic,IcMax);
  printf(" %lf < %lf (In) < %lf \n",InMin,In,InMax);
  printf(" %lf > %lf (T) \n",TMax,T);
  printf(" %lf < %lf (kappa) < %lf \n",kappaMin,kappa,kappaMax);
  printf(" %lf < %lf (gamma) < %lf \n",gammaMin,gamma,gammaMax);*/

  rc=1;
  if(Ic<IcMin)rc=0;
  if(Ic>IcMax)rc=0;
  if(In<InMin)rc=0;
  if(In>InMax)rc=0;
  if(T>TMax)rc=0;
  if(kappa<kappaMin)rc=0;
  if(kappa>kappaMax)rc=0;
  if(gamma<gammaMin)rc=0;
  if(gamma>gammaMax)rc=0;

  return rc;
 }

void PendulaFourFree(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"PendulaFourFree"};

  free(data);
  return;
 }

#ifdef __cplusplus
}
#endif
