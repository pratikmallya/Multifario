/* 
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   November 2, 1999
 */

static char *id="@(#) $Id: MFVolumeOfChart.c,v 1.5 2011/07/21 17:42:46 mhender Exp $";

static char MFVolumeOfChartErrorMsg[256]="";

#include <math.h>
#include <MFAtlas.h>
#include <MFKVector.h>
#include <MFNVector.h>
#include <MFChart.h>
#include <MFNRegion.h>
#include <MFErrorHandler.h>
#include <MFPrint.h>
#include <stdlib.h>
#include <string.h>
#include <MFEnumPolytope.h>
#include <MFEnumDualPolytope.h>
#include <MFFortran.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

int MFCheckHC(MFChart,MFNRegion,int,int,double*,double*,MFErrorHandler);
int MFCheckHCAllPolytopes(MFAtlas,int,int,double*,double*,MFErrorHandler);
void MFSubdivideHC(int,int,int*,int*,double***,double***,MFErrorHandler);
MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);
void Vol(MFEnumPolytope,int,int,int*,int,int*,double*,MFErrorHandler);
double MFVolumeOfPolytope(MFPolytope,MFErrorHandler);
double MFVolumeOfAllPolytopes(MFAtlas,MFErrorHandler);

static int Subset(int,int*,int,int*,MFErrorHandler);
static int Element(int,int,int*,MFErrorHandler);
int MFEnumPolytopeVerticesOnCell(MFEnumPolytope,int,int,int,int**,MFErrorHandler);
static int AddElement(int,int,int*,int**,MFErrorHandler);
double MFVolumeDet(MFEnumPolytope,int,int*,MFErrorHandler);
static int Fact(int,MFErrorHandler);
#include <MFFortran.h>


MFKVector MFVolumeTestPtK=NULL;
MFNVector MFVolumeTestPtN=NULL;
int      *MFVolumeIndex=NULL;

double MFVolumeOfChart(MFChart chart, MFNRegion Omega, MFErrorHandler e)
 {
  static char RoutineName[]={"MFVolumeOfChart"};
  int i,j;
  int n,k;
  int rc;
  int nHC,mHC;
  double **HCMin=NULL;
  double **HCMax=NULL;
  double volume,volHC;
  int nIncluded;
  double refVolHC,refVolBall;
  double testVolume;
  double R;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n","MFVolumeOfChart");fflush(stdout);}
#endif

  mHC=100;
  nHC=0;
  HCMin=(double**)malloc(mHC*sizeof(double*));

#ifndef MFNOSAFETYNET
  if(HCMin==NULL)
   {
    sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",mHC*sizeof(double));
    MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0.;
   }
#endif

  HCMax=(double**)malloc(mHC*sizeof(double*));

#ifndef MFNOSAFETYNET
  if(HCMax==NULL)
   {
    sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",mHC*sizeof(double));
    MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0.;
   }
#endif

  for(i=0;i<mHC;i++)
   {
    HCMin[i]=NULL;
    HCMax[i]=NULL;
   }

  k=MFChartK(chart,e);
  n=MFChartN(chart,e);

  if(MFVolumeTestPtK==NULL)MFVolumeTestPtK=MFCreateKVector(k,e);
  if(MFVolumeTestPtN==NULL)MFVolumeTestPtN=MFCreateNVector(n,e);

  if(MFVolumeIndex==NULL)
   {
    MFVolumeIndex=(int*)malloc(k*sizeof(int));

#ifndef MFNOSAFETYNET
    if(MFVolumeIndex==NULL)
     {
      sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(int));
      MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0.;
     }
#endif

   }
  
  volume=0.;
  testVolume=0.;
  nIncluded=0;
  HCMin[0]=(double*)(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(HCMin[0]==NULL)
   {
    sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0.;
   }
#endif

  HCMax[0]=(double*)(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(HCMax[0]==NULL)
   {
    sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0.;
   }
#endif

  R=MFChartRadius(chart,e);
  refVolHC=1.;
  for(i=0;i<k;i++)
   {
    (HCMin[0])[i]=-R;
    (HCMax[0])[i]=R;
    refVolHC=refVolHC*2*R;
   }
  nHC=1;

  if(k==1)
    refVolBall=2*R;
   else if(k==2)
    refVolBall=3.1415926*R*R;
   else if(k==3)
    refVolBall=4./3.*3.1415926*R*R*R;
   else if(k%2==0)
    {
     refVolBall=3.1415926*R*R;
     for(i=4;i<=k;i+=2)refVolBall=2*3.1415926*R*R/i;
    }else{
     refVolBall=4./3.*3.1415926*R*R*R;
     for(i=5;i<=k;i+=2)refVolBall=2*3.1415926*R*R/i;
    }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Initial List of Hypercubes %d\n",nHC);fflush(stdout);
    for(i=0;i<nHC;i++)
     {
      printf("   %d ",i);
      printf("(%lf",(HCMin[i])[0]);
      for(j=1;j<k;j++)printf(",%lf",(HCMin[i])[j]);
      printf(")<->");
      printf("(%lf",(HCMax[i])[0]);
      for(j=1;j<k;j++)printf(",%lf",(HCMax[i])[j]);
      printf(")\n");fflush(stdout);
     }
   }
#endif

  MFSubdivideHC(k,0,&nHC,&mHC,&HCMin,&HCMax,e);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Initial Subdivided List of Hypercubes %d\n",nHC);fflush(stdout);
    for(i=0;i<nHC;i++)
     {
      printf("   %d ",i);
      printf("(%lf",(HCMin[i])[0]);
      for(j=1;j<k;j++)printf(",%lf",(HCMin[i])[j]);
      printf(")<->");
      printf("(%lf",(HCMax[i])[0]);
      for(j=1;j<k;j++)printf(",%lf",(HCMax[i])[j]);
      printf(")\n");fflush(stdout);
     }
   }
#endif

  while(nHC>0)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("Checking HyperCube\n");fflush(stdout);}
#endif

    rc=MFCheckHC(chart,Omega,k,n,HCMin[0],HCMax[0],e);
    volHC=1.;
    for(i=0;i<k;i++)volHC*=(HCMax[0])[i]-(HCMin[0])[i];
    if(rc<0)
     {
      testVolume+=volHC;

#ifdef MFALLOWVERBOSE
      if(verbose){printf("All vertices are out, discarding\n");fflush(stdout);}
#endif

      free(HCMin[0]);
      free(HCMax[0]);
      HCMin[0]=NULL;
      HCMax[0]=NULL;
      if(nHC>1)
       {
        HCMin[0]=HCMin[nHC-1];
        HCMax[0]=HCMax[nHC-1];
       }
      nHC--;
     }else if(rc>0){
      testVolume+=volHC;

#ifdef MFALLOWVERBOSE
      if(verbose){printf("All vertices are in, volume is %lf\n",volHC);fflush(stdout);}
#endif

      volume+=volHC;
      nIncluded++;
      free(HCMin[0]);
      free(HCMax[0]);
      HCMin[0]=NULL;
      HCMax[0]=NULL;
      if(nHC>1)
       {
        HCMin[0]=HCMin[nHC-1];
        HCMax[0]=HCMax[nHC-1];
       }
      nHC--;
     }else{

#ifdef MFALLOWVERBOSE
      if(verbose){printf("Some vertices are in, some are out, volume is %lf\n",volHC);fflush(stdout);}
#endif

      if(volHC>1.e-10)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  Subdividing it.\n");fflush(stdout);}
#endif

        MFSubdivideHC(k,0,&nHC,&mHC,&HCMin,&HCMax,e);
       }else{
        testVolume+=volHC;

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  Too small, discarding it.\n");fflush(stdout);}
#endif

        free(HCMin[0]);
        free(HCMax[0]);
        HCMin[0]=NULL;
        HCMax[0]=NULL;
        if(nHC>1)
         {
          HCMin[0]=HCMin[nHC-1];
          HCMax[0]=HCMax[nHC-1];
         }
        nHC--;
       }
     }

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("List of Hypercubes %d\n",nHC);fflush(stdout);
      for(i=0;i<nHC;i++)
       {
        printf("   %d ",i);
        printf("(%lf",(HCMin[i])[0]);
        for(j=1;j<k;j++)printf(",%lf",(HCMin[i])[j]);
        printf(")<->");
        printf("(%lf",(HCMax[i])[0]);
        for(j=1;j<k;j++)printf(",%lf",(HCMax[i])[j]);
        printf(")\n");fflush(stdout);
       }
     }
#endif

   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("done %s\n","MFVolumeOfChart");
    printf("    volume of chart=%lf, nCubesInside=%d\n",volume,nIncluded);
    printf("    volume of ball =%lf,\n",refVolBall);
    printf("    volume of Cube =%lf\n",refVolHC);
    printf("    test volume    =%lf\n",testVolume);
    fflush(stdout);
   }
#endif

  return volume;
 }

int MFCheckHC(MFChart chart, MFNRegion Omega, int k, int n, double *HCMin, double *HCMax, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCheckHC"};
  int i;
  int rc;
  int allIn,allOut;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s k=%d, nHC=%d\n","MFCheckHC",k,n);fflush(stdout);
    printf("(%lf",HCMin[0]);
    for(i=1;i<k;i++)printf(",%lf",HCMin[i]);
    printf(")<->");
    printf("(%lf",HCMax[0]);
    for(i=1;i<k;i++)printf(",%lf",HCMax[i]);
    printf(")\n");fflush(stdout);
   }
#endif

  if(MFVolumeTestPtK==NULL)MFVolumeTestPtK=MFCreateKVector(k,e);
  if(MFVolumeTestPtN==NULL)MFVolumeTestPtN=MFCreateNVector(n,e);

  if(MFVolumeIndex==NULL)
   {
    MFVolumeIndex=(int*)malloc(k*sizeof(int));

#ifndef MFNOSAFETYNET
    if(MFVolumeIndex==NULL)
     {
      sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(int));
      MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

   }

  allIn=1;
  allOut=1;

  for(i=0;i<k;i++)MFVolumeIndex[i]=0;

  while(1)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf(" index (%d",MFVolumeIndex[k-1]);
      for(i=1;i<k;i++)printf(",%d",MFVolumeIndex[k-1-i]);
      printf(") s=");fflush(stdout);
     }
#endif

    for(i=0;i<k;i++)
     MFKVSetC(MFVolumeTestPtK,i,HCMin[i]+MFVolumeIndex[i]*(HCMax[i]-HCMin[i]),e);

#ifdef MFALLOWVERBOSE
    if(verbose)MFPrintKVector(stdout,MFVolumeTestPtK,e);
#endif

    rc=MFChartInterior(chart,MFVolumeTestPtK,e);
    if(rc)
     {

#ifdef MFALLOWVERBOSE
      if(verbose)printf(" In chart, u=");
#endif

      MFChartPointInTangentSpace(chart,MFVolumeTestPtK,MFVolumeTestPtN,e);

#ifdef MFALLOWVERBOSE
      if(verbose)MFPrintNVector(stdout,MFVolumeTestPtN,e);
#endif

      rc=MFNRegionInterior(Omega,MFVolumeTestPtN,e);

#ifdef MFALLOWVERBOSE
      if(verbose)
       {
        if(rc)printf(" In Omega.\n");
         else printf(" Not in Omega.\n");
       }
#endif

      if(rc)allOut=0;
       else allIn=0;
     }else{
       allIn=0;

#ifdef MFALLOWVERBOSE
       if(verbose)printf(" Not in chart.\n");
#endif
     }

    i=0;
    MFVolumeIndex[i]++;
    while(MFVolumeIndex[i]==2)
     {
      MFVolumeIndex[i]=0;
      i++;
      if(i>=k)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("done %s allIn=%d, allOut=%d\n","MFCheckHC",allIn,allOut);fflush(stdout);}
#endif

        if(allIn)return 1;
         else if(allOut)return -1;
         else return 0;
       }
      MFVolumeIndex[i]++;
     }
   }
 }

void MFSubdivideHC(int k,int iHC,int *nHC,int *mHC,double ***HCMin,double ***HCMax, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSubdivideHC"};
  int i;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, k=%d, iHC=%d\n","MFSubdivideHC",k, iHC);fflush(stdout);}
#endif

  if(MFVolumeIndex==NULL)
   {
    MFVolumeIndex=(int*)malloc(k*sizeof(int));

#ifndef MFNOSAFETYNET
    if(MFVolumeIndex==NULL)
     {
      sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(int));
      MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

   }

  for(i=0;i<k;i++)MFVolumeIndex[i]=0;

  while(1)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf(" index (%d",MFVolumeIndex[k-1]);
      for(i=1;i<k;i++)printf(",%d",MFVolumeIndex[k-1-i]);
      printf(")\n");fflush(stdout);
     }
#endif

    if(*nHC+1>*mHC)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("Resizing HCMin and HCMax\n");fflush(stdout);}
#endif

      (*mHC)+=100;
      *HCMin=(double**)realloc((void*)(*HCMin),(*mHC)*sizeof(double*));

#ifndef MFNOSAFETYNET
      if(*HCMin==NULL)
       {
        sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",(*mHC)*sizeof(double*));
        MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      *HCMax=(double**)realloc((void*)(*HCMax),(*mHC)*sizeof(double*));

#ifndef MFNOSAFETYNET
      if(*HCMax==NULL)
       {
        sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",(*mHC)*sizeof(double*));
        MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
        MFErrorHandlerOutOfMemory(e);
        return;
       }
#endif

      if(*HCMin==NULL||*HCMax==(double**)NULL)
       {
        printf("Out of memory in SubdivideHC\n");fflush(stdout);
        exit(12);
       }
      for(i=0;i<100;i++)
       {
        (*HCMin)[*mHC-101+i]=NULL;
        (*HCMax)[*mHC-101+i]=NULL;
       }
     }

    (*nHC)++;
    (*HCMin)[*nHC-1]=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
    if(*HCMax[*nHC-1]==NULL)
     {
      sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
      MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    (*HCMax)[*nHC-1]=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
    if(*HCMax[*nHC-1]==NULL)
     {
      sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
      MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    for(i=0;i<k;i++)
     {
      if(MFVolumeIndex[i]==0)
       {
        ((*HCMin)[*nHC-1])[i]=((*HCMin)[iHC])[i];
        ((*HCMax)[*nHC-1])[i]=((*HCMin)[iHC])[i]+.5*(((*HCMax)[iHC])[i]-((*HCMin)[iHC])[i]);
       }else{
        ((*HCMin)[*nHC-1])[i]=((*HCMin)[iHC])[i]+.5*(((*HCMax)[iHC])[i]-((*HCMin)[iHC])[i]);
        ((*HCMax)[*nHC-1])[i]=((*HCMax)[iHC])[i];
       }
     }

    i=0;
    MFVolumeIndex[i]++;
    while(MFVolumeIndex[i]==2)
     {
      MFVolumeIndex[i]=0;
      i++;
      if(i>=k)
       {
        if(*nHC>1)
         {
          free((*HCMin)[iHC]);
          free((*HCMax)[iHC]);
          (*HCMin)[0]=NULL;
          (*HCMax)[0]=NULL;
          if(*nHC>1)
           {
            (*HCMin)[iHC]=(*HCMin)[(*nHC)-1];
            (*HCMax)[iHC]=(*HCMax)[(*nHC)-1];
            (*HCMin)[(*nHC)-1]=NULL;
            (*HCMax)[(*nHC)-1]=NULL;
           }
          (*nHC)--;
         }

#ifdef MFALLOWVERBOSE
        if(verbose){printf("done %s %d\n","MFSubdivideHC",*nHC);fflush(stdout);}
#endif

        return;
       }
      MFVolumeIndex[i]++;
     }
   }
 }

double MFVolumeOfAtlas(MFAtlas A, MFNRegion Omega, MFErrorHandler e)
 {
  static char RoutineName[]={"MFVolumeOfAtlas"};
  int i,n;
  double volume;
  double volumeP;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n","MFVolumeOfAtlas");fflush(stdout);}
#endif

  volume=0.;
  volumeP=0.;
  n=MFAtlasNumberOfCharts(A,e);
  for(i=0;i<n;i++)
   {
/*  volume+=MFVolumeOfChart(MFAtlasChart(A,i,e),Omega,e);*/
    volumeP=MFVolumeOfPolytope(MFChartPolytope(MFAtlasChart(A,i,e),e),e);

#ifdef MFALLOWVERBOSE
    if(verbose)printf("  volume of polytope[%d]=%lf, radius of chart=%f\n",i,volumeP,MFChartRadius(MFAtlasChart(A,i,e),e));
#endif

    volume+=volumeP;
   }
/*volumeP=MFVolumeOfAllPolytopes(A,e);

#ifdef MFALLOWVERBOSE
  if(verbose)printf("  volume of all polytopes=%lf\n",volumeP);
#endif

*/

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n","MFVolumeOfAtlas");fflush(stdout);}
#endif

  return volume;
 }

double MFVolumeOfPolytope(MFPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFVolumeOfPolytope"};
  double volume;
  MFEnumPolytope Q;
  int nlast,*last;
  int nS,*S;
  int i;
  int d;
  int j,iv,nv,*v;
  MFKVector Pt;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("MFVolumeOfPolytope\n");fflush(stdout);}
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" Enumerate\n");fflush(stdout);}
#endif

  Q=MFEnumeratePolytope(P,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" done Enumerate\n");fflush(stdout);}
#endif

  volume=0;
  d=MFEnumPolytopeDimension(Q,e);
  nlast=MFEnumPolytopeNumberOfVertices(Q,e);
  last=(int*)malloc(nlast*sizeof(int));

#ifndef MFNOSAFETYNET
  if(last==NULL)
   {
    sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",nlast*sizeof(int));
    MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0.;
   }
#endif

  for(i=0;i<nlast;i++)last[i]=i;
  nS=0;
  S=NULL;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Enumerated Polytope:\n");

    printf("  vertices:\n");
    for(i=0;i<MFEnumPolytopeNumberOfVertices(Q,e);i++)
     {
      Pt=MFEnumPolytopeVertex(Q,i,e);
      printf("  %d [%lf",i,MFKV_C(Pt,0,e));
      for(j=1;j<d;j++)printf(",%lf",MFKV_C(Pt,j,e));
      printf("]\n");fflush(stdout);
     }
  
    for(i=1;i<=d;i++)
     {
      printf("  %d-cells:\n",i);
      for(j=0;j<MFEnumPolytopeNumberOfCells(Q,i,e);j++)
       {
        nv=0;
        v=NULL;
        nv=MFEnumPolytopeVerticesOnCell(Q,i,j,nv,&v,e);
        printf("  %d",j);
        if(nv>0)
         {
          printf(" [%d",v[0]);
          for(iv=1;iv<nv;iv++)
           {
            printf(",%d",v[iv]);
           }
          printf("]\n");fflush(stdout);
         }else{
          printf(" [empty set]\n");fflush(stdout);
         }
        free(v);
       }
     }

    printf(" calc Vol\n");fflush(stdout);
   }
#endif

  Vol(Q,d,nlast,last,nS,S,&volume,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" done Vol, vol=%lf/%d!=%lf/%d\n",volume,d,volume,Fact(d,e));fflush(stdout);}
#endif

  volume=volume/Fact(d,e);

  MFFreeEnumPolytope(Q,e);
  free(last);

  return volume;
 }

void Vol(MFEnumPolytope Q, int d, int nlast, int *last, int nS, int *S, double *volume, MFErrorHandler e)
 {
  static char RoutineName[]={"Vol"};
  int i;
  int n;
  int nPsi,*psi;
  int nSnew,*Snew=NULL;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf(" Vol, dimension %d",d);fflush(stdout);
    if(nlast>0)
     {
      printf(", last=[%d",last[0]);
      for(i=1;i<nlast;i++)printf(",%d",last[i]);
      printf("]");fflush(stdout);
     }else{
      printf(", last=empty");
     }
    if(nS>0)
     {
      printf(", S=[%d",S[0]);
      for(i=1;i<nS;i++)printf(",%d",S[i]);
      printf("]");fflush(stdout);
     }else{
      printf(", S=empty");
     }
    printf("\n");fflush(stdout);
   }
#endif

  n=MFEnumPolytopeNumberOfCells(Q,d,e);

  if(d>-1)
   {
    for(i=0;i<n;i++)
     {
      psi=NULL;
      nPsi=0;
      nPsi=MFEnumPolytopeVerticesOnCell(Q,d,i,nPsi,&psi,e);
      if(Subset(nPsi,psi,nlast,last,e) && !Element(psi[0],nS,S,e))
       {
        nSnew=AddElement(psi[0],nS,S,&Snew,e);
        Vol(Q,d-1,nPsi,psi,nSnew,Snew,volume,e);
        free(Snew);
       }
      free(psi);
     }
   }else{

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("Simplex: ");
      printf("[%d",S[0]);
      for(i=1;i<nS;i++)printf(",%d",S[i]);
      printf("]\n");fflush(stdout);
     }
#endif

    *volume+=fabs(MFVolumeDet(Q,nS,S,e));
   }
  return;
 }

int Subset(int n1,int *S1,int n2,int *S2, MFErrorHandler e)
 {
  static char RoutineName[]={"Subset"};
  static int i,j;
  static int result;
  static int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("    Subset, is");
    if(n1>0)
     {
      printf(" [%d",S1[0]);
      for(i=1;i<n1;i++)printf(",%d",S1[i]);
      printf("]");
     }else{
      printf(" the empty set");
     }
    printf(" a subset of ");
    if(n2>0)
     {
      printf(" [%d",S2[0]);
      for(i=1;i<n2;i++)printf(",%d",S2[i]);
      printf("]");
     }else{
      printf(" the empty set");
     }
    printf("?");fflush(stdout);
   }
#endif

  for(i=0;i<n1;i++)
   {
    result=0;
    for(j=0;j<n2;j++)if(S1[i]==S2[j])result=1;

#ifdef MFALLOWVERBOSE
    if(!result){if(verbose){printf(" No\n");fflush(stdout);}return 0;}
#endif

   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" Yes\n");fflush(stdout);}
#endif

  return 1;
 }

int Element(int e,int n,int *S, MFErrorHandler eee)
 {
  static char RoutineName[]={"Element"};
  static int i;
  static int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("    Element, is %d an element of ");
    if(n>0)
     {
      printf(" [%d",S[0]);
      for(i=1;i<n;i++)printf(",%d",S[i]);
      printf("]");
     }else{
      printf(" the empty set");
     }
    printf("?");fflush(stdout);
   }
#endif

  for(i=0;i<n;i++)

#ifdef MFALLOWVERBOSE
   if(e==S[i]){if(verbose){printf(" Yes.\n");fflush(stdout);}return 1;}
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" No.\n");fflush(stdout);}
#endif

  return 0;
 }

int AddElement(int e,int n,int *S,int **Snew, MFErrorHandler err)
 {
  static char RoutineName[]={"AddElement"};
  static int i;
  static int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("       AddElement %d to",err); 
    if(n>0)
     {
      printf(" [%d",S[0]);
      for(i=1;i<n;i++)printf(",%d",S[i]);
      printf("]");
     }else{
      printf(" [empty list]");
     }
   }
#endif
  
  *Snew=(int*)malloc((n+1)*sizeof(int));

#ifndef MFNOSAFETYNET
  if(*Snew==NULL)
   {
    sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",(n+1)*sizeof(int));
    MFSetError(err,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(err);
    return -1;
   }
#endif

  i=0;
  while(i<n && S[i]<e)
   {
    (*Snew)[i]=S[i];
    i++;
   }
  (*Snew)[i]=e;
  i++;
  while(i<n+1)
   {
    (*Snew)[i]=S[i-1];
    i++;
   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("=[%d",(*Snew)[0]);
    for(i=1;i<n+1;i++)printf(",%d",(*Snew)[i]);
    printf("]\n");fflush(stdout);
   }
#endif

  return n+1;
 }

double MFVolumeDet(MFEnumPolytope Q,int n,int *S, MFErrorHandler e)
 {
  static char RoutineName[]={"MFVolumeDet"};
  static int i,j;
  static int d;
  static double *A=NULL;
  static MFKVector row;
  static MFKVector lastRow;
  static int *ipivot=NULL;
  static int info;
  static double det;
  static int verbose=0;

  d=MFEnumPolytopeDimension(Q,e);

  A=(double*)realloc((void*)A,d*d*sizeof(double));

#ifndef MFNOSAFETYNET
  if(A==NULL)
   {
    sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",d*d*sizeof(double));
    MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0.;
   }
#endif

  lastRow=MFEnumPolytopeVertex(Q,S[d],e);
  if(n!=d+1)
   {
    printf(" problem in MFVolumeDet, n=%d, while dimension is %d\n",n,d);
    fflush(stdout);
    abort();
   }

#ifdef MFALLOWVERBOSE
  if(verbose)printf(" MFVolumeDet:\n");
#endif

  for(i=0;i<d;i++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)printf("    [");
#endif

    row=MFEnumPolytopeVertex(Q,S[i],e);
    for(j=0;j<d;j++)
     {
      A[i+d*j]=MFKV_C(row,j,e)-MFKV_C(lastRow,j,e);

#ifdef MFALLOWVERBOSE
      if(verbose)
       {
        if(j>0)printf(" ");
        printf("%lf",A[i+d*j]);
       }
#endif

     }

#ifdef MFALLOWVERBOSE
    if(verbose){printf("]\n");fflush(stdout);}
#endif

   }

  if(d==1)
   {
    det=A[0];
   }else if(d==2){
    det=A[0]*A[3]-A[1]*A[2];
   }else if(d==3){
    det=A[0]*(A[4]*A[8]-A[5]*A[7])
       -A[1]*(A[3]*A[8]-A[5]*A[6])
       +A[2]*(A[3]*A[7]-A[4]*A[6]);
   }else if(d==4){
    det=A[0]*( 
               A[5]*(A[10]*A[15]-A[11]*A[14])
              -A[6]*(A[ 9]*A[15]-A[11]*A[13])
              +A[7]*(A[ 9]*A[14]-A[10]*A[13])
             )
       -A[1]*(
               A[4]*(A[10]*A[15]-A[11]*A[14])
              -A[6]*(A[ 8]*A[15]-A[11]*A[12])
              +A[7]*(A[ 8]*A[14]-A[10]*A[12])
             )
       +A[2]*(
               A[4]*(A[ 9]*A[15]-A[11]*A[13])
              -A[5]*(A[ 8]*A[15]-A[11]*A[12])
              +A[7]*(A[ 8]*A[13]-A[ 9]*A[12])
             )
       -A[3]*(
               A[4]*(A[ 9]*A[14]-A[10]*A[13])
              -A[5]*(A[ 8]*A[14]-A[10]*A[12])
              +A[6]*(A[ 8]*A[13]-A[ 9]*A[12])
             );
   }else{
/*  ipivot=realloc((void*)ipivot,d*sizeof(int));

#ifndef MFNOSAFETYNET
    if(ipivot==NULL)
     {
      sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",d*sizeof(int));
      MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0.;
     }
#endif

    CALLDGETRF(&d,&d,A,&d,ipivot,&info);
    CALLDGETRS(&trans,&n,&one,A,&n,ipvt,b,&n,&info);*/
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("   MFVolumeDet: %lf\n",det);fflush(stdout);}
#endif
  return det;
 }

int Fact(int n, MFErrorHandler e)
 {
  static char RoutineName[]={"Fact"};

  if(n>2)return n*Fact(n-1,e);
   else return 2;
 }

double MFVolumeOfAllPolytopes(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFVolumeOfAllPolytopes"};
  int i,j;
  int n,k;
  int rc;
  int nHC,mHC;
  double **HCMin=NULL;
  double **HCMax=NULL;
  double volume,volHC;
  int nIncluded;
  MFChart chart;
  int c;
  double R;
  MFNVector u;
  double t;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n","MFVolumeOfAllPolytopes");fflush(stdout);}
#endif

  mHC=100;
  nHC=0;
  HCMin=(double**)malloc(mHC*sizeof(double*));

#ifndef MFNOSAFETYNET
  if(HCMin==NULL)
   {
    sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",mHC*sizeof(double*));
    MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0.;
   }
#endif

  HCMax=(double**)malloc(mHC*sizeof(double*));

#ifndef MFNOSAFETYNET
  if(HCMax==NULL)
   {
    sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",mHC*sizeof(double*));
    MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0.;
   }
#endif

  for(i=0;i<mHC;i++)
   {
    HCMin[i]=NULL;
    HCMax[i]=NULL;
   }

  k=MFAtlasK(A,e);
  n=MFAtlasN(A,e);

  if(k!=n)return -1.;

  if(MFVolumeTestPtK==NULL)MFVolumeTestPtK=MFCreateKVector(k,e);
  if(MFVolumeTestPtN==NULL)MFVolumeTestPtN=MFCreateNVector(n,e);
  if(MFVolumeIndex==NULL)
   {
    MFVolumeIndex=(int*)malloc(k*sizeof(int));

#ifndef MFNOSAFETYNET
    if(HCMax==NULL)
     {
      sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(int));
      MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0.;
     }
#endif

   }
  
  volume=0.;
  nIncluded=0;
  HCMin[0]=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(HCMin[0]==NULL)
   {
    sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0.;
   }
#endif

  HCMax[0]=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(HCMax[0]==NULL)
   {
    sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return 0.;
   }
#endif


  for(c=0;c<MFAtlasNumberOfCharts(A,e);c++)
   {
    chart=MFAtlasChart(A,c,e);
    R=MFChartRadius(chart,e);
    u=MFChartCenter(chart,e);
    for(i=0;i<k;i++)
     {
      if(c==0)
       {
        (HCMin[0])[i]=-1.05*R+MFNV_C(u,i,e);
        (HCMax[0])[i]=1.05*R+MFNV_C(u,i,e);
       }else{
        t=-1.05*R+MFNV_C(u,i,e);
        if((HCMin[0])[i]>t)(HCMin[0])[i]=t;
        t= 1.05*R+MFNV_C(u,i,e);
        if((HCMax[0])[i]<t)(HCMax[0])[i]=t;
       }
     }
   }

  nHC=1;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Initial List of Hypercubes %d\n",nHC);fflush(stdout);
    for(i=0;i<nHC;i++)
     {
      printf("   %d ",i);
      printf("(%lf",(HCMin[i])[0]);
      for(j=1;j<k;j++)printf(",%lf",(HCMin[i])[j]);
      printf(")<->");
      printf("(%lf",(HCMax[i])[0]);
      for(j=1;j<k;j++)printf(",%lf",(HCMax[i])[j]);
      printf(")\n");fflush(stdout);
     }
   }
#endif

  MFSubdivideHC(k,0,&nHC,&mHC,&HCMin,&HCMax,e);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Initial Subdivided List of Hypercubes %d\n",nHC);fflush(stdout);
    for(i=0;i<nHC;i++)
     {
      printf("   %d ",i);
      printf("(%lf",(HCMin[i])[0]);
      for(j=1;j<k;j++)printf(",%lf",(HCMin[i])[j]);
      printf(")<->");
      printf("(%lf",(HCMax[i])[0]);
      for(j=1;j<k;j++)printf(",%lf",(HCMax[i])[j]);
      printf(")\n");fflush(stdout);
     }
   }
#endif

  while(nHC>0)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("Checking HyperCube\n");fflush(stdout);}
#endif

    rc=MFCheckHCAllPolytopes(A,k,n,HCMin[0],HCMax[0],e);
    volHC=1.;
    for(i=0;i<k;i++)volHC*=(HCMax[0])[i]-(HCMin[0])[i];
    if(rc<0)
     {

#ifdef MFALLOWVERBOSE
      if(verbose){printf("All vertices are out, discarding\n");fflush(stdout);}
#endif

      free(HCMin[0]);
      free(HCMax[0]);
      HCMin[0]=NULL;
      HCMax[0]=NULL;
      if(nHC>1)
       {
        HCMin[0]=HCMin[nHC-1];
        HCMax[0]=HCMax[nHC-1];
       }
      nHC--;
     }else if(rc>0){

#ifdef MFALLOWVERBOSE
      if(verbose){printf("All vertices are in, volume is %lf\n",volHC);fflush(stdout);}
#endif

      volume+=volHC;
      nIncluded++;
      free(HCMin[0]);
      free(HCMax[0]);
      HCMin[0]=NULL;
      HCMax[0]=NULL;
      if(nHC>1)
       {
        HCMin[0]=HCMin[nHC-1];
        HCMax[0]=HCMax[nHC-1];
       }
      nHC--;
     }else{

#ifdef MFALLOWVERBOSE
      if(verbose){printf("Some vertices are in, some are out, volume is %lf\n",volHC);fflush(stdout);}
#endif

      if(volHC>1.e-10)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  Subdividing it.\n");fflush(stdout);}
#endif

        MFSubdivideHC(k,0,&nHC,&mHC,&HCMin,&HCMax,e);
       }else{

#ifdef MFALLOWVERBOSE
        if(verbose){printf("  Too small, discarding it.\n");fflush(stdout);}
#endif

        free(HCMin[0]);
        free(HCMax[0]);
        HCMin[0]=NULL;
        HCMax[0]=NULL;
        if(nHC>1)
         {
          HCMin[0]=HCMin[nHC-1];
          HCMax[0]=HCMax[nHC-1];
         }
        nHC--;
       }
     }

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf("List of Hypercubes %d\n",nHC);fflush(stdout);
      for(i=0;i<nHC;i++)
       {
        printf("   %d ",i);
        printf("(%lf",(HCMin[i])[0]);
        for(j=1;j<k;j++)printf(",%lf",(HCMin[i])[j]);
        printf(")<->");
        printf("(%lf",(HCMax[i])[0]);
        for(j=1;j<k;j++)printf(",%lf",(HCMax[i])[j]);
        printf(")\n");fflush(stdout);
       }
     }
#endif

   }

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    if(verbose)printf("done %s\n","MFVolumeOfAllPolytopes");
    printf("    volume of atlas=%lf, nCubesInside=%d\n",volume,nIncluded);
    fflush(stdout);
   }
#endif

  return volume;
 }

int MFCheckHCAllPolytopes(MFAtlas A,int k,int n,double *HCMin,double *HCMax, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCheckHCAllPolytopes"};
  int i;
  int rc;
  int allIn,allOut;
  int verbose=0;
  int c;
  MFChart chart;
  double R;
  MFNVector u;
  int in;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s\n","MFCheckHCAllPolytopes");fflush(stdout);
    printf("(%lf",HCMin[0]);
    for(i=1;i<k;i++)printf(",%lf",HCMin[i]);
    printf(")<->");
    printf("(%lf",HCMax[0]);
    for(i=1;i<k;i++)printf(",%lf",HCMax[i]);
    printf(")\n");fflush(stdout);
   }
#endif

  if(MFVolumeTestPtK==NULL)MFVolumeTestPtK=MFCreateKVector(k,e);
  if(MFVolumeTestPtN==NULL)MFVolumeTestPtN=MFCreateNVector(n,e);
  if(MFVolumeIndex==NULL)
   {
    MFVolumeIndex=(int*)malloc(k*sizeof(int));

#ifndef MFNOSAFETYNET
    if(MFVolumeIndex==NULL)
     {
      sprintf(MFVolumeOfChartErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(int));
      MFSetError(e,12,RoutineName,MFVolumeOfChartErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

   }

  allIn=1;
  allOut=1;

  for(i=0;i<k;i++)MFVolumeIndex[i]=0;

  while(1)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)
     {
      printf(" index (%d",MFVolumeIndex[k-1]);
      for(i=1;i<k;i++)printf(",%d",MFVolumeIndex[k-1-i]);
      printf(") s=");fflush(stdout);
      i=0;
      printf("[%lf",HCMin[i]+MFVolumeIndex[i]*(HCMax[i]-HCMin[i]));
      for(i=1;i<k;i++)printf(",%lf",HCMin[i]+MFVolumeIndex[i]*(HCMax[i]-HCMin[i]));
      printf("]");fflush(stdout);
     }
#endif

    rc=0;
    for(c=0;c<MFAtlasNumberOfCharts(A,e);c++)
     {
      chart=MFAtlasChart(A,c,e);
      R=MFChartRadius(chart,e);
      u=MFChartCenter(chart,e);
      in=1;
      for(i=0;i<k;i++)
       {
        if(HCMin[i]+MFVolumeIndex[i]*(HCMax[i]-HCMin[i])-MFNV_C(u,i,e)<-1.05*R)in=0;
        if(HCMin[i]+MFVolumeIndex[i]*(HCMax[i]-HCMin[i])-MFNV_C(u,i,e)> 1.05*R)in=0;
       }
      if(in)rc=1;
     }

    if(rc)
     {

#ifdef MFALLOWVERBOSE
      if(verbose)printf(" In a Chart.\n");
#endif

      if(rc)allOut=0;
       else allIn=0;
     }else{
       allIn=0;

#ifdef MFALLOWVERBOSE
       if(verbose)printf(" Not in a chart.\n");
#endif

     }

    i=0;
    MFVolumeIndex[i]++;
    while(MFVolumeIndex[i]==2)
     {
      MFVolumeIndex[i]=0;
      i++;
      if(i>=k)
       {

#ifdef MFALLOWVERBOSE
        if(verbose){printf("done %s allIn=%d, allOut=%d\n","MFCheckHCAllPolytopes",allIn,allOut);fflush(stdout);}
#endif

        if(allIn)return 1;
         else if(allOut)return -1;
         else return 0;
       }
      MFVolumeIndex[i]++;
     }
   }
 }

#ifdef __cplusplus
}
#endif
