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

static char *id="@(#) $Id: MFComputePointsOn3dPolygon.c,v 1.4 2011/07/21 17:42:46 mhender Exp $";

static char MFComputePointsOn3dPolygonErrorMsg[256]="";

#include <MFAtlas.h>
#include <MFAtlasFriends.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <MFEnumDualPolytope.h>
#include <math.h>
#include <sh.h>
#include <MFDraw.h>
#include <MFPrint.h>

#ifdef __cplusplus
 extern "C" {
#endif

int MFComputePointsOn3dPolygon(double r,int nV,double *V,int nPt0,double *Pt0,double *R0,double **Pt,double **R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFComputePointsOn3dPolygon"};
  MFImplicitMF M;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0;
  MFNVector u;
  MFNKMatrix Phi;
  MFChart C;
  double Rad;
  int n;
  int i,j;
  double delta=0;
  int chart;
  int kkk;
  MFPolytope P;
  double rho,Ri,Rj;
  MFEnumDualPolytope Q;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("---begin-------------------------------------------------------------------\n");
    printf("MFComputePointsOn3dPolygon(r=%lf,nV=%d,V,nPt0=%d,Pt0,Pt);\n",r,nV,nPt0);
   }
#endif

  Omega=MFNRegionCreatePolygonal3dRegion(nV,V,e);

  M=MFIMFCreatePolygonIn3SpaceWithRadius(nV,V,r,e);

  S=MFCreateAtlas(M,e);

  kkk=0;
  for(i=0;i<nPt0;i++)
   {
    u0=MFCreateNVector(3,e);
    MFNVSetC(u0,0,Pt0[3*i  ],e);
    MFNVSetC(u0,1,Pt0[3*i+1],e);
    MFNVSetC(u0,2,Pt0[3*i+2],e);
    Phi=MFIMFTangentSpace(M,u0,e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("%d) Add initial point to Atlas ",MFAtlasNumberOfCharts(S,e)+1);MFPrintNVector(stdout,u0,e);printf("\n");fflush(stdout);}
#endif

    C=MFCreateChartWithCubeSize(M,u0,Phi,R0[i],10.,e);
    MFAtlasAddChart(S,C,e);

/*  MFAtlasAddChartWithAll(S,u0,Phi,R0[i],e);*/

    kkk++;
    MFFreeNVector(u0,e);
    MFFreeNKMatrix(Phi,e);
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Done adding inital points\n");fflush(stdout);}
#endif

  u=MFCreateNVector(3,e);
  while((chart=MFAtlasPointOnBoundaryInsideRegion(S,Omega,u,&Phi,&delta,e))>-1)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("%d) Adding point ",MFAtlasNumberOfCharts(S,e)+1);MFPrintNVector(stdout,MFAtlasCenterOfChart(S,chart,e),e);printf("\n");fflush(stdout);}
#endif

    kkk++;
    Rad=MFIMFScale(MFAtlasMF(S,e),u,Phi,e);

    C=MFCreateChartWithCubeSize(M,u,Phi,Rad,10.,e);
    MFAtlasAddChart(S,C,e);
 
/*  MFAtlasAddChartWithAll(S,u,Phi,Rad,e); */
    MFFreeNVector(u,e);
    MFFreeNKMatrix(Phi,e);
/*  MFFreeChart(C,e);*/
    u=MFCreateNVector(3,e);
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("No more points on boundary\n");fflush(stdout);}
#endif

/*MFDrawClear(e);
  MFDrawAtlasOnce(S,e);
  Q=MFEnumDualOfAtlas(S,e);
  MFDrawEnumDualPolytope(Q,e);
  MFFreeEnumDualPolytope(Q,e);*/

/*DrawEdges(,e);
  shSetOutputFormat("tiff",e);
  shSetOutputFilename("Filled",e);
  MFDrawDisplay();MFDrawClear(,e);*/

  MFFreeNVector(u,e);

  n=MFAtlasNumberOfCharts(S,e)-nPt0;
  if(n>0)
   {
    (*Pt)=(double*)realloc((void*)(*Pt),3*n*sizeof(double));

#ifndef MFNOSAFETYNET
    if(*Pt==NULL)
     {
      sprintf(MFComputePointsOn3dPolygonErrorMsg,"Out of memory, trying to allocate %d bytes",3*n*sizeof(double));
      MFSetError(e,12,RoutineName,MFComputePointsOn3dPolygonErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    (*R)=(double*)realloc((void*)(*R),n*sizeof(double));

#ifndef MFNOSAFETYNET
    if(*R==NULL)
     {
      sprintf(MFComputePointsOn3dPolygonErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
      MFSetError(e,12,RoutineName,MFComputePointsOn3dPolygonErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return 0;
     }
#endif

    for(i=0;i<n;i++)
     {
      u=MFAtlasCenterOfChart(S,i+nPt0,e);
      (*Pt)[3*i  ]=MFNV_C(u,0,e);
      (*Pt)[3*i+1]=MFNV_C(u,1,e);
      (*Pt)[3*i+2]=MFNV_C(u,2,e);
      (*R)[i]=0.;
      P=MFChartPolytope(MFAtlasChart(S,i+nPt0,e),e);
      Ri=MFAtlasChartRadius(S,i,e);
      for(j=0;j<MFPolytopeNumberOfVertices(P,e);j++)
       {
        Rj=MFPolytopeRadiusOfVertex(P,j,e);
        rho=sqrt(2*Rj*Rj-Ri*Ri);
        if(rho>(*R)[i])(*R)[i]=rho;
       }
     }
   }

  MFFreeAtlas(S,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);

#ifdef MFALLOWVERBOSE
  if(verbose)printf("----end--------------------------------------------------------------------\n");
#endif

  return n;
 }

#ifdef __cplusplus
}
#endif
