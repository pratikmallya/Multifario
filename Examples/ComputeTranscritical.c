/* 
    %W%
    %D% %T%
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <MFAtlas.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <math.h>
#include <MFMultifariosMethod.h>

int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  int n;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0[2];
  FILE *fid;
  MFContinuationMethod H;
  MFErrorHandler e;
  MFNVector ur,ll;

  e=MFCreateErrorHandler();

/*M=MFIMFCreateAlgebraicExpression("[x,y,z]","[(z-x)*(x+y+z)]",e);*/
  M=MFIMFCreateAlgebraicExpression("[x,y,z]","[(z-x-.5*(y-1)*(y-1))*(x+y+z)]",e);
  n=MFIMF_N(M,e);

  ur=MFCreateNVector(n,e);
  MFNVSetC(ur,0,0.6,e);
  MFNVSetC(ur,1,0.6,e);
  MFNVSetC(ur,2,1.6,e);

  ll=MFCreateNVector(n,e);
  MFNVSetC(ll,0,-0.6,e);
  MFNVSetC(ll,1,-0.6,e);
  MFNVSetC(ll,2,-0.6,e);

/*Omega=MFNRegionCreateHyperCubeByCorners(n,ll,ur,e);*/
  Omega=MFNRegionCreateHyperCube(n,1.1,e);

  u0[0]=MFIMFVectorFactory(M,e);
  MFNVSetC(u0[0],0, 0.,e);
  MFNVSetC(u0[0],1, 1.,e);
  MFNVSetC(u0[0],2,-1.,e);

  u0[1]=MFIMFVectorFactory(M,e);
  MFNVSetC(u0[1],0, 1.,e);
  MFNVSetC(u0[1],1, 1.,e);
  MFNVSetC(u0[1],2, 1.,e);

  MFNVSetC(u0[0],0,0.0,e);
  MFNVSetC(u0[0],1, .1,e);
  MFNVSetC(u0[0],2,-.1,e);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.1,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",8000,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",1,e);
  MFMultifarioSetIntegerParameter(H,"pageEvery",100000,e);

  MFMultifarioSetRealParameter(H,"dotmin",.3,e);

  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  MFMultifarioSetIntegerParameter(H,"branchSwitch",2,e);
  MFMultifarioSetFilename(H,"Transcritical",e);

  S=MFComputeAtlasMultiple(H,M,Omega,1,u0,e);


  if(0){
   int i;
   int j;
   int n;
   int m;
   int vMark;
   int prt;

   MFChart C;
   MFPolytope P;

   n=MFAtlasNumberOfCharts(S,e);
   printf("\n\n\nThere are %d charts\n",n);fflush(stdout);
   for(i=0;i<n;i++)
    {
     C=MFAtlasChart(S,i,e);
     P=MFChartPolytope(C,e);
     m=MFPolytopeNumberOfVertices(P,e);

#if 0
     {
      double x,y,z;
      double v0,v1;
      char sign0,sign1;

      MFNVector u;

      u=MFChartCenter(C,e);
      x=MFNV_C(u,0,e);
      y=MFNV_C(u,1,e);
      z=MFNV_C(u,2,e);

      v0=z-x-.5*(y-1)*(y-1);
      sign0='+';if(v0<0.)sign0='-';if(fabs(v0)<1.e-7)sign0='0';

      v1=x+y+z;
      sign1='+';if(v1<1.)sign1='-';if(fabs(v1)<1.e-7)sign1='0';

      printf("  chart %d, index1=%d, index2=%d,  sign(v0)=%c, sign(v1)=%c\n",i,MFNVGetIndex(MFChartCenter(C,e),e),MFNVGetIndex2(MFChartCenter(C,e),e),sign0,sign1);fflush(stdout);
     }
#endif

     if(MFNVGetIndex2(MFChartCenter(C,e),e)>-1)
      {
       prt=0;
       for(j=0;j<m;j++)
        {
         double x,y,z;
         double v0;
         int sign0;

         MFNVector u;
         MFKVector s;

         if(MFPolytopeRadiusOfVertex(P,j,e)/MFChartRadius(C,e)>1.&& MFPolytopeGetVertexMark(P,j,e)==0)prt=1;

         u=MFCreateNVector(MFAtlasN(S,e),e);
         s=MFCreateKVector(MFAtlasK(S,e),e);

         MFPolytopeVertex(P,j,s,e);
         MFChartEvaluate(C,s,u,e);

         x=MFNV_C(u,0,e);
         y=MFNV_C(u,1,e);
         z=MFNV_C(u,2,e);

         MFFreeNVector(u,e);
         MFFreeKVector(s,e);

         v0=z-x-.5*(y-1)*(y-1);
         sign0=1;if(v0<0.)sign0=-1;if(fabs(v0)<1.e-7)sign0=0;
         if(sign0==1)prt=1;
        }

       if(prt)
        {
         double x,y,z;
         double v0;
         int sign0;

         MFNVector u;
         MFKVector s;

         u=MFCreateNVector(MFAtlasN(S,e),e);
         s=MFCreateKVector(MFAtlasK(S,e),e);

         printf("  chart %d, index1=%d, index2=%d,  polytope has %d vertices\n",i,MFNVGetIndex(MFChartCenter(C,e),e),MFNVGetIndex2(MFChartCenter(C,e),e),m);fflush(stdout);
         for(j=0;j<m;j++)
          {
           vMark=MFPolytopeGetVertexMark(P,j,e);

           MFPolytopeVertex(P,j,s,e);
           MFChartEvaluate(C,s,u,e);
  
           x=MFNV_C(u,0,e);
           y=MFNV_C(u,1,e);
           z=MFNV_C(u,2,e);

           v0=z-x-.5*(y-1)*(y-1);
           sign0=1;if(v0<0.)sign0=-1;if(fabs(v0)<1.e-7)sign0=0;

           printf("     vertex %d has mark %d, sign %2d, vertexRadius/chartRadius %lf (%lf,%lf,%lf)\n",j,vMark,sign0,MFPolytopeRadiusOfVertex(P,j,e)/MFChartRadius(C,e),x,y,z);fflush(stdout);
          }

         MFFreeNVector(u,e);
         MFFreeKVector(s,e);
        }
      }
    }
  }

  MFCloseAtlas(H,S,e);

  MFFreeAtlas(S,e);

  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0[0],e);
  MFFreeNVector(u0[1],e);
  MFFreeContinuationMethod(H,e);
  MFFreeErrorHandler(e);

  return;
 }
