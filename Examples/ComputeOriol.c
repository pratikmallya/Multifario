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

static char *id="@(#) $Id: ComputeOriol.c,v 1.2 2011/07/21 17:43:45 mhender Exp $";

#include <MFAtlas.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <math.h>
#include <MFFortran.h>
#include <MFMultifariosMethod.h>

int DrawOriol(MFNVector,double*,void*,MFErrorHandler);

static void F(int *n,double *z,int *m,double *f, void *data, MFErrorHandler e)
 {
	 
  double L1 = 0.4;
  double L2 = 0.2;
  double L3 = 0.1;
  double x,y,c1,s1,c2,s2,c3,s3,t,xi1,xi2,xi3,xi4,xi5;
  x = z[0];
  y = z[1];
  c1 = z[2];
  s1 = z[3];
  c2 = z[4];
  s2 = z[5];
  c3 = z[6];
  s3 = z[7];
  t = z[8];
  xi1 = z[9]; 
  xi2 = z[10];
  xi3 = z[11];
  xi4 = z[12];
  xi5 = z[13];
	 
  f[0] = L1*c1 + L2*c2 + L3*c3-x;
  f[1] = L1*s1 + L2*s2 + L3*s3-y;
  f[2] = c1*c1 + s1*s1 - 1;
  f[3] = c2*c2 + s2*s2 - 1;
  f[4] = c3*c3 + s3*s3 - 1;
  f[5] = c1 - t*t - 0.5;
  f[6] = L1*xi1 + 2*c1*xi3 + xi5;
  f[7] = L1*xi2 + 2*s1*xi3;
  f[8] = L2*xi1 + 2*c2*xi4;
  f[9] = L2*xi2 + 2*s2*xi4;
  f[10] = -2*t*xi5;
  f[11] = xi1*xi1 + xi2*xi2 + xi3*xi3 + xi4*xi4 + xi5*xi5 - 1;
	 
  return;
 }

static void dF(int *n,double *z,int *m,double *df, void *data, MFErrorHandler e)
 {
  double L1 = 0.4;
  double L2 = 0.2;
  double L3 = 0.1;
  int i,j;
  int N = 12;
  int V = 14;
  double x,y,c1,s1,c2,s2,c3,s3,t,xi1,xi2,xi3,xi4,xi5;
  x = z[0];
  y = z[1];
  c1 = z[2];
  s1 = z[3];
  c2 = z[4];
  s2 = z[5];
  c3 = z[6];
  s3 = z[7];
  t = z[8];
  xi1 = z[9]; 
  xi2 = z[10];
  xi3 = z[11];
  xi4 = z[12];
  xi5 = z[13];

  for (i=0;i<N;i++) {
	  for(j=0;j<V;j++)df[i+N*j]=0;
  }
	 
  df[0+N*0] = -1;
  df[0+N*2] = L1;	 
  df[0+N*4] = L2;
  df[0+N*6] = L3;
  
  df[1+N*1] = -1;
  df[1+N*3] = L1;	 
  df[1+N*5] = L2;
  df[1+N*7] = L3;
  
  df[2+N*2] = 2*c1;
  df[2+N*3] = 2*s1;
	 
  df[3+N*4] = 2*c2;
  df[3+N*5] = 2*s2;
	 
  df[4+N*6] = 2*c3;
  df[4+N*7] = 2*s3;
	 
  df[5+N*2] = 1;
  df[5+N*8] = -2*t;
	 
  df[6+N*2] = 2*xi3;	
  df[7+N*3] = 2*xi3;
  df[8+N*4] = 2*xi4;
  df[9+N*5] = 2*xi4;
  df[10+N*8] = -2*xi5;
	 
  df[6+N*9] = L1;
  df[8+N*9] = L2;
	 
  df[7+N*10] = L1;
  df[9+N*10] = L2;
	 
  df[6+N*11] = 2*c1;
  df[7+N*11] = 2*s1;
	 
  df[8+N*12] = 2*c2;
  df[9+N*12] = 2*s2;
	 
  df[6+N*13] = 1;
  df[10+N*13] = -2*t;
	 
  df[11+N*9] = 2*xi1;
  df[11+N*10] = 2*xi2;
  df[11+N*11] = 2*xi3;
  df[11+N*12] = 2*xi4;
  df[11+N*13] = 2*xi5;

  return;
 }

static void ddF(int *n,double *z,int *m, double *ddf, void *data, MFErrorHandler e)
 {
	 
	 double L1 = 0.4;
	 double L2 = 0.2;
	 double L3 = 0.1;
	 int i,j,k;
	 int N = 12;
	 int V = 14;
	 double x,y,c1,s1,c2,s2,c3,s3,t,xi1,xi2,xi3,xi4,xi5;
	 x = z[0];
	 y = z[1];
	 c1 = z[2];
	 s1 = z[3];
	 c2 = z[4];
	 s2 = z[5];
	 c3 = z[6];
	 s3 = z[7];
	 t = z[8];
	 xi1 = z[9]; 
	 xi2 = z[10];
	 xi3 = z[11];
	 xi4 = z[12];
	 xi5 = z[13];
	 
	 for (k=0;k<V;k++) {
		 for (i=0;i<N;i++) {
			 for(j=0;j<V;j++)ddf[i+N*(j+V*k)]=0;
		 }
	 }

	 
  ddf[2+N*(2+V*2)]= 2;
  ddf[2+N*(3+V*3)]= 2;
	 
  ddf[3+N*(4+V*4)]= 2;
  ddf[3+N*(5+V*5)]= 2;
	
  ddf[4+N*(6+V*6)]= 2;
  ddf[4+N*(7+V*7)]= 2;
	 
  ddf[5+N*(8+V*8)]= -2;
	 
  ddf[6+N*(2+V*11)]= 2;
  ddf[7+N*(3+V*11)]= 2;
  ddf[8+N*(4+V*12)]= 2;
  ddf[9+N*(5+V*12)]= 2;
  ddf[10+N*(8+V*13)]= -2;
	 
  ddf[6+N*(11+V*2)]= 2;
  ddf[7+N*(11+V*3)]= 2;
  
  ddf[8+N*(12+V*4)]= 2;
  ddf[9+N*(12+V*5)]= 2;
	 
  ddf[10+N*(13+V*8)]= -2;
  
  ddf[11+N*(9+V*9)]= 2;
  ddf[11+N*(10+V*10)]= 2;
  ddf[11+N*(11+V*11)]= 2;
  ddf[11+N*(12+V*12)]= 2;
  ddf[11+N*(13+V*13)]= 2;

  return;
 }

int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  int n,k;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0[6];
  MFContinuationMethod H;
  MFErrorHandler e;
  MFNVector ll,ur;

  e=MFCreateErrorHandler();

  n=14;
  k=2;

  M=MFIMFCreateAlgebraicSubroutine(n,k,F,dF,ddF,NULL,e);
  MFIMFSetProjectForDraw(M,DrawOriol,e);
  MFIMFSetR(M,.17,e);
	 
	 ll=MFIMFVectorFactory(M,e);
	 ur=MFIMFVectorFactory(M,e);
	 MFNVSetC(ll,0,-0.8,e);MFNVSetC(ur,0,0.8,e);
	 MFNVSetC(ll,1,-0.8,e);MFNVSetC(ur,1,0.8,e);
	 MFNVSetC(ll,2,-1.1,e);MFNVSetC(ur,2,1.1,e);
	 MFNVSetC(ll,3,-1.1,e);MFNVSetC(ur,3,1.1,e); 
	 MFNVSetC(ll,4,-1.1,e);MFNVSetC(ur,4,1.1,e);
	 MFNVSetC(ll,5,-1.1,e);MFNVSetC(ur,5,1.1,e);
	 MFNVSetC(ll,6,-1.1,e);MFNVSetC(ur,6,1.1,e);
	 MFNVSetC(ll,7,-1.1,e);MFNVSetC(ur,7,1.1,e);
	 MFNVSetC(ll,8,-1.1,e);MFNVSetC(ur,8,1.1,e);
	 MFNVSetC(ll,9,-1.1,e);MFNVSetC(ur,9,1.1,e);
	 MFNVSetC(ll,10,-1.1,e);MFNVSetC(ur,10,1.1,e);
	 MFNVSetC(ll,11,-1.1,e);MFNVSetC(ur,11,1.1,e);
	 MFNVSetC(ll,12,-1.1,e);MFNVSetC(ur,12,1.1,e);
	 MFNVSetC(ll,13,-1.1,e);MFNVSetC(ur,13,1.1,e);

         MFNVSetC(ll,6,-0.91,e);MFNVSetC(ur,6,1.1,e);  /* Cuts off long polygons */

	 Omega=MFNRegionCreateHyperCubeByCorners(n,ll,ur,e);
	 
	 u0[0]=MFIMFVectorFactory(M,e);
	 MFNVSetC(u0[0],0,0.403995608135 ,e);
	 MFNVSetC(u0[0],1,-0.126455721958 ,e);
	 MFNVSetC(u0[0],2,0.5 ,e);
	 MFNVSetC(u0[0],3,-0.866025403784 ,e);
	 MFNVSetC(u0[0],4,0.675821357395 ,e);
	 MFNVSetC(u0[0],5,0.737065460383 ,e);
	 MFNVSetC(u0[0],6,0.688313366558 ,e);
	 MFNVSetC(u0[0],7,0.725413474797 ,e);
	 MFNVSetC(u0[0],8,2.73130897604e-24 ,e);
	 MFNVSetC(u0[0],9,0.608612150262 ,e);
	 MFNVSetC(u0[0],10,0.663765638387 ,e);
	 MFNVSetC(u0[0],11,0.153290108001 ,e);
	 MFNVSetC(u0[0],12,-0.0900551815361 ,e);
	 MFNVSetC(u0[0],13,-0.396734968105 ,e);
	 
	 /*
	 u0[1]=MFIMFVectorFactory(M,e);
	 MFNVSetC(u0[1],0,0.403923565905,e);
	 MFNVSetC(u0[1],1,0.184753746631,e);
	 MFNVSetC(u0[1],2,0.5,e);
	 MFNVSetC(u0[1],3,0.866025403784,e);
	 MFNVSetC(u0[1],4,0.521636629833,e);
	 MFNVSetC(u0[1],5,-0.853167759832,e);
	 MFNVSetC(u0[1],6,0.99596239938,e);
	 MFNVSetC(u0[1],7,0.0897713708364,e);
	 MFNVSetC(u0[1],8,-1.29479472549e-23,e);
	 MFNVSetC(u0[1],9,-0.473549944913,e);
	 MFNVSetC(u0[1],10,0.774519124164,e);
	 MFNVSetC(u0[1],11,-0.178867529931,e);
	 MFNVSetC(u0[1],12,0.0907815743432,e);
	 MFNVSetC(u0[1],13,0.368287507897,e);
	 
	 u0[2]=MFIMFVectorFactory(M,e);
	 MFNVSetC(u0[2],0,0.0757147689592,e);
	 MFNVSetC(u0[2],1,0.470719378113,e);
	 MFNVSetC(u0[2],2,0.5,e);
	 MFNVSetC(u0[2],3,0.866025403784,e);
	 MFNVSetC(u0[2],4,-0.965851888917,e);
	 MFNVSetC(u0[2],5,0.259094825645,e);
	 MFNVSetC(u0[2],6,0.688851467426,e);
	 MFNVSetC(u0[2],7,0.724902514705,e);
	 MFNVSetC(u0[2],8,-1.43285554384e-24,e);
	 MFNVSetC(u0[2],9,0.877086012116,e);
	 MFNVSetC(u0[2],10,-0.235282914484,e);
	 MFNVSetC(u0[2],11,0.0543362616053,e);
	 MFNVSetC(u0[2],12,0.090809576725,e);
	 MFNVSetC(u0[2],13,-0.405170666452,e);
	 
	 u0[3]=MFIMFVectorFactory(M,e);
	 MFNVSetC(u0[3],0,0.643962044525,e);
	 MFNVSetC(u0[3],1,-0.11714133967,e);
	 MFNVSetC(u0[3],2,0.951167303112,e);
	 MFNVSetC(u0[3],3,-0.308675819415,e);
	 MFNVSetC(u0[3],4,0.951167303112,e);
	 MFNVSetC(u0[3],5,-0.308675819415,e);
	 MFNVSetC(u0[3],6,0.732616626579,e);
	 MFNVSetC(u0[3],7,0.680641519789,e);
	 MFNVSetC(u0[3],8,0.671689886117,e);
	 MFNVSetC(u0[3],9,-0.928244240493,e);
	 MFNVSetC(u0[3],10,0.301236754684,e);
	 MFNVSetC(u0[3],11,0.19518001459,e);
	 MFNVSetC(u0[3],12,0.0975900072949,e);
	 MFNVSetC(u0[3],13,8.05386342117e-23,e);
	 
	 u0[4]=MFIMFVectorFactory(M,e);
	 MFNVSetC(u0[4],0,0.0997769115265 ,e);
	 MFNVSetC(u0[4],1,-0.0302080352413 ,e);
	 MFNVSetC(u0[4],2,0.880415758677 ,e);
	 MFNVSetC(u0[4],3,-0.474202585267 ,e);
	 MFNVSetC(u0[4],4,-0.880415758677 ,e);
	 MFNVSetC(u0[4],5,0.474202585267 ,e);
	 MFNVSetC(u0[4],6,-0.763062402089 ,e);
	 MFNVSetC(u0[4],7,0.646324818121 ,e);
	 MFNVSetC(u0[4],8,0.616778532925 ,e);
	 MFNVSetC(u0[4],9,-0.859197803118 ,e);
	 MFNVSetC(u0[4],10,0.462774337554 ,e);
	 MFNVSetC(u0[4],11,0.19518001459 ,e);
	 MFNVSetC(u0[4],12,-0.0975900072949 ,e);
	 MFNVSetC(u0[4],13,-7.56094778662e-24 ,e);
	 
	 u0[5]=MFIMFVectorFactory(M,e);
	 MFNVSetC(u0[5],0,-0.0515163256928,e);
	 MFNVSetC(u0[5],1,-0.2688494880480,e);
	 MFNVSetC(u0[5],2,0.5000000000000,e);
	 MFNVSetC(u0[5],3,-0.8660254037840,e);
	 MFNVSetC(u0[5],4,-0.9991894614680,e);
	 MFNVSetC(u0[5],5,-0.0402544418940,e);
	 MFNVSetC(u0[5],6,-0.5167843339930,e);
	 MFNVSetC(u0[5],7,0.8561156184410,e);
	 MFNVSetC(u0[5],8,0.0000000000000,e);
	 MFNVSetC(u0[5],9,0.9208651533680,e);
	 MFNVSetC(u0[5],10,0.0370989829637,e);
	 MFNVSetC(u0[5],11,0.0085676431203,e);
	 MFNVSetC(u0[5],12,0.0921612155532,e);
	 MFNVSetC(u0[5],13,-0.3769137044670,e);
	 
	 u0[6]=MFIMFVectorFactory(M,e);
	 MFNVSetC(u0[6],0,0.05 ,e);
	 MFNVSetC(u0[6],1,0.086602540378 ,e);
	 MFNVSetC(u0[6],2,0.5 ,e);
	 MFNVSetC(u0[6],3,0.8660254037840 ,e);
	 MFNVSetC(u0[6],4,-0.5 ,e);
	 MFNVSetC(u0[6],5,-0.8660254037840 ,e);
	 MFNVSetC(u0[6],6,-0.5 ,e);
	 MFNVSetC(u0[6],7,-0.8660254037840 ,e);
	 MFNVSetC(u0[6],8,0.0 ,e);
	 MFNVSetC(u0[6],9,0.487950036474 ,e);
	 MFNVSetC(u0[6],10,0.845154254729 ,e);
	 MFNVSetC(u0[6],11,-0.195180014590 ,e);
	 MFNVSetC(u0[6],12,0.097590007295 ,e);
	 MFNVSetC(u0[6],13,0.0 ,e);
	  
	 u0[7]=MFIMFVectorFactory(M,e);
	 MFNVSetC(u0[7],0,0.35 ,e);
	 MFNVSetC(u0[7],1,0.606217782649 ,e);
	 MFNVSetC(u0[7],2,0.5 ,e);
	 MFNVSetC(u0[7],3,0.8660254037840 ,e);
	 MFNVSetC(u0[7],4,0.5 ,e);
	 MFNVSetC(u0[7],5,0.8660254037840 ,e);
	 MFNVSetC(u0[7],6,0.5 ,e);
	 MFNVSetC(u0[7],7,0.8660254037840 ,e);
	 MFNVSetC(u0[7],8,0.0 ,e);
	 MFNVSetC(u0[7],9,0.487950036474 ,e);
	 MFNVSetC(u0[7],10,0.845154254729 ,e);
	 MFNVSetC(u0[7],11,-0.195180014590 ,e);
	 MFNVSetC(u0[7],12,-0.097590007295 ,e);
	 MFNVSetC(u0[7],13,0.0 ,e);
	 */
             
  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.05,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",0,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  MFMultifarioSetIntegerParameter(H,"branchSwitch",0,e);
  MFMultifarioSetRealParameter(H,"DotMin",0.9,e);
  MFMultifarioSetFilename(H,"Oriol",e);

  S=MFComputeAtlasMultiple(H,M,Omega,1,u0,e);

  MFCloseAtlas(H,S,e);
  printf("Done computating Atlas\n");fflush(stdout);

  MFFreeAtlas(S,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0[0],e);
/*
	 MFFreeNVector(u0[1],e);
	 MFFreeNVector(u0[2],e);
	 MFFreeNVector(u0[3],e);
	 MFFreeNVector(u0[4],e);
	 MFFreeNVector(u0[5],e);
 */
	 /*MFFreeNVector(u0[6],e);*/
  MFFreeNVector(ll,e);
  MFFreeNVector(ur,e);
  MFFreeContinuationMethod(H,e);

  MFFreeErrorHandler(e);

  return 0;
 }

int DrawOriol(MFNVector v,double *u,void *d, MFErrorHandler e)
 {
  if(u==(double*)NULL||v==(MFNVector)NULL) return 3;

  u[0]=MFNV_C(v,0,e);
  u[1]=MFNV_C(v,1,e);
  u[2]=atan2(MFNV_C(v,7,e),MFNV_C(v,6,e));
/*u[2]=MFNV_C(v,2,e);*/

  return 0;
 }
