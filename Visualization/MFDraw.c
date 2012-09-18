/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *  author: Mike Henderson mhender@watson.ibm.com
 */

static char *id="@(#) $Id: MFDraw.c,v 1.5 2011/07/21 17:44:29 mhender Exp $";

static char MFDrawErrorMsg[256]="";

#include <MFErrorHandler.h>
#include <MFChart.h>
#include <MFImplicitMF.h>
#include <MFErrorHandler.h>
#include <MFNVector.h>
#include <MFKVector.h>
#include <MFEnumPolytope.h>
#include <math.h>
#include <stdlib.h>
#include <sh.h>
#include <MFPrint.h>
#include <MFDraw.h>
#include <MFDrawClippedSphere.h>
#include <string.h>
#include <MFDrawChartState.h>
#include <MFEnumPolytope.h>
#include <MFEnumDualPolytope.h>

#ifdef __cplusplus
 extern "C" {
#endif

MFAtlas MFMFA;
static int MFPendulaPeriodicFlag=0;

#define PERIODIC
#define sqrt3 1.73205080
#define MFPI 3.14159265358979323846264338327950288

void MFDrawLineOnChart(MFChart,MFKVector,MFKVector,MFChartState,MFErrorHandler);
void MFDrawLineOnChartTS(MFChart,MFKVector,MFKVector,MFErrorHandler);
void MFDrawLine(MFNVector,MFNVector,MFChartState,MFErrorHandler);
void MFMarkPoint(MFChart,MFNVector,int,MFChartState,MFErrorHandler);
void MFMarkPointOnChart(MFChart,MFKVector,int,MFChartState,MFErrorHandler);
void MFDrawProjectTo3d(MFChart,MFNVector,float*,float*,float*,MFErrorHandler);
void MFDrawPolytope(MFPolytope,MFChart,MFChartState,MFErrorHandler);
/* void MFDrawEnumPolytope(MFEnumPolytope,MFChart,MFChartState,MFErrorHandler); */
/* void MFDrawEnumDualPolytope(MFEnumDualPolytope,MFErrorHandler); */
/* void MFDrawEnumDualPolytopeEdges(MFEnumDualPolytope,MFErrorHandler); */
MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);
void MFDrawMakeListOfTriangles(MFChart,int*,int*,double**,double,MFErrorHandler);
double *MFNV_CStar(MFNVector,MFErrorHandler);
double *MFKV_CStar(MFKVector,MFErrorHandler);
void MFDrawChartTS(MFChart,MFErrorHandler);
void MFDraw1dChartTS(MFChart,MFErrorHandler);
void MFDraw2dChartTS(MFChart,MFErrorHandler);
void MFDrawChartBoundaryTS(MFChart,MFErrorHandler);
void MFDraw1dChartBoundaryTS(MFChart,MFErrorHandler);
void MFDraw2dChartBoundaryTS(MFChart,MFErrorHandler);
void MFDrawPolytopeTS(MFPolytope,MFChart,MFErrorHandler);
void MFDrawEnumPolytopeTS(MFEnumPolytope,MFChart,MFErrorHandler);

MFKVector MFTmpK=NULL;
MFNVector MFTmpN=NULL;

int MFDrawSetExt=1;
double MFDrawXMax=-1.e20;
double MFDrawXMin= 1.e20;
double MFDrawYMax=-1.e20;
double MFDrawYMin= 1.e20;
double MFDrawZMax=-1.e20;
double MFDrawZMin= 1.e20;

double MFDrawCubeXMax=1.;
double MFDrawCubeXMin=-1.;
double MFDrawCubeYMax=1.;
double MFDrawCubeYMin=-1.;
double MFDrawCubeZMax=1.;
double MFDrawCubeZMin=-1.;
int MFDrawCube=0;

void MFTestClippedSphere(MFChart,double,double,double,MFErrorHandler);

#ifdef NEW  /* Not written, but will be needed someday */
void MFDrawSlicedChartTS(MFChart,int,double*,double,MFErrorHandler);
void MFDrawSliced1dChartTS(MFChart,int,double*,double,MFErrorHandler);
void MFDrawSliced2dChartTS(MFChart,int,double*,double,MFErrorHandler);
void MFCreateSectionOfChart(MFChart chart,int m,double *n,double *on,int *knew, double **s0, souble **Psi,MFErrorHandler);
#endif

void MFDrawInitialize(float alpha, float beta, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawInitialize"};
  int grey=100;
  int full=255;
  float dist,xmin,xmax;

/*shinit(&grey,&grey,&grey);*/
  shinit(&full,&full,&full);
  dist=100.;
  xmin=-1.;
  xmax= 1.;
  shview(&dist,&alpha,&beta,&xmin,&xmax,&xmin,&xmax,&xmin,&xmax);
  MFDrawCube=0;
  return;
 }

void MFDrawInitializeFromFile(char *name, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawInitializeFromFile"};
  float alpha=0.;
  float beta=0.;
  float xmin=0.;
  float xmax=0.;
  float ymin=0.;
  float ymax=0.;
  float zmin=0.;
  float zmax=0.;
  int cube=0;
  FILE *fid=NULL;
  char *fullname=NULL;
  int zero=0;
  int one=1;
  float x=0.;
  float y=0.;
  float z=0.;
  float rone=1.;
  float rzero=0.;
  printf("MFDrawInitializeFromFile\n");fflush(stdout);

  fullname=(char*)malloc((strlen(name)+1+5)*sizeof(char));

#ifndef MFNOSAFETYNET
  if(fullname==NULL)
   {
    sprintf(MFDrawErrorMsg,"Out of memory trying to allocate %d bytes.",(strlen(name)+1+5)*sizeof(char));
    MFSetError(e,12,RoutineName,MFDrawErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  strcpy(fullname,name);
  strcat(fullname,".view");
  fid=fopen(fullname,"r");
  free(fullname);
  if(fid==NULL)
   {
    alpha=0.;
    beta=0.;
    xmin=-1.;
    xmax= 1.;
    ymin=-1.;
    ymax= 1.;
    zmin=-1.;
    zmax= 1.;
    cube= 0;
    printf("Couldn't find view file \"%s\", using defaults.\n",fullname);
   }else{
    fscanf(fid,"%f %f %f %f %f %f %f %f %d",&alpha,&beta,&xmin,&xmax,&ymin,&ymax,&zmin,&zmax,&cube);
    fclose(fid);
   }
  printf("%f %f %f %f %f %f %f %f %d\n",alpha,beta,xmin,xmax,ymin,ymax,zmin,zmax,cube);fflush(stdout);

  x=0.;
  y=0.;
  z=0.;

  if(!cube)MFDrawInitializeNoCube(alpha,beta,xmin,xmax,ymin,ymax,zmin,zmax,e);
   else MFDrawInitializeCube(alpha,beta,xmin,xmax,ymin,ymax,zmin,zmax,e);

  MFDrawCube=cube;

  return;
 }

void MFDrawInitializeCube(float alpha, float beta, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawInitializeCube"};
  int grey=100;
  int full=255;
  int zero=0;
  float dist;

/*shinit(&grey,&grey,&grey);*/
  shinit(&full,&full,&full);
  dist=34.*(xmax-xmin+ymax-ymin+zmax-zmin);
  shview(&dist,&alpha,&beta,&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);
  shlinc(&zero,&zero,&zero);
  shcube(&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);
  MFDrawCubeXMax=xmax;
  MFDrawCubeXMin=xmin;
  MFDrawCubeYMax=ymax;
  MFDrawCubeYMin=ymin;
  MFDrawCubeZMax=zmax;
  MFDrawCubeZMin=zmin;
  MFDrawCube=1;
  return;
 }

void MFDrawInitializeNoCube(float alpha, float beta, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawInitializeCube"};
  int grey=100;
  int full=255;
  int zero=0;
  float dist;

/*shinit(&grey,&grey,&grey);*/
  shinit(&full,&full,&full);
  dist=100.;
  shview(&dist,&alpha,&beta,&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);
  shlinc(&zero,&zero,&zero);
  MFDrawCube=0;
  return;
 }

void MFDrawClose(MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawClose"};
  shend();
  MFFreeDrawCharts(e);
  if(MFTmpK!=NULL)MFFreeKVector(MFTmpK,e);
  if(MFTmpN!=NULL)MFFreeNVector(MFTmpN,e);
  return;
 }

void MFDrawClear(MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawClear"};

  shclr();

  if(MFDrawCube)
   {
    float xmin,xmax,ymin,ymax,zmin,zmax;
    int zero=0;
    shlinc(&zero,&zero,&zero);
    xmin=MFDrawCubeXMin;
    xmax=MFDrawCubeXMax;
    ymin=MFDrawCubeYMin;
    ymax=MFDrawCubeYMax;
    zmin=MFDrawCubeZMin;
    zmax=MFDrawCubeZMax;
    shcube(&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);
   }

  return;
 }

void MFDrawDisplay(MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawDisplay"};
  shpause();
  return;
 }

void MFDrawInterval(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawInterval"};
  return;
 }

void MFDrawPolygon(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawPolygon"};
  return;
 }

void MFDrawPolyhedron(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawPolyhedron"};
  return;
 }

void MFDrawPolytope(MFPolytope P,MFChart chart,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawPolytope"};
  MFEnumPolytope EP;

  if(MFPolytopeDimension(P,e)>3)return;

  EP=MFEnumeratePolytope(P,e);
  MFDrawEnumPolytope(EP,chart,state,e);
  MFFreeEnumPolytope(EP,e);

  return;
 }

void MFDrawPolytopeTS(MFPolytope P,MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawPolytopeTS"};
  MFEnumPolytope EP;

  if(MFPolytopeDimension(P,e)>3)return;

  EP=MFEnumeratePolytope(P,e);
  MFDrawEnumPolytopeTS(EP,chart,e);
  MFFreeEnumPolytope(EP,e);

  return;
 }

void MFDraw1dChart(MFChart chart,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDraw1dChart"};
  int color;
  int zero=0;
  int half=128;
  int full=255;
  MFPolytope P;
  MFKVector s0,s1;
  int n,k;

  P=MFChartPolytope(chart,e);
  k=MFChartK(chart,e);
  n=MFChartN(chart,e);
  if(MFPolytopeNumberOfVertices(P,e)<2)return;
  s0=MFCreateKVector(k,e);
  s1=MFCreateKVector(k,e);

  MFPolytopeVertex(P,0,s0,e);
  MFPolytopeVertex(P,1,s1,e);

  MFDrawLineOnChart(chart,s0,s1,state,e);
  MFMarkPointOnChart(chart,s0,2,state,e);
  MFMarkPointOnChart(chart,s1,2,state,e);
  MFMarkPoint(chart,MFChartCenter(chart,e),3,state,e);

  MFFreeKVector(s0,e);
  MFFreeKVector(s1,e);

  return;
 }

void MFDraw2dChart(MFChart chart,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDraw2dChart"};
  int color;
  int zero=0;
  int half=128;
  int full=255;
  double *tri;
  int ntri,mtri;
  static float x[3],y[3],z[3];
  double t;
  MFKVector s;
  MFNVector u;
  int j,n;
  int rc;

/* The interior of the circle */

  tri=(double*)malloc(6*6*sizeof(double));

#ifndef MFNOSAFETYNET
  if(tri==NULL)
   {
    sprintf(MFDrawErrorMsg,"Out of memory trying to allocate %d bytes.",6*6*sizeof(double));
    MFSetError(e,12,RoutineName,MFDrawErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  t=2./sqrt3*MFChartRadius(chart,e);

  tri[0+6*0]=0.;
  tri[1+6*0]=0.;
  tri[2+6*0]=t;
  tri[3+6*0]=0.;
  tri[4+6*0]=.5*t;
  tri[5+6*0]=.5*sqrt3*t;

  tri[0+6*1]=0.;
  tri[1+6*1]=0.;
  tri[2+6*1]=.5*t;
  tri[3+6*1]=.5*sqrt3*t;
  tri[4+6*1]=-.5*t;
  tri[5+6*1]=.5*sqrt3*t;

  tri[0+6*2]=0.;
  tri[1+6*2]=0.;
  tri[2+6*2]=-.5*t;
  tri[3+6*2]=.5*sqrt3*t;
  tri[4+6*2]=-t;
  tri[5+6*2]=0.;

  tri[0+6*3]=0.;
  tri[1+6*3]=0.;
  tri[2+6*3]=-t;
  tri[3+6*3]=0.;
  tri[4+6*3]=-.5*t;
  tri[5+6*3]=-.5*sqrt3*t;

  tri[0+6*4]=0.;
  tri[1+6*4]=0.;
  tri[2+6*4]=-.5*t;
  tri[3+6*4]=-.5*sqrt3*t;
  tri[4+6*4]= .5*t;
  tri[5+6*4]=-.5*sqrt3*t;

  tri[0+6*5]=0.;
  tri[1+6*5]=0.;
  tri[2+6*5]= .5*t;
  tri[3+6*5]=-.5*sqrt3*t;
  tri[4+6*5]=t;
  tri[5+6*5]=0.;

  ntri=6;
  mtri=6;
  t=.2*MFChartRadius(chart,e);
  if(t<.01)t=.01;
  MFDrawMakeListOfTriangles(chart,&ntri,&mtri,&tri,t,e);

  s=MFCreateKVector(MFChartK(chart,e),e);
  u=MFCreateNVector(MFChartN(chart,e),e);

  for(j=0;j<ntri;j++)
   {
    MFKVSetC(s,0,tri[0+6*j],e);
    MFKVSetC(s,1,tri[1+6*j],e);
    rc=MFChartEvaluate(chart,s,u,e);
    MFDrawProjectTo3d(chart,u,x,y,z,e);
    MFKVSetC(s,0,tri[2+6*j],e);
    MFKVSetC(s,1,tri[3+6*j],e);
    rc=rc && MFChartEvaluate(chart,s,u,e);
    MFDrawProjectTo3d(chart,u,x+1,y+1,z+1,e);
    MFKVSetC(s,0,tri[4+6*j],e);
    MFKVSetC(s,1,tri[5+6*j],e);
    rc = rc && MFChartEvaluate(chart,s,u,e);
    MFDrawProjectTo3d(chart,u,x+2,y+2,z+2,e);

#ifdef PERIODIC
    if(MFPendulaPeriodicFlag)
     {
      rc=rc&&fabs(x[0]-x[1])<.1;
      rc=rc&&fabs(x[1]-x[2])<.1;
      rc=rc&&fabs(x[2]-x[0])<.1;
      rc=rc&&fabs(y[0]-y[1])<.1;
      rc=rc&&fabs(y[1]-y[2])<.1;
      rc=rc&&fabs(y[2]-y[0])<.1;
     }
#endif
    n=3;
    if(rc&&state!=NULL)MFChartStateAddPolygon(state,n,x,y,z,0,e);
    if(rc&&state==NULL)shpg(&n,x,y,z);
   }
  free(tri);
  MFFreeKVector(s,e);
  MFFreeNVector(u,e);

  MFDrawPolytope(MFChartPolytope(chart,e),chart,state,e);

  return;
 }

void MFDraw3dChart(MFChart chart,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDraw3dChart"};
  int full=255;

  if(MFChartK(chart,e)>2)shtric(&full,&full,&full);

  MFDrawPolytope(MFChartPolytope(chart,e),chart,state,e);

  return;
 }

void MFDraw1dChartBoundary(MFChart chart,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDraw1dChartBoundary"};
  MFKVector s;
  int in0;
  int zero=0;
  int full=255;

  s=MFCreateKVector(MFChartK(chart,e),e);
  MFKVSetC(s,0,MFChartRadius(chart,e),e);

  in0=MFChartInterior(chart,s,e);
  if(in0)MFMarkPointOnChart(chart,s,0,state,e);

  MFKVSetC(s,0,-MFChartRadius(chart,e),e);

  in0=MFChartInterior(chart,s,e);
  if(in0)MFMarkPointOnChart(chart,s,1,state,e);

  MFFreeKVector(s,e);
  return;
 }

void MFDraw1dChartBoundaryTS(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDraw1dChartBoundary"};
  MFKVector s;
  MFNVector u;
  int in0;
  int zero=0;
  int full=255;

  s=MFCreateKVector(MFChartK(chart,e),e);
  MFKVSetC(s,0,MFChartRadius(chart,e),e);

  in0=MFChartInterior(chart,s,e);
  if(in0)
   {
    u=MFCreateNVector(MFChartN(chart,e),e);
    MFChartPointInTangentSpace(chart,s,u,e);
    MFMarkPoint(chart,u,0,NULL,e);
    MFFreeNVector(u,e);
   }

  MFKVSetC(s,0,-MFChartRadius(chart,e),e);

  in0=MFChartInterior(chart,s,e);
  if(in0)
   {
    u=MFCreateNVector(MFChartN(chart,e),e);
    MFChartPointInTangentSpace(chart,s,u,e);
    MFMarkPoint(chart,u,1,NULL,e);
    MFFreeNVector(u,e);
   }

  MFFreeKVector(s,e);
  return;
 }

void MFDraw2dChartBoundary(MFChart chart,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDraw2dChartBoundary"};
  MFKVector s;
  MFNVector u;
  int zero=0;
  int full=255;
  int i;
  double r;
  int in0,in1;
  double pt[4]={0.,0.,0.,0.};
  float xa=0;
  float ya=0;
  float za=0.;
  float xb=0;
  float yb=0;
  float zb=0.;
  float x[2]={0.,0.};
  float y[2]={0.,0.};
  float z[2]={0.,0.};
  int rc;
  int n;
  double dt;
  int rr,gg,bb,rs,gs,bs;

/* Draw yellow circle for boundary of ball */

  s=MFCreateKVector(MFChartK(chart,e),e);
  u=MFCreateNVector(MFChartN(chart,e),e);

  r=MFChartRadius(chart,e);

  n=400;
  if(r<.05)n=6;
  dt=2./n*3.1415926;
  shqlinc(&rs,&gs,&bs);
  rr=255;
  gg=24;
  bb=47;
  shlinc(&rr,&gg,&bb);

  for(i=0;i<n;i++)
   {
    pt[0]=.999999*r*cos(dt*i);
    pt[1]=.999999*r*sin(dt*i);
    MFKVSetC(s,0,pt[0],e);
    MFKVSetC(s,1,pt[1],e);
    in0=MFChartInterior(chart,s,e);
    pt[2]=.999999*r*cos(dt*(i+1));
    pt[3]=.999999*r*sin(dt*(i+1));
    MFKVSetC(s,0,pt[2],e);
    MFKVSetC(s,1,pt[3],e);
    in1=MFChartInterior(chart,s,e);

    if(in0&&in1)
     {
      MFKVSetC(s,0,pt[0],e);
      MFKVSetC(s,1,pt[1],e);
      rc=MFChartEvaluate(chart,s,u,e);
      MFDrawProjectTo3d(chart,u,&xa,&ya,&za,e);

      MFKVSetC(s,0,pt[2],e);
      MFKVSetC(s,1,pt[3],e);
      rc=rc && MFChartEvaluate(chart,s,u,e);

      MFDrawProjectTo3d(chart,u,&xb,&yb,&zb,e);

      x[0]=xa;
      y[0]=ya;
      z[0]=za;
      x[1]=xb;
      y[1]=yb;
      z[1]=zb;
#ifdef PERIODIC
      if(MFPendulaPeriodicFlag)
       {
        rc=rc&&fabs(xa-xb)<.1;
        rc=rc&&fabs(ya-yb)<.1;
       }
#endif
      if(rc&&state!=NULL)MFChartStateAddLine(state,x[0],y[0],z[0],x[1],y[1],z[1],e);
      if(rc&&state==NULL)shline(x,y,z,x+1,y+1,z+1);
     }
   }
  MFFreeKVector(s,e);
  MFFreeNVector(u,e);
  shlinc(&rs,&gs,&bs);

  return;
 }

void MFDraw2dChartBoundaryTS(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDraw2dChartBoundaryTS"};
  MFKVector s;
  MFNVector u;
  int zero=0;
  int full=255;
  int i;
  double r;
  int in0,in1;
  double pt[4]={0.,0.,0.,0.};
  float xa=0;
  float ya=0;
  float za=0.;
  float xb=0;
  float yb=0;
  float zb=0.;
  float x[2]={0.,0.};
  float y[2]={0.,0.};
  float z[2]={0.,0.};
  int rc;
  int n;
  double dt;
  int rr,gg,bb,rs,gs,bs;

/* Draw yellow circle for boundary of ball */

  s=MFCreateKVector(MFChartK(chart,e),e);
  u=MFCreateNVector(MFChartN(chart,e),e);

  r=MFChartRadius(chart,e);

  n=400;
  if(r<.05)n=6;
  dt=2./n*3.1415926;

  shqlinc(&rs,&gs,&bs);
  rr=255;
  gg=24;
  bb=47;
  shlinc(&rr,&gg,&bb);

  for(i=0;i<n;i++)
   {
    pt[0]=.999999*r*cos(dt*i);
    pt[1]=.999999*r*sin(dt*i);
    MFKVSetC(s,0,pt[0],e);
    MFKVSetC(s,1,pt[1],e);
    in0=MFChartInterior(chart,s,e);
    pt[2]=.999999*r*cos(dt*(i+1));
    pt[3]=.999999*r*sin(dt*(i+1));
    MFKVSetC(s,0,pt[2],e);
    MFKVSetC(s,1,pt[3],e);
    in1=MFChartInterior(chart,s,e);

    if(in0&&in1)
     {
      MFKVSetC(s,0,pt[0],e);
      MFKVSetC(s,1,pt[1],e);
      MFChartPointInTangentSpace(chart,s,u,e);
      MFDrawProjectTo3d(chart,u,&xa,&ya,&za,e);

      MFKVSetC(s,0,pt[2],e);
      MFKVSetC(s,1,pt[3],e);
      MFChartPointInTangentSpace(chart,s,u,e);

      MFDrawProjectTo3d(chart,u,&xb,&yb,&zb,e);

      x[0]=xa;
      y[0]=ya;
      z[0]=za;
      x[1]=xb;
      y[1]=yb;
      z[1]=zb;
      rc=1;
#ifdef PERIODIC
      if(MFPendulaPeriodicFlag)
       {
        rc=rc&&fabs(xa-xb)<.1;
        rc=rc&&fabs(ya-yb)<.1;
       }
#endif
      if(rc)shline(x,y,z,x+1,y+1,z+1);
     }
   }
  MFFreeKVector(s,e);
  MFFreeNVector(u,e);
  shlinc(&rs,&gs,&bs);

  return;
 }

void MFDraw3dChartBoundary(MFChart chart,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDraw3dChartBoundary"};
  int zero=0;
  int full=255;
  double e1,e2,e3;

  e1=sqrt(2.)/8.*MFChartRadius(chart,e)/2.;
  e2=e1/2.;
  e3=e2/50.;
  MFDrawClippedSphere(chart,e1,e2,e3,state,e);

  return;
 }

void MFTest3dChartBoundary(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTest3dChartBoundary"};
  int zero=0;
  int full=255;
  double e1,e2,e3;

  e1=sqrt(2.)/8.*MFChartRadius(chart,e)/2.;
  e2=e1/2.;
  e3=e2/50.;
  MFTestClippedSphere(chart,e1,e2,e3,e);

  return;
 }

void MFDrawChart(MFChart chart,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawChart"};

  if(chart==NULL)return;
  if(MFChartK(chart,e)==1)MFDraw1dChart(chart,state,e);
   else if(MFChartK(chart,e)==2)MFDraw2dChart(chart,state,e);
   else if(MFChartK(chart,e)==3)MFDraw3dChart(chart,state,e);

  return;
 }

void MFDrawChartBoundaryTS(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawChartBoundaryTS"};
  if(chart==NULL)return;
  if(MFChartK(chart,e)==1)MFDraw1dChartBoundaryTS(chart,e);
   else if(MFChartK(chart,e)==2)MFDraw2dChartBoundaryTS(chart,e);
   else if(MFChartK(chart,e)==3)MFDraw3dChartBoundary(chart,NULL,e);

  return;
 }

void MFDrawChartBoundary(MFChart chart,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawChartBoundary"};
  if(chart==NULL)return;
  if(MFChartK(chart,e)==1)MFDraw1dChartBoundary(chart,state,e);
   else if(MFChartK(chart,e)==2)MFDraw2dChartBoundary(chart,state,e);
   else if(MFChartK(chart,e)==3)MFDraw3dChartBoundary(chart,state,e);

  return;
 }

void MFDrawAtlas(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawAtlas"};
  int chart,bchart;
  MFChartState state;
  int rb=243;
  int gb=223;
  int bb=192;
  int verbose=0;

  MFDrawSetExt=1;
  MFDrawXMin=0.;
  MFDrawXMax=0.;
  MFDrawYMin=0.;
  MFDrawYMax=0.;
  MFDrawZMin=0.;
  MFDrawZMax=0.;


  for(chart=0;chart<MFAtlasNumberOfCharts(A,e);chart++)
   {
    printf("Draw Chart %d\n",chart);fflush(stdout);
    if(chart<MFNChartStates(e)) /* !!! */
      state=MFGetChartState(chart,e);
     else{
      state=MFAddChartState(MFAtlasChart(A,chart,e),e);
      MFDrawChart(MFAtlasChart(A,chart,e),state,e);
      MFMFA=A;
      MFDrawChartBoundary(MFAtlasChart(A,chart,e),state,e);
     }
    MFSetChartStateColor(state,0,180,180,205,e);
    MFSetChartStateColor(state,1,180,180,205,e);
   }

  for(chart=0;chart<MFAtlasNumberOfCharts(A,e);chart++)
   {
    if(MFAtlasIsChartNearBoundary(A,chart,e))
     {
      state=MFGetChartState(chart,e);
      MFSetChartStateColor(state,0,100,100,150,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("This one is near the boundary, %d\n",chart);fflush(stdout);}
#endif

     }
   }

  for(bchart=0;bchart<MFAtlasNumberOfChartsWithBoundary(A,e);bchart++)
   {
    chart=MFAtlasChartWithBoundary(A,bchart,e);
    if(bchart<MFNChartStates(e))
     {
      state=MFGetChartState(chart,e);
      MFSetChartStateColor(state,0,243,223,192,e);
      MFSetChartStateColor(state,1,243,223,192,e);
     }
   }

  for(chart=0;chart<MFAtlasNumberOfCharts(A,e);chart++)
   {
    state=MFGetChartState(chart,e);
    if(MFChartStateChanged(state,e))
     {
      MFClearChartState(state,e);
      MFDrawChart(MFAtlasChart(A,chart,e),state,e);
      MFDrawChartBoundary(MFAtlasChart(A,chart,e),state,e);
     }
    MFDrawChartFromState(state,e);
   }

  printf("at end of MFDrawAtlas, x in (%lf %lf %lf %lf %lf %lf)\n",MFDrawXMin,MFDrawXMax,MFDrawYMin,MFDrawYMax,MFDrawZMin,MFDrawZMax);fflush(stdout);

  return;
 }

void MFDrawAtlasOnce(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawAtlasOnce"};
  int chart,bchart;
  int verbose=0;
  int ri=62;
  int gi=158;
  int bi=255;
  int zero=0;
  int full=255;
  int rb=255;
  int gb=211;
  int bb=62;
  int r=0;
  int g=0;
  int b=0;
  int rnb=100;
  int gnb=100;
  int bnb=150;
  int rsb=150;
  int gsb= 80;
  int bsb= 80;

  shlinc(&zero,&zero,&zero);
  for(chart=0;chart<MFAtlasNumberOfCharts(A,e);chart++)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("%s chart %d/%d\n",RoutineName,chart,MFAtlasNumberOfCharts(A,e));fflush(stdout);}
#endif

    shtric(&ri,&gi,&bi);
    r=ri;g=gi;b=bi;
    if(MFAtlasIsChartNearBoundary(A,chart,e))
     {
      shtric(&rnb,&gnb,&bnb);
      r=rnb;g=gnb;b=bnb;

#ifdef MFALLOWVERBOSE
      if(verbose){printf("    This one is not near the boundary, %d\n",chart);fflush(stdout);}
#endif

     }
    if(MFAtlasIsChartSingular(A,chart,e))
     {
      shtric(&rsb,&gsb,&bsb);
      r=rsb;g=gsb;b=bsb;

#ifdef MFALLOWVERBOSE
      if(verbose){printf("    This one is singular, %d\n",chart);fflush(stdout);}
#endif

     }
    for(bchart=0;bchart<MFAtlasNumberOfChartsWithBoundary(A,e);bchart++)
     {
      if(chart==MFAtlasChartWithBoundary(A,bchart,e))
       {
        shtric(&rb,&gb,&bb);
        r=rb;g=gb;b=bb;

#ifdef MFALLOWVERBOSE
        if(verbose){printf("    This one is on the boundary, %d\n",chart);fflush(stdout);}
#endif

       }
     }
    MFDrawChart(MFAtlasChart(A,chart,e),NULL,e);
    shtric(&r,&g,&b);
    if(MFAtlasK(A,e)<3)shtric(&full,&full,&zero);
    MFDrawChartBoundary(MFAtlasChart(A,chart,e),NULL,e);
   }

  return;
 }

/* ---------------------------------------------------------------------------------*/

void MFDrawMakeListOfTriangles(MFChart chart,int *ntri,int *mtri,double **ptri, double d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawMakeListOfTriangles"};
  int i,j,n;
  int in0,in1,in2;
  int allIn;
  int allOut;
  double x0,y0;
  double x1,y1;
  double x2,y2;
  double delta;
  int changed;
  int subdivide;
  int remove;
  double *tri;
  MFKVector pt0;
  MFKVector pt1;
  MFKVector pt2;

  tri=*ptri;

  pt0=MFCreateKVector(MFChartK(chart,e),e);
  pt1=MFCreateKVector(MFChartK(chart,e),e);
  pt2=MFCreateKVector(MFChartK(chart,e),e);

  changed=1;
  while(changed)
   {
    changed=0;
    n=*ntri;
    for(j=0;j<n;j++)
     {
      MFKVSetC(pt0,0,tri[  6*j],e);
      MFKVSetC(pt0,1,tri[1+6*j],e);
      MFKVSetC(pt1,0,tri[2+6*j],e);
      MFKVSetC(pt1,1,tri[3+6*j],e);
      MFKVSetC(pt2,0,tri[4+6*j],e);
      MFKVSetC(pt2,1,tri[5+6*j],e);

      in0=MFChartInterior(chart,pt0,e);
      in1=MFChartInterior(chart,pt1,e);
      in2=MFChartInterior(chart,pt2,e);

      allOut=0;
      allIn=0;
      if(in0&&in1&&in2)allIn=1;
      if(!in0&&!in1&&!in2)allOut=1;

      x0=tri[0+6*j];
      y0=tri[1+6*j];
      x1=tri[2+6*j];
      y1=tri[3+6*j];
      x2=tri[4+6*j];
      y2=tri[5+6*j];

      delta=(sqrt(pow(x1-x0,2)+pow(y1-y0,2))
            +sqrt(pow(x2-x0,2)+pow(y2-y0,2))
            +sqrt(pow(x2-x1,2)+pow(y2-y1,2)))*.3333333;

      subdivide=0;
      remove=0;
      if(delta>d)subdivide=1;
      if(!allOut&&!allIn)subdivide=1;
      if(delta<.25*d)subdivide=0;
      if(allOut)remove=1;

      if(remove && !subdivide)
       {
        changed=1;
        tri[0+6*j]=tri[0+6*(*ntri-1)];
        tri[1+6*j]=tri[1+6*(*ntri-1)];
        tri[2+6*j]=tri[2+6*(*ntri-1)];
        tri[3+6*j]=tri[3+6*(*ntri-1)];
        tri[4+6*j]=tri[4+6*(*ntri-1)];
        tri[5+6*j]=tri[5+6*(*ntri-1)];
        (*ntri)--;
        if(n>*ntri)n--;
       }

      if(subdivide)
       {
        changed=1;
        if(*ntri+3>=*mtri)
         {
          *mtri+=3*(*ntri);
          *ptri=(double*)realloc((void*)(*ptri),6*(*mtri)*sizeof(double));

#ifndef MFNOSAFETYNET
          if(*ptri==NULL)
           {
            sprintf(MFDrawErrorMsg,"Out of memory trying to allocate %d bytes.",6*(*mtri)*sizeof(double));
            MFSetError(e,12,RoutineName,MFDrawErrorMsg,__LINE__,__FILE__);
            MFErrorHandlerOutOfMemory(e);
            return;
           }
#endif

          tri=*ptri;
         }
        tri[0+6*j]=.5*(x0+x1);
        tri[1+6*j]=.5*(y0+y1);
        tri[2+6*j]=.5*(x1+x2);
        tri[3+6*j]=.5*(y1+y2);
        tri[4+6*j]=.5*(x2+x0);
        tri[5+6*j]=.5*(y2+y0);
        tri[0+6*(*ntri)]=x0;
        tri[1+6*(*ntri)]=y0;
        tri[2+6*(*ntri)]=.5*(x0+x1);
        tri[3+6*(*ntri)]=.5*(y0+y1);
        tri[4+6*(*ntri)]=.5*(x2+x0);
        tri[5+6*(*ntri)]=.5*(y2+y0);
        tri[0+6*(*ntri+1)]=x1;
        tri[1+6*(*ntri+1)]=y1;
        tri[2+6*(*ntri+1)]=.5*(x1+x2);
        tri[3+6*(*ntri+1)]=.5*(y1+y2);
        tri[4+6*(*ntri+1)]=.5*(x0+x1);
        tri[5+6*(*ntri+1)]=.5*(y0+y1);
        tri[0+6*(*ntri+2)]=x2;
        tri[1+6*(*ntri+2)]=y2;
        tri[2+6*(*ntri+2)]=.5*(x2+x0);
        tri[3+6*(*ntri+2)]=.5*(y2+y0);
        tri[4+6*(*ntri+2)]=.5*(x1+x2);
        tri[5+6*(*ntri+2)]=.5*(y1+y2);
        *ntri+=3;
       }
     }
   }

  MFFreeKVector(pt0,e);
  MFFreeKVector(pt1,e);
  MFFreeKVector(pt2,e);
  return;
 }

/* Use a project for each NSpace? */

int MFTPBVPGetNX(MFImplicitMF,MFErrorHandler);
int MFTPBVPGetNU(MFImplicitMF,MFErrorHandler);
int MFTPBVPGetNP(MFImplicitMF,MFErrorHandler);

void MFDrawProjectTo3d(MFChart chart,MFNVector u, float *x,float *y,float *z, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawProjectTo3d"};
  int i,j,n;
  int nx,nu,np;
  int verbose=0;

  n=MFNV_NC(u,e);

#ifdef MFALLOWVERBOSE
  if(verbose)printf("MFDrawProjectTo3d, n=%d\n",n);
#endif

  if(n<4)
   {

#ifdef MFALLOWVERBOSE
    if(verbose)printf("  n<4\n");
#endif

    *x=MFNV_C(u,0,e);
    *y=MFNV_C(u,1,e);
    if(n>2)
      *z=MFNV_C(u,2,e);
     else
      *z=0.;
   }else if(n==4){

#ifdef MFALLOWVERBOSE
    if(verbose)printf("  n==4\n");
#endif

    *x=.5*(MFNV_C(u,0,e)+MFNV_C(u,1,e));
    *y=MFNV_C(u,2,e);
    *z=MFNV_C(u,3,e);
   }else{

#ifdef MFALLOWVERBOSE
    if(verbose)printf("n>4, looking for a way %s\n",MFImplicitMFId(MFChartGetManifold(chart,e),e));
#endif

    if(chart!=NULL && !strcmp(MFImplicitMFId(MFChartGetManifold(chart,e),e),"Pendula"))
     {

#ifdef MFALLOWVERBOSE
      if(verbose)printf("  Pendula\n");
#endif

      *x=MFNV_C(u,n-1,e); /* Ic/(k pi) */
#ifdef PERIODIC
      if(MFPendulaPeriodicFlag)
       {
        while(*x<-1.e-7)*x+=1.;
        while(*x>.1*MFPI+1.e-7)*x-=1.;
       }
#endif
      *y=.2*MFNV_C(u,n-2,e); /* In */
      *z=.1*2.*MFNV_C(u,n-3,e); /* T  */

      *x=2*MFNV_C(u,n-1,e); /* Ic/(k pi) */
      *y=MFNV_C(u,n-2,e); /* In */
      *z=MFNV_C(u,n-3,e); /* T  */
     }else if(chart!=NULL && !strcmp(MFImplicitMFId(MFChartGetManifold(chart,e),e),"TPBVP"))
     {

#ifdef MFALLOWVERBOSE
      if(verbose)printf("  TPBVP\n");
#endif

      nu=MFTPBVPGetNU(MFChartGetManifold(chart,e),e);
      nx=MFTPBVPGetNX(MFChartGetManifold(chart,e),e);
      np=MFTPBVPGetNP(MFChartGetManifold(chart,e),e);
      if(np==2)
       {
        *x=0.;
        for(i=0;i<nx;i++)
         {
          for(j=0;j<nu;j++)
           {
            *x+=(MFNV_C(u,i*nu+j,e)+MFNV_C(u,(i+1)*nu+j,e))*(MFNV_C(u,nx*nu+np+i+1,e)-MFNV_C(u,nx*nu+np+i,e))/2.;
           }
         }
        *x=*x/nu;
        *y=MFNV_C(u,nx*nu,e);
        *z=MFNV_C(u,nx*nu+1,e);
       }else{
        *x=MFNV_C(u,nx*nu,e);
        *y=MFNV_C(u,nx*nu+1,e);
        *z=MFNV_C(u,nx*nu+2,e);
       }
     }else{

#ifdef MFALLOWVERBOSE
      if(verbose)printf("MF is not TPBVP or others\n");
#endif

      *x=MFNV_C(u,n-3,e);
      *y=MFNV_C(u,n-2,e);
      *z=MFNV_C(u,n-1,e);
     }
   }

  if(MFDrawSetExt||*x<MFDrawXMin)MFDrawXMin=*x;
  if(MFDrawSetExt||*x>MFDrawXMax)MFDrawXMax=*x;
  if(MFDrawSetExt||*y<MFDrawYMin)MFDrawYMin=*y;
  if(MFDrawSetExt||*y>MFDrawYMax)MFDrawYMax=*y;
  if(MFDrawSetExt||*z<MFDrawZMin)MFDrawZMin=*z;
  if(MFDrawSetExt||*z>MFDrawZMax)MFDrawZMax=*z;
  MFDrawSetExt=0;
/*printf("%s, ",RoutineName);MFPrintNVector(stdout,u,e);printf(" ->(%f,%f,%f)\n",x,y,z);fflush(stdout);*/

  return;
 }

void MFDrawEnumPolytope(MFEnumPolytope EP,MFChart chart,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawEnumPolytope"};
  float x0,y0,z0;
  float x1,y1,z1;
  int v[100]={0};
  int v0,v1,v2;
  int e0,e1,e2;
  int full=255;
  int zero=0;
  int i,j,l,n,k;
  int verbose=0;
  float *x,*y,*z;
  MFNVector u;
  int rc;

  int r=180;
  int g=180;
  int b=205;
  shpgc(&r,&g,&b);
  shbgc(&r,&g,&b);
  shpec(&zero,&zero,&zero);
  shbec(&zero,&zero,&zero);

  k=MFChartK(chart,e);
  n=MFChartN(chart,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("MFDrawEnumPolytope base dim =%d, target dim=%d\n",k,n);fflush(stdout);}
#endif

  if(k<1||k>3)return;

  u=MFCreateNVector(n,e);

/* Vertices */

  shlinc(&zero,&zero,&zero);
  if(0)for(i=0;i<MFEnumPolytopeNumberOfCells(EP,0,e);i++)
   {
    rc=MFChartEvaluate(chart,MFEnumPolytopeVertex(EP,i,e),u,e);
    if(rc)MFMarkPoint(chart,u,3,state,e);
   }
  
/* Edges */

  for(i=0;i<MFEnumPolytopeNumberOfCells(EP,1,e);i++)
   {
    if(MFEnumPolytopeNumberOfCellFaces(EP,1,i,e)>1)
     {
      v0=MFEnumPolytopeCellFace(EP,1,i,0,e);
      v1=MFEnumPolytopeCellFace(EP,1,i,1,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("edge %d joins vertices %d and %d\n",i,v0,v1);fflush(stdout);}
#endif

      MFDrawLineOnChart(chart,MFEnumPolytopeVertex(EP,v0,e),MFEnumPolytopeVertex(EP,v1,e),state,e);
     }
   }

/* faces */

  shtric(&full,&full,&full);
  if(k==3 && n==3)
   {
    x=NULL;
    y=NULL;
    z=NULL;
    for(i=0;i<MFEnumPolytopeNumberOfCells(EP,2,e);i++)
     {
      n=MFEnumPolytopeNumberOfCellFaces(EP,2,i,e);
  
      if(n>2)
       {
        e0=MFEnumPolytopeCellFace(EP,2,i,0,e);
        v[0]=MFEnumPolytopeCellFace(EP,1,e0,0,e);
        v[1]=MFEnumPolytopeCellFace(EP,1,e0,1,e);
        for(j=1;j<n;j++)
         {
          for(l=0;l<n;l++)
           {
            e1=MFEnumPolytopeCellFace(EP,2,i,l,e);
            if(e1!=e0&&MFEnumPolytopeCellFace(EP,1,e1,0,e)==v[j])
             {
              v[j+1]=MFEnumPolytopeCellFace(EP,1,e1,1,e);
              l=n;
             }else if(e1!=e0&&MFEnumPolytopeCellFace(EP,1,e1,1,e)==v[j])
             {
              v[j+1]=MFEnumPolytopeCellFace(EP,1,e1,0,e);
              l=n;
             }
           }
          e0=e1;
         }

        x=(float*)realloc((void*)x,n*sizeof(float));

#ifndef MFNOSAFETYNET
        if(x==NULL)
         {
          sprintf(MFDrawErrorMsg,"Out of memory trying to allocate %d bytes.",n*sizeof(float));
          MFSetError(e,12,RoutineName,MFDrawErrorMsg,__LINE__,__FILE__);
          MFErrorHandlerOutOfMemory(e);
          return;
         }
#endif

        y=(float*)realloc((void*)y,n*sizeof(float));

#ifndef MFNOSAFETYNET
        if(y==NULL)
         {
          sprintf(MFDrawErrorMsg,"Out of memory trying to allocate %d bytes.",n*sizeof(float));
          MFSetError(e,12,RoutineName,MFDrawErrorMsg,__LINE__,__FILE__);
          MFErrorHandlerOutOfMemory(e);
          return;
         }
#endif

        z=(float*)realloc((void*)z,n*sizeof(float));

#ifndef MFNOSAFETYNET
        if(z==NULL)
         {
          sprintf(MFDrawErrorMsg,"Out of memory trying to allocate %d bytes.",n*sizeof(float));
          MFSetError(e,12,RoutineName,MFDrawErrorMsg,__LINE__,__FILE__);
          MFErrorHandlerOutOfMemory(e);
          return;
         }
#endif

        rc=1;
        for(j=0;j<n;j++)
         {
          rc=rc&&MFChartEvaluate(chart,MFEnumPolytopeVertex(EP,v[j],e),u,e);
          MFDrawProjectTo3d(chart,u,x+j,y+j,z+j,e);
#ifdef PERIODIC
          if(MFPendulaPeriodicFlag&&j>0)
           {
            rc=rc&&fabs(x[j-1]-x[j])<.1;
            rc=rc&&fabs(y[j-1]-y[j])<.1;
           }
#endif
         }
#ifdef PERIODIC
        if(MFPendulaPeriodicFlag)
         {
          rc=rc&&fabs(x[n-1]-x[0])<.1;
          rc=rc&&fabs(y[n-1]-y[0])<.1;
         }
#endif
        if(rc&&MFEnumPolytopeCellIndex(EP,2,i,0,e)>7)
         {
          if(state!=NULL)MFChartStateAddPolygon(state,n,x,y,z,0,e);
            else shpg(&n,x,y,z);
         }
       }
     }
    free(x);
    free(y);
    free(z);
   }

  MFFreeNVector(u,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done MFDrawEnumPolytope\n");fflush(stdout);}
#endif

  return;
 }

void MFMarkPoint(MFChart chart, MFNVector u, int i,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawMarkPoint"};
  float x=0.;
  float y=0.;
  float z=0.;
  double xt,yt,zt;
  double xn,yn,zn;

  float x0,y0,z0;
  float x1,y1,z1;

  MFDrawProjectTo3d(chart,u,&x,&y,&z,e);

  xt=1.;
  yt=0.;
  zt=0.;

  xn=0.;
  yn=1.;
  zn=0.;

  switch(i)
   {
    case 0:  /* [ */
     x0=x-.01*xt-.02*xn;
     y0=y-.01*yt-.02*yn;
     z0=z-.01*zt-.02*zn;
     x1=x       -.02*xn;
     y1=y       -.02*yn;
     z1=z       -.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x       +.02*xn;
     y1=y       +.02*yn;
     z1=z       +.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x-.01*xt+.02*xn;
     y1=y-.01*yt+.02*yn;
     z1=z-.01*zt+.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     break;
    case 1: /* ]  */
     x0=x+.01*xt-.02*xn;
     y0=y+.01*yt-.02*yn;
     z0=z+.01*zt-.02*zn;
     x1=x       -.02*xn;
     y1=y       -.02*yn;
     z1=z       -.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x       +.02*xn;
     y1=y       +.02*yn;
     z1=z       +.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x+.01*xt+.02*xn;
     y1=y+.01*yt+.02*yn;
     z1=z+.01*zt+.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     break;
    case 2: /*  O */
     x0=x;
     y0=y;
     z0=z;
     x1=x+.02*xt-.02*xn;
     y1=y+.02*yt-.02*yn;
     z1=z+.02*zt-.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x+.02*xt+.02*xn;
     y1=y+.02*yt+.02*yn;
     z1=z+.02*zt+.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x-.02*xt+.02*xn;
     y1=y-.02*yt+.02*yn;
     z1=z-.02*zt+.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x-.02*xt-.02*xn;
     y1=y-.02*yt-.02*yn;
     z1=z-.02*zt-.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x+.02*xt-.02*xn;
     y1=y+.02*yt-.02*yn;
     z1=z+.02*zt-.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     break;
    default: /*  +  */
     x0=x       -.02*xn;
     y0=y       -.02*yn;
     z0=z       -.02*zn;
     x1=x       +.02*xn;
     y1=y       +.02*yn;
     z1=z       +.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x       -.02*xt;
     y0=y       -.02*yt;
     z0=z       -.02*zt;
     x1=x       +.02*xt;
     y1=y       +.02*yt;
     z1=z       +.02*zt;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     break;
   }
  return;
 }

void MFMarkPointOnChart(MFChart chart,MFKVector s, int ms,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawMarkPointOnChart"};
  MFNVector u;
  float x=0.;
  float y=0.;
  float z=0.;
  float xt=0.;
  float yt=0.;
  float zt=0.;
  float xn,yn,zn;

  float x0,y0,z0;
  float x1,y1,z1;

  u=MFCreateNVector(MFChartN(chart,e),e);
  if(!MFChartEvaluate(chart,s,u,e))return;
  MFDrawProjectTo3d(chart,u,&x,&y,&z,e);
  MFFreeNVector(u,e);

  u=MFMColumn(MFChartTangentSpace(chart,e),0,e);
  MFDrawProjectTo3d(chart,u,&xt,&yt,&zt,e);
  MFFreeNVector(u,e);

  xn=-yt;
  yn= xt;
  zn= 0.;

  switch(ms)
   {
    case 0:  /* [ */
     x0=x-.01*xt-.02*xn;
     y0=y-.01*yt-.02*yn;
     z0=z-.01*zt-.02*zn;
     x1=x       -.02*xn;
     y1=y       -.02*yn;
     z1=z       -.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x       +.02*xn;
     y1=y       +.02*yn;
     z1=z       +.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x-.01*xt+.02*xn;
     y1=y-.01*yt+.02*yn;
     z1=z-.01*zt+.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     break;
    case 1: /* ]  */
     x0=x+.01*xt-.02*xn;
     y0=y+.01*yt-.02*yn;
     z0=z+.01*zt-.02*zn;
     x1=x       -.02*xn;
     y1=y       -.02*yn;
     z1=z       -.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x       +.02*xn;
     y1=y       +.02*yn;
     z1=z       +.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x+.01*xt+.02*xn;
     y1=y+.01*yt+.02*yn;
     z1=z+.01*zt+.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     break;
    case 2: /*  O */
     x0=x;
     y0=y;
     z0=z;
     x1=x+.02*xt-.02*xn;
     y1=y+.02*yt-.02*yn;
     z1=z+.02*zt-.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x+.02*xt+.02*xn;
     y1=y+.02*yt+.02*yn;
     z1=z+.02*zt+.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x-.02*xt+.02*xn;
     y1=y-.02*yt+.02*yn;
     z1=z-.02*zt+.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x-.02*xt-.02*xn;
     y1=y-.02*yt-.02*yn;
     z1=z-.02*zt-.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x1; y0=y1; z0=z1;
     x1=x+.02*xt-.02*xn;
     y1=y+.02*yt-.02*yn;
     z1=z+.02*zt-.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     break;
    default: /*  +  */
     x0=x       -.02*xn;
     y0=y       -.02*yn;
     z0=z       -.02*zn;
     x1=x       +.02*xn;
     y1=y       +.02*yn;
     z1=z       +.02*zn;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     x0=x       -.02*xt;
     y0=y       -.02*yt;
     z0=z       -.02*zt;
     x1=x       +.02*xt;
     y1=y       +.02*yt;
     z1=z       +.02*zt;
     if(state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
     if(state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);
     break;
   }
  return;
 }

void MFDrawLineOnChart(MFChart chart,MFKVector s0,MFKVector s1,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawLineOnChart"};
  float x0=0.;
  float y0=0.;
  float z0=0.;
  float x1=0.;
  float y1=0.;
  float z1=0.;
  double dx,dy,dz;
  double ox,oy,oz;
  int k;
  double t,dt;
  int rc,rc0;

  if(MFTmpK==NULL)MFTmpK=MFCreateKVector(MFChartK(chart,e),e);
  if(MFTmpN==NULL)MFTmpN=MFCreateNVector(MFChartN(chart,e),e);

  k=MFChartK(chart,e);
  ox=MFKV_C(s0,0,e);dx=MFKV_C(s1,0,e)-ox;
  dy=0.;
  oy=0.;
  if(1<k){oy=MFKV_C(s0,1,e);dy=MFKV_C(s1,1,e)-oy;}
  dz=0.;
  oz=0.;
  if(2<k){oz=MFKV_C(s0,2,e);dz=MFKV_C(s1,2,e)-oz;}

  rc0=MFChartEvaluate(chart,s0,MFTmpN,e);
  MFDrawProjectTo3d(chart,MFTmpN,&x0,&y0,&z0,e);
  rc=MFChartEvaluate(chart,s1,MFTmpN,e);
  MFDrawProjectTo3d(chart,MFTmpN,&x1,&y1,&z1,e);

  dt=.05;
  if(fabs(x0-x1)+fabs(y0-y1)+fabs(z0-z1)<dt)dt=1.;
  for(t=dt;t<1.+.5*dt;t+=dt)
   {
    MFKVSetC(MFTmpK,0,ox+t*dx,e);
    if(1<k)MFKVSetC(MFTmpK,1,oy+t*dy,e);
    if(2<k)MFKVSetC(MFTmpK,2,oz+t*dz,e);
    rc=MFChartEvaluate(chart,MFTmpK,MFTmpN,e);
    MFDrawProjectTo3d(chart,MFTmpN,&x1,&y1,&z1,e);
#ifdef PERIODIC
    if(MFPendulaPeriodicFlag)
     {
      rc=rc&&fabs(x1-x0)<.1;
      rc=rc&&fabs(y1-y0)<.1;
     }
#endif

    if(rc&&rc0&&state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
    if(rc&&rc0&&state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);

    x0=x1;
    y0=y1;
    z0=z1;
    rc0=rc;
   }

  return;
 }

void MFDrawEnumDualPolytope(MFEnumDualPolytope EP, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawEnumDualPolytope"};
  float x0=0.;
  float y0=0.;
  float z0=0.;
  float x1=0.;
  float y1=0.;
  float z1=0.;
  int v[100]={0};
  int v0,v1;
  int e0,e1;
  int full=255;
  int zero=0;
  int i,n,j,l,k;
  int verbose=0;
  float x[100]={0.};
  float y[100]={0.};
  float z[100]={0.};
  float xm[100]={0.};
  float ym[100]={0.};
  float zm[100]={0.};
  MFNVector vert;

  int r=180;
  int g=180;
  int b=205;
  int rc;

  shpgc(&r,&g,&b);
  shbgc(&r,&g,&b);
  shpec(&zero,&zero,&zero);
  shbec(&zero,&zero,&zero);

  n=MFEnumDualPolytopeDimension(EP,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("MFDrawEnumDualPolytope dim =%d\n",n);fflush(stdout);}
#endif

  if(n<1||n>3)return;

/* Vertices */

  shlinc(&zero,&zero,&zero);
  if(0)
   {
    for(i=0;i<MFEnumDualPolytopeNumberOfCells(EP,0,e);i++)
     {
      MFMarkPoint(NULL,MFEnumDualPolytopeVertex(EP,i,e),3,NULL,e);
     }
   }

/* Edges */

  k=MFEnumDualPolytopeDimension(EP,e);
  if(1)for(i=0;i<MFEnumDualPolytopeNumberOfCells(EP,1,e);i++)
   {
    if(k==2 && MFEnumDualPolytopeNumberOfFaceCells(EP,1,i,e)<2)
      shlinc(&full,&zero,&zero);
     else
      shlinc(&zero,&zero,&zero);

    if(MFEnumDualPolytopeNumberOfCellFaces(EP,1,i,e)>1)
     {
      v0=MFEnumDualPolytopeCellFace(EP,1,i,0,e);
      v1=MFEnumDualPolytopeCellFace(EP,1,i,1,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("edge %d joins vertices %d and %d\n",i,v0,v1);fflush(stdout);}
#endif

      if(MFEnumDualPolytopeVertex(EP,v0,e)!=NULL
       &&MFEnumDualPolytopeVertex(EP,v1,e)!=NULL)
       {
        MFDrawProjectTo3d(NULL,MFEnumDualPolytopeVertex(EP,v0,e),&x0,&y0,&z0,e);
        MFDrawProjectTo3d(NULL,MFEnumDualPolytopeVertex(EP,v1,e),&x1,&y1,&z1,e);
        shline(&x0,&y0,&z0,&x1,&y1,&z1);
       }
     }
   }

/* faces */
  if(1)for(i=0;i<MFEnumDualPolytopeNumberOfCells(EP,2,e);i++)
   {
    n=MFEnumDualPolytopeNumberOfCellFaces(EP,2,i,e);
    if(1){

    if(n>3){printf("MFDrawEnumDualPolytope: Here's a 2-cell (%d) which has more sides than a triangle (%d)\n",i,n);fflush(stdout);}

/*  if(n>2 || ( k==3&&MFEnumDualPolytopeNumberOfFaceCells(EP,2,i,e)<2))*/
    if(n>2)
     {
      e0=MFEnumDualPolytopeCellFace(EP,2,i,0,e);
      v[0]=MFEnumDualPolytopeCellFace(EP,1,e0,0,e);
      v[1]=MFEnumDualPolytopeCellFace(EP,1,e0,1,e);
      if(n>4)printf("Edge %d is %d [%d,%d]\n",i,e0,v[0],v[1]);
      for(j=1;j<n;j++)
       {
        if(n>4)printf("Edge %d is %d [%d,%d]\n",i,MFEnumDualPolytopeCellFace(EP,2,i,j,e),MFEnumDualPolytopeCellFace(EP,1,MFEnumDualPolytopeCellFace(EP,2,i,j,e),0,e),MFEnumDualPolytopeCellFace(EP,1,MFEnumDualPolytopeCellFace(EP,2,i,j,e),1,e));
        for(l=0;l<n;l++)
         {
          e1=MFEnumDualPolytopeCellFace(EP,2,i,l,e);
/*        if(MFEnumDualPolytopeNumberOfFaceCells(EP,1,e1,e)>1)*/
           {
            if(e1!=e0&&MFEnumDualPolytopeCellFace(EP,1,e1,0,e)==v[j])
             {
              v[j+1]=MFEnumDualPolytopeCellFace(EP,1,e1,1,e);
              l=n;
             }else if(e1!=e0&&MFEnumDualPolytopeCellFace(EP,1,e1,1,e)==v[j])
             {
              v[j+1]=MFEnumDualPolytopeCellFace(EP,1,e1,0,e);
              l=n;
             }
           }
         }
        e0=e1;
       }

      for(j=0;j<n;j++)
       {
        vert=MFEnumDualPolytopeVertex(EP,v[j],e);
        if(vert!=NULL)
         {
          x[j]=0;
          y[j]=0;
          z[j]=0;
          if(MFNV_NC(vert,e)>0)
            MFDrawProjectTo3d(NULL,MFEnumDualPolytopeVertex(EP,v[j],e),x+j,y+j,z+j,e);
         }
       }
      x[n]=x[0];
      y[n]=y[0];
      z[n]=z[0];

      rc=1;
      for(j=0;j<n;j++)
       {
        xm[j]=.5*(x[j]+x[j+1]);
        ym[j]=.5*(y[j]+y[j+1]);
        zm[j]=.5*(z[j]+z[j+1]);
#ifdef PERIODIC
        if(MFPendulaPeriodicFlag&&j>0)
         {
          rc=rc&&fabs(x[j-1]-x[j])<.1;
          rc=rc&&fabs(y[j-1]-y[j])<.1;
         }
#endif
       }
#ifdef PERIODIC
      if(MFPendulaPeriodicFlag)
       {
        rc=rc&&fabs(x[n-1]-x[0])<.1;
        rc=rc&&fabs(y[n-1]-y[0])<.1;
       }
#endif

      if(rc&&n<5)
       {
        if(n<4)shtric(&r,&g,&b);
         else  shtric(&full,&full,&zero);
        if(k==3 && MFEnumDualPolytopeNumberOfFaceCells(EP,2,i,e)<2)
         {
          shtric(&full,&zero,&zero);
          shpec(&zero,&zero,&zero);
          shbec(&zero,&zero,&zero);
         }else{
          shtric(&r,&g,&b);
          shpec(&zero,&zero,&zero);
          shbec(&zero,&zero,&zero);
         }
        if(1)shpg(&n,x,y,z);
         else shpg(&n,xm,ym,zm);
       }else if(rc)
       {
        shlinc(&zero,&full,&zero);
        shpl(&n,x,y,z);
       }
     }    /* if(n>2 && MFEnumDualPolytopeNumberOfFaceCells(EP,2,i,e)<2) */
   }      /* if(FALSE) */
   }      /* for(i=0,..,MFEnumDualPolytopeNumberOfCells(EP,2,e)) */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done MFDrawEnumDualPolytope\n");fflush(stdout);}
#endif

  return;
 }

void MFDrawLine(MFNVector u0,MFNVector u1,MFChartState state, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawLine"};
  float x0=0.;
  float y0=0.;
  float z0=0.;
  float x1=0.;
  float y1=0.;
  float z1=0.;
  int k;
  int rc;

  MFDrawProjectTo3d(NULL,u0,&x0,&y0,&z0,e);
  MFDrawProjectTo3d(NULL,u1,&x1,&y1,&z1,e);

  rc=1;
#ifdef PERIODIC
  if(MFPendulaPeriodicFlag)
   {
    rc=rc&&fabs(x0-x1)<.1;
    rc=rc&&fabs(y0-y1)<.1;
   }
#endif

  if(rc&&state!=NULL)MFChartStateAddLine(state,x0,y0,z0,x1,y1,z1,e);
  if(rc&&state==NULL)shline(&x0,&y0,&z0,&x1,&y1,&z1);

  return;
 }

void MFDrawChartTS(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawChartTS"};

  if(chart==NULL)return;
  if(MFChartK(chart,e)==1)MFDraw1dChartTS(chart,e);
   else if(MFChartK(chart,e)==2)MFDraw2dChartTS(chart,e);
/* else if(MFChartK(chart,e)==3)MFDraw3dChartTS(chart,e);*/

  return;
 }

void MFDrawAtlasTS(MFAtlas A, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawAtlasTS"};
  int chart,bchart;
/*int ri=180;
  int gi=180;
  int bi=205;*/
  int ri=62;
  int gi=158;
  int bi=255;
  int zero=0;
  int full=255;
/*int rb=243;
  int gb=223;
  int bb=192;*/
  int rb=255;
  int gb=211;
  int bb=62;
  int N;
  int ibnd;
  int verbose=0;

  MFDrawSetExt=1;
  MFDrawXMin=0.;
  MFDrawXMax=0.;
  MFDrawYMin=0.;
  MFDrawYMax=0.;
  MFDrawZMin=0.;
  MFDrawZMax=0.;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("In %s\n",RoutineName);fflush(stdout);}
#endif

  shlinc(&zero,&zero,&zero);
  N=MFAtlasNumberOfCharts(A,e);
  for(chart=0;chart<N;chart++)
   {
    ibnd=0;
    shtric(&ri,&gi,&bi);
    for(bchart=0;bchart<MFAtlasNumberOfChartsWithBoundary(A,e);bchart++)
     {
      if(chart==MFAtlasChartWithBoundary(A,bchart,e))
       {
        shtric(&rb,&gb,&bb);
        ibnd=1;
       }
     }
    MFDrawChartTS(MFAtlasChart(A,chart,e),e);
    MFDrawPolytopeTS(MFChartPolytope(MFAtlasChart(A,chart,e),e),MFAtlasChart(A,chart,e),e);
    if(ibnd)
      MFDrawChartBoundaryTS(MFAtlasChart(A,chart,e),e);

#ifdef MFALLOWVERBOSE
    if(verbose&&chart%1000==0&&chart!=0){printf("     %d/%d\n",chart,N);fflush(stdout);}
#endif

   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  printf("at end of MFDrawAtlas, x in ([%lf,%lf],[%lf,%lf],[%lf,%lf])\n",MFDrawXMin,MFDrawXMax,MFDrawYMin,MFDrawYMax,MFDrawZMin,MFDrawZMax);fflush(stdout);
  return;
 }

void MFDraw1dChartTS(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDraw1dChartTS"};
  int color;
  int zero=0;
  int half=128;
  int full=255;
  MFPolytope P;
  MFKVector s0,s1;
  MFNVector u0,u1;
  int n,k;

  P=MFChartPolytope(chart,e);
  k=MFChartK(chart,e);
  n=MFChartN(chart,e);
  if(MFPolytopeNumberOfVertices(P,e)<2)return;
  s0=MFCreateKVector(k,e);
  s1=MFCreateKVector(k,e);
  u0=MFCreateNVector(n,e);
  u1=MFCreateNVector(n,e);

  MFPolytopeVertex(P,0,s0,e);
  MFChartPointInTangentSpace(chart,s0,u0,e);
  MFPolytopeVertex(P,1,s1,e);
  MFChartPointInTangentSpace(chart,s1,u1,e);

  MFMarkPoint(NULL,u0,2,NULL,e);
  MFMarkPoint(NULL,u1,2,NULL,e);
  MFMarkPoint(NULL,MFChartCenter(chart,e),3,NULL,e);

  MFDrawLine(u0,u1,NULL,e);

  MFFreeKVector(s0,e);
  MFFreeKVector(s1,e);
  MFFreeNVector(u0,e);
  MFFreeNVector(u1,e);

  return;
 }

void MFDraw2dChartTS(MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDraw2dChartTS"};
  int color;
  int zero=0;
  int half=128;
  int full=255;
  double *tri;
  int ntri,mtri;
  static float x[3],y[3],z[3];
  double t;
  MFKVector s;
  MFNVector u;
  int j,n;
  int rc;

/* The interior of the circle */

  tri=(double*)malloc(6*6*sizeof(double));

#ifndef MFNOSAFETYNET
  if(tri==NULL)
   {
    sprintf(MFDrawErrorMsg,"Out of memory trying to allocate %d bytes.",6*6*sizeof(double));
    MFSetError(e,12,RoutineName,MFDrawErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  t=2./sqrt3*MFChartRadius(chart,e);

  tri[0+6*0]=0.;
  tri[1+6*0]=0.;
  tri[2+6*0]=t;
  tri[3+6*0]=0.;
  tri[4+6*0]=.5*t;
  tri[5+6*0]=.5*sqrt3*t;

  tri[0+6*1]=0.;
  tri[1+6*1]=0.;
  tri[2+6*1]=.5*t;
  tri[3+6*1]=.5*sqrt3*t;
  tri[4+6*1]=-.5*t;
  tri[5+6*1]=.5*sqrt3*t;

  tri[0+6*2]=0.;
  tri[1+6*2]=0.;
  tri[2+6*2]=-.5*t;
  tri[3+6*2]=.5*sqrt3*t;
  tri[4+6*2]=-t;
  tri[5+6*2]=0.;

  tri[0+6*3]=0.;
  tri[1+6*3]=0.;
  tri[2+6*3]=-t;
  tri[3+6*3]=0.;
  tri[4+6*3]=-.5*t;
  tri[5+6*3]=-.5*sqrt3*t;

  tri[0+6*4]=0.;
  tri[1+6*4]=0.;
  tri[2+6*4]=-.5*t;
  tri[3+6*4]=-.5*sqrt3*t;
  tri[4+6*4]= .5*t;
  tri[5+6*4]=-.5*sqrt3*t;

  tri[0+6*5]=0.;
  tri[1+6*5]=0.;
  tri[2+6*5]= .5*t;
  tri[3+6*5]=-.5*sqrt3*t;
  tri[4+6*5]=t;
  tri[5+6*5]=0.;

  ntri=6;
  mtri=6;
  t=.2*MFChartRadius(chart,e);
/*if(t<.01)t=.01;*/
  MFDrawMakeListOfTriangles(chart,&ntri,&mtri,&tri,t,e);

  s=MFCreateKVector(MFChartK(chart,e),e);
  u=MFCreateNVector(MFChartN(chart,e),e);
  for(j=0;j<ntri;j++)
   {
    MFKVSetC(s,0,tri[0+6*j],e);
    MFKVSetC(s,1,tri[1+6*j],e);
    MFChartPointInTangentSpace(chart,s,u,e);
    MFDrawProjectTo3d(chart,u,x,y,z,e);
    MFKVSetC(s,0,tri[2+6*j],e);
    MFKVSetC(s,1,tri[3+6*j],e);
    MFChartPointInTangentSpace(chart,s,u,e);
    MFDrawProjectTo3d(chart,u,x+1,y+1,z+1,e);
    MFKVSetC(s,0,tri[4+6*j],e);
    MFKVSetC(s,1,tri[5+6*j],e);
    MFChartPointInTangentSpace(chart,s,u,e);
    MFDrawProjectTo3d(chart,u,x+2,y+2,z+2,e);
    rc=1;
#ifdef PERIODIC
    if(MFPendulaPeriodicFlag)
     {
      rc=rc&&fabs(x[0]-x[1])<.1;
      rc=rc&&fabs(x[1]-x[2])<.1;
      rc=rc&&fabs(x[2]-x[0])<.1;
      rc=rc&&fabs(y[0]-y[1])<.1;
      rc=rc&&fabs(y[1]-y[2])<.1;
      rc=rc&&fabs(y[2]-y[0])<.1;
     }
#endif
    n=3;
    if(rc)shpg(&n,x,y,z);
   }
  free(tri);
  MFFreeKVector(s,e);
  MFFreeNVector(u,e);

  return;
 }

void MFDrawEnumPolytopeTS(MFEnumPolytope EP,MFChart chart, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawEnumPolytopeTS"};
  float x0,y0,z0;
  float x1,y1,z1;
  int v[100]={0};
  int v0,v1,v2;
  int e0,e1,e2;
  int full=255;
  int zero=0;
  int i,j,l,n,k;
  int verbose=0;
  float *x,*y,*z;
  MFNVector u;
  int rc;

  int r=180;
  int g=180;
  int b=205;
  shpgc(&r,&g,&b);
  shbgc(&r,&g,&b);
  shpec(&zero,&zero,&zero);
  shbec(&zero,&zero,&zero);

  k=MFChartK(chart,e);
  n=MFChartN(chart,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("MFDrawEnumPolytope base dim =%d, target dim=%d\n",k,n);fflush(stdout);}
#endif

  if(k<1||k>3)return;

  u=MFCreateNVector(n,e);

/* Vertices */

  shlinc(&zero,&zero,&zero);
  if(0)
   {
    for(i=0;i<MFEnumPolytopeNumberOfCells(EP,0,e);i++)
     {
      MFChartPointInTangentSpace(chart,MFEnumPolytopeVertex(EP,i,e),u,e);
      MFMarkPoint(chart,u,3,NULL,e);
     }
   }
  
/* Edges */

  for(i=0;i<MFEnumPolytopeNumberOfCells(EP,1,e);i++)
   {
    if(MFEnumPolytopeNumberOfCellFaces(EP,1,i,e)>1)
     {
      v0=MFEnumPolytopeCellFace(EP,1,i,0,e);
      v1=MFEnumPolytopeCellFace(EP,1,i,1,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("edge %d joins vertices %d and %d\n",i,v0,v1);fflush(stdout);}
#endif
      MFDrawLineOnChartTS(chart,MFEnumPolytopeVertex(EP,v0,e),MFEnumPolytopeVertex(EP,v1,e),e);
     }
   }

/* faces */

  if(k==3 && n==3)
   {
    x=NULL;
    y=NULL;
    z=NULL;
    for(i=0;i<MFEnumPolytopeNumberOfCells(EP,2,e);i++)
     {
      n=MFEnumPolytopeNumberOfCellFaces(EP,2,i,e);
  
      if(n>2)
       {
        e0=MFEnumPolytopeCellFace(EP,2,i,0,e);
        v[0]=MFEnumPolytopeCellFace(EP,1,e0,0,e);
        v[1]=MFEnumPolytopeCellFace(EP,1,e0,1,e);
        for(j=1;j<n;j++)
         {
          for(l=0;l<n;l++)
           {
            e1=MFEnumPolytopeCellFace(EP,2,i,l,e);
            if(e1!=e0&&MFEnumPolytopeCellFace(EP,1,e1,0,e)==v[j])
             {
              v[j+1]=MFEnumPolytopeCellFace(EP,1,e1,1,e);
              l=n;
             }else if(e1!=e0&&MFEnumPolytopeCellFace(EP,1,e1,1,e)==v[j])
             {
              v[j+1]=MFEnumPolytopeCellFace(EP,1,e1,0,e);
              l=n;
             }
           }
          e0=e1;
         }

        x=(float*)realloc((void*)x,n*sizeof(float));

#ifndef MFNOSAFETYNET
        if(x==NULL)
         {
          sprintf(MFDrawErrorMsg,"Out of memory trying to allocate %d bytes.",n*sizeof(float));
          MFSetError(e,12,RoutineName,MFDrawErrorMsg,__LINE__,__FILE__);
          MFErrorHandlerOutOfMemory(e);
          return;
         }
#endif

        y=(float*)realloc((void*)y,n*sizeof(float));

#ifndef MFNOSAFETYNET
        if(y==NULL)
         {
          sprintf(MFDrawErrorMsg,"Out of memory trying to allocate %d bytes.",n*sizeof(float));
          MFSetError(e,12,RoutineName,MFDrawErrorMsg,__LINE__,__FILE__);
          MFErrorHandlerOutOfMemory(e);
          return;
         }
#endif

        z=(float*)realloc((void*)z,n*sizeof(float));

#ifndef MFNOSAFETYNET
        if(z==NULL)
         {
          sprintf(MFDrawErrorMsg,"Out of memory trying to allocate %d bytes.",n*sizeof(float));
          MFSetError(e,12,RoutineName,MFDrawErrorMsg,__LINE__,__FILE__);
          MFErrorHandlerOutOfMemory(e);
          return;
         }
#endif

        rc=1;
        for(j=0;j<n;j++)
         {
          MFChartPointInTangentSpace(chart,MFEnumPolytopeVertex(EP,v[j],e),u,e);
          MFDrawProjectTo3d(chart,u,x+j,y+j,z+j,e);
#ifdef PERIODIC
          if(MFPendulaPeriodicFlag&&j>0)
           {
            rc=rc&&fabs(x[j-1]-x[j])<.1;
            rc=rc&&fabs(y[j-1]-y[j])<.1;
           }
#endif
         }
#ifdef PERIODIC
        if(MFPendulaPeriodicFlag)
         {
          rc=rc&&fabs(x[n-1]-x[0])<.1;
          rc=rc&&fabs(y[n-1]-y[0])<.1;
         }
#endif
        if(rc&&MFEnumPolytopeCellIndex(EP,2,i,0,e)>7)shpg(&n,x,y,z);
       }
     }
    free(x);
    free(y);
    free(z);
   }

  MFFreeNVector(u,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return;
 }

void MFDrawLineOnChartTS(MFChart chart,MFKVector s0,MFKVector s1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawLineOnChartTS"};
  float x0=0.;
  float y0=0.;
  float z0=0.;
  float x1=0.;
  float y1=0.;
  float z1=0.;
  int rc;

  if(MFTmpK==NULL)MFTmpK=MFCreateKVector(MFChartK(chart,e),e);
  if(MFTmpN==NULL)MFTmpN=MFCreateNVector(MFChartN(chart,e),e);

  MFChartPointInTangentSpace(chart,s0,MFTmpN,e);
  MFDrawProjectTo3d(chart,MFTmpN,&x0,&y0,&z0,e);
  MFChartPointInTangentSpace(chart,s1,MFTmpN,e);
  MFDrawProjectTo3d(chart,MFTmpN,&x1,&y1,&z1,e);

  rc=1;
#ifdef PERIODIC
  if(MFPendulaPeriodicFlag)
   {
    rc=rc&&fabs(x0-x1)<.1;
    rc=rc&&fabs(y0-y1)<.1;
   }
#endif

  if(rc)shline(&x0,&y0,&z0,&x1,&y1,&z1);

  return;
 }

void MFPendulaPeriodic(int flag, MFErrorHandler e)
 {
  MFPendulaPeriodicFlag=flag;
 }

void MFDrawEnumDualPolytopeEdges(MFEnumDualPolytope EP, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawEnumDualPolytopeEdges"};
  float x0=0.;
  float y0=0.;
  float z0=0.;
  float x1=0.;
  float y1=0.;
  float z1=0.;
  int v[100]={0};
  int v0,v1;
  int e0,e1;
  int full=255;
  int zero=0;
  int i,n,j,l,k;
  int verbose=0;
  float x[100]={0.};
  float y[100]={0.};
  float z[100]={0.};
  float xm[100]={0.};
  float ym[100]={0.};
  float zm[100]={0.};
  MFNVector vert;

  int r=180;
  int g=180;
  int b=205;
  int rc;

  shpgc(&r,&g,&b);
  shbgc(&r,&g,&b);
  shpec(&zero,&zero,&zero);
  shbec(&zero,&zero,&zero);

  n=MFEnumDualPolytopeDimension(EP,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("MFDrawEnumDualPolytope dim =%d\n",n);fflush(stdout);}
#endif

  if(n<1||n>3)return;

/* Vertices */

  shlinc(&zero,&zero,&zero);
  if(0)
   {
    for(i=0;i<MFEnumDualPolytopeNumberOfCells(EP,0,e);i++)
     {
      MFMarkPoint(NULL,MFEnumDualPolytopeVertex(EP,i,e),3,NULL,e);
     }
   }

/* Edges */

  k=MFEnumDualPolytopeDimension(EP,e);
  if(1)for(i=0;i<MFEnumDualPolytopeNumberOfCells(EP,1,e);i++)
   {
    if(k==2 && MFEnumDualPolytopeNumberOfFaceCells(EP,1,i,e)<2)
      shlinc(&full,&zero,&zero);
     else
      shlinc(&zero,&zero,&zero);

    if(MFEnumDualPolytopeNumberOfCellFaces(EP,1,i,e)>1)
     {
      v0=MFEnumDualPolytopeCellFace(EP,1,i,0,e);
      v1=MFEnumDualPolytopeCellFace(EP,1,i,1,e);

#ifdef MFALLOWVERBOSE
      if(verbose){printf("edge %d joins vertices %d and %d\n",i,v0,v1);fflush(stdout);}
#endif

      if(MFEnumDualPolytopeVertex(EP,v0,e)!=NULL
       &&MFEnumDualPolytopeVertex(EP,v1,e)!=NULL)
       {
        MFDrawProjectTo3d(NULL,MFEnumDualPolytopeVertex(EP,v0,e),&x0,&y0,&z0,e);
        MFDrawProjectTo3d(NULL,MFEnumDualPolytopeVertex(EP,v1,e),&x1,&y1,&z1,e);
        shline(&x0,&y0,&z0,&x1,&y1,&z1);
       }
     }
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return;
 }

#include <MFFortran.h>
double K(double);

int MFDrawGetData(MFChart chart,MFNVector u, double *l, MFErrorHandler e)
 {
  static char RoutineName[]={"MFDrawGetData"};
  int i,j,n;
  int nx,nu,np;
  double a,c,t,x1,x2,x3,y1,y2,y3;
  double kp,kp0;
  int min,max;
  double d1x,d1y,d1z;
  double mx0,mz,f,m3,theta0;
  double delta,up,um,u0;
  double m,T,Omega,R,Wr;

  if(l==NULL)
   {
    if(chart!=NULL && !strcmp(MFImplicitMFId(MFChartGetManifold(chart,e),e),"TPBVP"))return 6;
     else return 0;
   }

  n=MFNV_NC(u,e);

#ifdef HAVE_DRF
  if(chart!=NULL && !strcmp(MFImplicitMFId(MFChartGetManifold(chart,e)),"TPBVP"))
   {
    nu=MFTPBVPGetNU(MFChartGetManifold(chart,e),e);
    nx=MFTPBVPGetNX(MFChartGetManifold(chart,e),e);
    np=MFTPBVPGetNP(MFChartGetManifold(chart,e),e);
    min=0;max=0;
    a=MFNV_C(u,nx*nu+0,e);
    c=MFNV_C(u,nx*nu+1,e);
    t=MFNV_C(u,nx*nu+2,e);
    x1=(MFNV_C(u,0,e)+MFNV_C(u,nu+0,e))/2;
    y1=(MFNV_C(u,1,e)+MFNV_C(u,nu+1,e))/2;
    x2=(MFNV_C(u,2,e)+MFNV_C(u,nu+2,e))/2;
    y2=(MFNV_C(u,3,e)+MFNV_C(u,nu+3,e))/2;
    x3=(MFNV_C(u,4,e)+MFNV_C(u,nu+4,e))/2;
    y3=(MFNV_C(u,5,e)+MFNV_C(u,nu+5,e))/2;
    kp0=t*x1*y1+(t*x2-c)*y2;

    mx0=c;
    mz=a;
    f=t;
    theta0=MFNV_C(u,nx*nu+3,e);
    m3=mx0*y1+mz*y3;

    u0=cos(theta0);
    delta=sqrt((mz*mz+mx0*mx0)*(mz*mz+mx0*mx0)+8*f*(2*f+u0*(mx0*mx0-mz*mz))
                 -16*mz*mx0*f*sin(theta0));
    up=(mz*mz+mx0*mx0+delta)/4/f;
    um=(mz*mz+mx0*mx0-delta)/4/f;

    nu=up-u0;
    if(um<u0)nu=up-um;
    m=fabs(u0-um)/nu;
    T=2*sqrt(2.)*K(m)/sqrt(f*nu);
    Omega=2*3.1415926/T;
    R=0.;
    Wr=0.;

    d1x=0.;
    d1y=1.;
    d1z=0.;
    for(i=1;i<nx-1;i++)
     {
      x1=(MFNV_C(u,nu*i+0,e)+MFNV_C(u,nu*(i+1)+0,e))/2;
      y1=(MFNV_C(u,nu*i+1,e)+MFNV_C(u,nu*(i+1)+1,e))/2;
      x2=(MFNV_C(u,nu*i+2,e)+MFNV_C(u,nu*(i+1)+2,e))/2;
      y2=(MFNV_C(u,nu*i+3,e)+MFNV_C(u,nu*(i+1)+3,e))/2;
      x3=(MFNV_C(u,nu*i+4,e)+MFNV_C(u,nu*(i+1)+4,e))/2;
      y3=(MFNV_C(u,nu*i+5,e)+MFNV_C(u,nu*(i+1)+5,e))/2;
      kp=t*x1*y1+(t*x2-c)*y2;
      if(kp*kp0<=0.)
       {
        if(i==1)
         {
          if(kp>0.)min+=1;
          if(kp<0.)max+=1;
         }else{
          if(kp>0.)min+=2;
          if(kp<0.)max+=2;
         }
       }
      kp0=kp;
     }
    l[0]=theta0;
    l[1]=x1*y1+x2*y2+x3*y3;
    l[2]=min;
    l[3]=max;
    l[4]=Omega;
    l[5]=y3;
/*  l[6]=Wr;
    l[7]=R;*/
   }else if(chart!=NULL)
#endif
   {
    l[0]=MFChartHasBoundary(chart,e);
    l[1]=MFChartIsSingular(chart,e);
    l[2]=0.;
   }


  return 0;
 }

#ifdef HAVE_DRF
#ifdef F77_FUNC
#define Call_DRF F77_FUNC(drf,DRF)
#else
#define Call_DRF FC_FUNC(drf,DRF)
#endif

double K(double x, MFErrorHandler e)
 {
  double zero=0.;
  double one=1.;
  int ier;

/*x=1-x;*/
  x=1-x*x;
  return Call_DRF(&zero,&x,&one,&ier);
 }
#endif

#ifdef __cplusplus
}
#endif
