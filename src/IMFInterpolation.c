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

static char *id="@(#) $Id: IMFInterpolation.c,v 1.11 2011/07/21 17:42:46 mhender Exp $";

static char IMFInterpolationErrorMsg[256];

#include <MFAtlas.h>
#include <MFChart.h>
#include <MFFortran.h>
#include <IMFFlow.h>
#include <IMFExpansion.h>
#include <IMFExpansionSpace.h>
#include <IMFExpansionPt.h>
#include <IMFSphereOnExpansion.h>
#include <IMFInterpolation.h>
#include <MFImplicitMF.h>
#include <MFErrorHandler.h>
#include <MFPrint.h>
#include <IMF.h>
#include <math.h>
#include <stdio.h>

#ifdef __cplusplus
 extern "C" {
#endif

/* #define DOINTERPANIM*/
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

FILE *IMFInterper=NULL;
int IMFNInterper=0;
FILE *IMFInterpee=NULL;
int IMFNInterpee=0;
FILE *IMFInterpT=NULL;
int IMFNInterpT=0;
FILE *IMFCircle=NULL;
int IMFNCircle=0;

static MFKVector IMFFlowP0;

MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);
int MFPolytopeVertexNumberOfIndices(MFPolytope,int,MFErrorHandler);

static MFNVector IMFInterpolateOnDualFace(MFAtlas,IMFFlow,MFKVector,MFNVector,int,int*,MFErrorHandler);
static MFNVector IMFInterpolateOnDualFace2(MFAtlas,IMFFlow,MFKVector,MFAtlas,MFKVector,MFNVector,int,int*,MFErrorHandler);
static void MFNullVectorSmall(int,double*,double*,MFErrorHandler);
static void MFSolveFullSmall(int,double*,double*,MFErrorHandler);
double *MFKV_CStar(MFKVector,MFErrorHandler);
double *MFNV_CStar(MFNVector,MFErrorHandler);
int MFChartChangedFlag(MFChart,MFErrorHandler);
void MFAtlasRemoveChartFromBoundaryList(MFAtlas,int,MFErrorHandler);
void MFChartResetChangedFlag(MFChart,MFErrorHandler);
int MFPolytopeIntersectIndexSets(MFPolytope,int,int,int*,MFErrorHandler);
static void MFSolveFullSmall(int,double*,double*,MFErrorHandler);
int  MFSolveFull(int,double*,double*,MFErrorHandler);
MFNVector IMFGetInterpolationPointOnList(MFAtlas,IMFFlow,MFKVector,MFAtlas,double,MFNRegion,int,int*,MFErrorHandler);
int MFAtlasHalfSpaceLeftChart(MFAtlas,int,MFErrorHandler);
int MFAtlasHalfSpaceRightChart(MFAtlas,int,MFErrorHandler);

int MFTPBVPGetNX(MFImplicitMF, MFErrorHandler);
int MFTPBVPGetNU(MFImplicitMF, MFErrorHandler);
int MFTPBVPGetNP(MFImplicitMF, MFErrorHandler);
int MFTPBVPGetNIC(MFImplicitMF, MFErrorHandler);
int MFTPBVPGetNBC(MFImplicitMF, MFErrorHandler);
int IMFExpansionNVGetPrev(MFNVector,MFErrorHandler);
void IMFExpansionNVSetPrev(MFNVector,int,MFErrorHandler);

#define kmax 10
static void TestEm(int,int,int,MFErrorHandler);
static void TestF(double r, int nu, double *u, int np, double *p, double *u0,double *l0,
                               void (*f)(double,int,double*,int,double*,double*,double*,double*,MFErrorHandler),
                               void (*fu)(double,int,double*,int,double*,double*,double*,double*,MFErrorHandler),
                               void (*fl)(double,int,double*,int,double*,double*,double*,double*,MFErrorHandler),MFErrorHandler);
static void TestA(int nbc, int nu, double *uL, double *uR, int np, double *p, double *u0L, double *u0R,double *l0,
                               void (*a)(int,int,double*,double*,int,double*,double*,double*,double*,double*,MFErrorHandler),
                               void (*au)(int,int,double*,double*,int,double*,double*,double*,double*,double*,MFErrorHandler),
                               void (*al)(int,int,double*,double*,int,double*,double*,double*,double*,double*,MFErrorHandler),MFErrorHandler);

MFNVector IMFGetInterpolationPointOnList(MFAtlas,IMFFlow,MFKVector,MFAtlas,double,MFNRegion,int,int*,MFErrorHandler);

/*! \fn MFNVector IMFGetInterpolationPoint(MFAtlas A,IMFFlow F, MFKVector p0, MFAtlas c, double tmax, MFNRegion Omega, MFErrorHandler e)
 *  \brief Finds a point on the boundary of the current charts in the atlas at which the flow is outward.
 *
 *  \param A The atlas.
 *  \param F The flow.
 *  \param p0 The parameters for the flow.
 *  \param c The manifold of initial conditions. Points on this have first priority.
 *  \param tmax The largest time allowed on the manifold.
 *  \param Omega A region which limits the computation.
 *  \param e An error handler.
 *  \returns A new point satisfying the requirements.
 */
MFNVector IMFGetInterpolationPoint(MFAtlas A,IMFFlow F, MFKVector p0, MFAtlas c, double tmax, MFNRegion Omega, MFErrorHandler e)
 {
/*
    edges whose index is all >m give points where edge crosses sphere
    clip points against all faces with index >m not in index set of edge
    mark intersection point
    index of edge is ind1,ind2
*/

  static char RoutineName[]={"IMFGetInterpolationPoint"};
  int verbose=0;

  MFPolytope P;
  int iBndChart,nBndCharts;
  int chart;
  int doit;

  int i,j,k,l;
  int n,m;
  int n1,n2;
  double d0,d1,t;
  double r;
  int good0;
  int good1;
  int inter[kmax];
  MFNVector result;
  double R;

  int ii,jj;
  int noncube;
  double a[(kmax-1)*kmax],b[kmax-1],aa[(kmax-1)*(kmax-1)],asave[(kmax-1)*kmax];
  int iface[kmax-1];
  double d[kmax];
  double org[kmax],nrm[kmax];
  double oo;
  double *nn;
  double *f;
  MFKVector vf;
  MFNVector z;
  MFNVector vp;
  double *p;
  MFNVector vF;
  double *pF;
  int nv;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

#ifdef DOINTERPANIM
  if(IMFInterper!=NULL)
   {
    fclose(IMFInterper);
    IMFInterper=NULL;
   }
  if(IMFInterper==NULL)
   {
    IMFInterper=fopen("Interp1.dx","w");
    IMFNInterper=0;
   }
  if(IMFInterpee==NULL)
   {
    IMFInterpee=fopen("Interp2.dx","w");
    IMFNInterpee=0;
   }
  if(IMFInterpT==NULL)
   {
    IMFInterpT=fopen("Interp3.dx","w");
    IMFNInterpT=0;
   }
  if(IMFCircle!=NULL)
   {
    fclose(IMFCircle);
    IMFCircle=NULL;
   }
  if(IMFCircle==NULL)
   {
    IMFCircle=fopen("IMFCircle.dx","w");
    IMFNCircle=0;
   }
#endif

  nBndCharts=MFAtlasNumberOfChartsWithBoundary(A,e);
  if(verbose){printf("There are %d charts on the boundary list\n",nBndCharts);fflush(stdout);}

  k=MFAtlasK(A,e);
  if(k>kmax)
   {
    fprintf(stderr,"%s, I statically allocated some things here for k<10! Complain.\n");
    fflush(stderr);
    abort();
   }

  m=1;
  for(i=0;i<k;i++)m=m*2;
  m=m-1;

  vf=MFCreateKVector(k,e);
  f=MFKV_CStar(vf,e);
  z=MFCreateNVector(MFAtlasN(A,e),e);
  vp=MFCreateNVector(MFAtlasN(A,e),e);
  p=MFNV_CStar(vp,e);
  vF=MFCreateNVector(MFAtlasN(A,e),e);
  pF=MFNV_CStar(vF,e);

/* Clean boundary list */

  for(iBndChart=0;iBndChart<nBndCharts;iBndChart++)
   {
    chart=MFAtlasChartWithBoundary(A,iBndChart,e);
    P=MFChartPolytope(MFAtlasChart(A,chart,e),e);

    d0=0.;
    for(i=0;i<MFPolytopeNumberOfVertices(P,e);i++)
     {
      d1=MFPolytopeRadiusOfVertex(MFChartPolytope(MFAtlasChart(A,chart,e),e),i,e);
      if(d1>d0)d0=d1;
     }

    R=MFChartRadius(MFAtlasChart(A,chart,e),e);
    if(MFChartIsSingular(MFAtlasChart(A,chart,e),e)
           ||d0<R
           ||IMFExpansionNVGetT(MFChartCenter(MFAtlasChart(A,chart,e),e),e)>tmax
           ||(Omega!=NULL&&!MFNRegionInterior(Omega,IMFExpansionNVGetSigma(MFChartCenter(MFAtlasChart(A,chart,e),e),e),e)) )
     {
      MFAtlasRemoveChartFromBoundaryList(A,iBndChart,e);
      iBndChart--;
      nBndCharts--;
     }
   }

#ifdef DOINTERPANIM
  nBndCharts=MFAtlasNumberOfChartsWithBoundary(A,e);
  nBndCharts=MFAtlasNumberOfCharts(A,e);
  for(iBndChart=0;iBndChart<nBndCharts;iBndChart++)
   {
    chart=MFAtlasChartWithBoundary(A,iBndChart,e);
    chart=iBndChart;
    P=MFChartPolytope(MFAtlasChart(A,chart,e),e);
    R=MFChartRadius(MFAtlasChart(A,chart,e),e);

    for(l=0;l<MFAtlasK(A,e);l++)(MFKV_CStar(vf,e))[l]=0;
    MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
    for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterpee," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);IMFNInterpee++;
    for(l=0;l<MFAtlasK(A,e);l++)(MFKV_CStar(vf,e))[l]=0;
    (MFKV_CStar(vf,e))[0]=1.2*R;
    MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
    for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterpee," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);IMFNInterpee++;

    for(l=0;l<MFAtlasK(A,e);l++)(MFKV_CStar(vf,e))[l]=0;
    MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
    for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterpee," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);IMFNInterpee++;
    for(l=0;l<MFAtlasK(A,e);l++)(MFKV_CStar(vf,e))[l]=0;
    (MFKV_CStar(vf,e))[1]=-1.2*R;
    MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
    for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterpee," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);IMFNInterpee++;

    if(0)
     { /* Draw Polygon */
      for(i=0;i<MFPolytopeNumberOfVertices(P,e);i++)
       {
        n1=MFPolytopeVertexNumberOfIndices(P,i,e);
        for(j=i+1;j<MFPolytopeNumberOfVertices(P,e);j++)
         {
          n2=MFPolytopeVertexNumberOfIndices(P,j,e);
          m=n1;if(n2>m)m=n2;

          if(m>0)
           {
            n=MFPolytopeIntersectIndexSets(P,i,j,inter,e);
           }else{
            n=0;
           }

          if(n==k-1)
           {
            MFPolytopeVertex(P,i,vf,e);
            MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
            for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterpee," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);
            IMFNInterpee++;
            MFPolytopeVertex(P,j,vf,e);
            MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
            for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterpee," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);
            IMFNInterpee++;
           }
         }
       }
     }
    if(1)
     { /* Draw Ball */
      for(i=0;i<301;i++)
       {
        f[0]=R*cos(2*3.1415926*i/300.);
        f[1]=R*sin(2*3.1415926*i/300.);
        good0=1;
        for(ii=0;ii<MFPolytopeNumberOfFaces(P,e);ii++)
         {
          oo=MFPolytopeFaceOrigin(P,ii,e);
          nn=MFKV_CStar(MFPolytopeFaceNormal(P,ii,e),e);
          d0=-oo;for(l=0;l<k;l++)d0+=f[l]*nn[l];
          if(d0>0)good0=0;
         }

        MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);

        f[0]=R*cos(2*3.1415926*(i+1)/300.);
        f[1]=R*sin(2*3.1415926*(i+1)/300.);
        for(ii=0;ii<MFPolytopeNumberOfFaces(P,e);ii++)
         {
          oo=MFPolytopeFaceOrigin(P,ii,e);
          nn=MFKV_CStar(MFPolytopeFaceNormal(P,ii,e),e);
          d0=-oo;for(l=0;l<k;l++)d0+=f[l]*nn[l];
          if(d0>0)good0=0;
         }

        if(good0)
         {
          for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFCircle," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFCircle,"\n");fflush(IMFCircle);IMFNCircle++;
          MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
          for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFCircle," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFCircle,"\n");fflush(IMFCircle);
          IMFNCircle++;
         }
       }
     }
   }

  if(0)
   {
    for(iBndChart=0;iBndChart<nBndCharts;iBndChart++)
     {
      chart=MFAtlasChartWithBoundary(A,iBndChart,e);
      P=MFChartPolytope(MFAtlasChart(A,chart,e),e);

      R=MFChartRadius(MFAtlasChart(A,chart,e),e);
      for(i=0;i<MFPolytopeNumberOfVertices(P,e);i++)
       {
        n1=MFPolytopeVertexNumberOfIndices(P,i,e);
        for(j=i+1;j<MFPolytopeNumberOfVertices(P,e);j++)
         {
          n2=MFPolytopeVertexNumberOfIndices(P,j,e);
          m=n1;if(n2>m)m=n2;

          if(m>0)
           {
            n=MFPolytopeIntersectIndexSets(P,i,j,inter,e);
           }else{
            n=0;
           }

          noncube=1;for(l=0;l<n;l++)noncube=noncube&&inter[l]>m;

          if(n==k-1&&noncube)
           {
            for(ii=0;ii<k-1;ii++)
             {
              iface[ii]=-1;
              for(l=0;l<MFPolytopeNumberOfFaces(P,e);l++)
                if(MFPolytopeFaceIndex(P,l,e)==inter[ii])iface[ii]=l;
             }
  
            for(ii=0;ii<k-1;ii++)
             {
              b[ii]=MFPolytopeFaceOrigin(P,iface[ii],e);
              nn=MFKV_CStar(MFPolytopeFaceNormal(P,iface[ii],e),e);
              d0;for(l=0;l<k;l++)d0+=nn[l]*nn[l];
              for(l=0;l<k;l++)
               {
                a[ii+(k-1)*l]=nn[l]/d0;
                asave[ii+(k-1)*l]=a[ii+(k-1)*l];
               }
              b[ii]=b[ii]/d0;
             }

            for(ii=0;ii<k-1;ii++)
             {
              for(jj=0;jj<k-1;jj++)
               {
                aa[ii+(k-1)*jj]=0;
                for(l=0;l<k;l++)aa[ii+(k-1)*jj]+=a[ii+(k-1)*l]*a[jj+(k-1)*l];
               }
             }

            MFNullVectorSmall(k,a,nrm,e);
            MFSolveFullSmall(k-1,aa,b,e);

            for(ii=0;ii<k;ii++)
             {
              org[ii]=0.;
              for(l=0;l<k-1;l++)org[ii]+=asave[l+(k-1)*ii]*b[l];
             }

            r=R*R;for(l=0;l<k;l++)r-=org[l]*org[l];
            if(r>0)
             {
              r=sqrt(r);

              good0=1;
              good1=1;
              for(ii=0;ii<MFPolytopeNumberOfFaces(P,e);ii++)
               {
                doit=1;for(l=0;l<k-1;l++)if(ii==iface[l])doit=0;
                if(doit)
                 {
                  oo=MFPolytopeFaceOrigin(P,ii,e);
                  nn=MFKV_CStar(MFPolytopeFaceNormal(P,ii,e),e);
                  d0=-oo;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l])*nn[l];
                  if(d0>0)good0=0;
                  d1=-oo;for(l=0;l<k;l++)d1+=(org[l]-r*nrm[l])*nn[l];
                  if(d1>0)good1=0;
                 }
               }

              if(good0)
               {
                for(l=0;l<k;l++)f[l]=org[l]+r*nrm[l];
                MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
                for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterper," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterper,"\n");fflush(IMFInterper);
                IMFNInterper++;
                if(1||inter[0]<inter[1])
                 {
/*                for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFCircle," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFCircle,"\n");fflush(IMFCircle);*/
/*                IMFNCircle++;*/
                  IMFEvaluateFlow(F,vp,p0,pF,e);
                  MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);
                  d0=0.;for(l=0;l<k;l++)d0+=f[l]*f[l];
                  for(l=0;l<k;l++)f[l]=org[l]+r*nrm[l]+.2*f[l]/sqrt(d0);
                  MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
/*                for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFCircle," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFCircle,"\n");fflush(IMFCircle);*/
/*                IMFNCircle++;*/
                 }
                for(l=0;l<k;l++)f[l]=org[l]+r*nrm[l];
                MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);

                IMFEvaluateFlow(F,vp,p0,pF,e);
                MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);
                d1=0.;for(l=0;l<k;l++)d1+=f[l]*nrm[l];
                t=-r/d1;
                if(t<0.)
                 {
                  good0=1;
                  for(ii=0;ii<k-1;ii++)
                   {
                    d0=0.;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l]+t*f[l])*asave[ii+(k-1)*l];
                    d1=0.;for(l=0;l<k;l++)d1+= org[l]                 *asave[ii+(k-1)*l];
                    if(d0<0 || d0>d1)good0=0;
                   }
                  if(good0)
                   {
                    for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterper," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterper,"\n");fflush(IMFInterper);
                    IMFNInterper++;
                   }
                 }
               }

              if(good1)
               {
                for(l=0;l<k;l++)f[l]=org[l]-r*nrm[l];
                MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
                for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterper," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterper,"\n");fflush(IMFInterper);
                IMFNInterper++;
                if(1||inter[0]<inter[1])
                 {
/*                for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFCircle," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFCircle,"\n");fflush(IMFCircle);
                  IMFNCircle++;*/
                  IMFEvaluateFlow(F,vp,p0,pF,e);
                  MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);
                  d0=0.;for(l=0;l<k;l++)d0+=f[l]*f[l];
                  for(l=0;l<k;l++)f[l]=org[l]-r*nrm[l]+.2*f[l]/sqrt(d0);
                  MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
/*                for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFCircle," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFCircle,"\n");fflush(IMFCircle);
                  IMFNCircle++;*/
                 }
                for(l=0;l<k;l++)f[l]=org[l]-r*nrm[l];
                MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);

                IMFEvaluateFlow(F,vp,p0,pF,e);
                MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);
                d1=0.;for(l=0;l<k;l++)d1+=f[l]*nrm[l];
                t=r/d1;
                if(t<0.)
                 {
                  good1=1;
                  for(ii=0;ii<k-1;ii++)
                   {
                    d0=0.;for(l=0;l<k;l++)d0+=(org[l]-r*nrm[l]+t*f[l])*asave[ii+(k-1)*l];
                    d1=0.;for(l=0;l<k;l++)d1+= org[l]                 *asave[ii+(k-1)*l];
                    if(d0<0 || d0>d1)good1=0;
                   }
                  if(good1)
                   {
                    for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterper," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterper,"\n");fflush(IMFInterper);
                    IMFNInterper++;
                   }
                 }
               }
             }
           }
         }
       }
     }
   }
    printf("*** done drawing balls\n");fflush(stdout);
#endif

  for(iBndChart=0;iBndChart<nBndCharts;iBndChart++)
   {
    chart=MFAtlasChartWithBoundary(A,iBndChart,e);
    P=MFChartPolytope(MFAtlasChart(A,chart,e),e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("The %dth chart on the boundary list is chart %d\n",iBndChart,chart);fflush(stdout);}
#endif

    nv=MFPolytopeNumberOfVertices(P,e);
/*  if(!MFChartChangedFlag(MFAtlasChart(A,chart,e),e))nv=0;*/

    R=MFChartRadius(MFAtlasChart(A,chart,e),e);
    for(i=0;i<nv;i++)
     {
      n1=MFPolytopeVertexNumberOfIndices(P,i,e);
      for(j=i+1;j<nv;j++)
       {
        n2=MFPolytopeVertexNumberOfIndices(P,j,e);
        m=n1;if(n2>m)m=n2;

        if(m>0)
         {
          n=MFPolytopeIntersectIndexSets(P,i,j,inter,e);
         }else{
          n=0;
         }

        noncube=1;for(l=0;l<n;l++)noncube=noncube&&inter[l]>m;

        if(n==k-1&&noncube)
         {

#ifdef MFALLOWVERBOSE
          if(verbose){printf("  Indices of this non-cube edge:\n");fflush(stdout);}
#endif

          for(ii=0;ii<k-1;ii++)
           {
            iface[ii]=-1;
            for(l=0;l<MFPolytopeNumberOfFaces(P,e);l++)
             {
              if(MFPolytopeFaceIndex(P,l,e)==inter[ii])iface[ii]=l;
             }

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  %d (face %d)\n",inter[ii],iface[ii]);fflush(stdout);}
#endif

           }

          for(ii=0;ii<k-1;ii++)
           {
            b[ii]=MFPolytopeFaceOrigin(P,iface[ii],e);
            nn=MFKV_CStar(MFPolytopeFaceNormal(P,iface[ii],e),e);
            d0;for(l=0;l<k;l++)d0+=nn[l]*nn[l];
            for(l=0;l<k;l++)
             {
              a[ii+(k-1)*l]=nn[l]/d0;
              asave[ii+(k-1)*l]=a[ii+(k-1)*l];
             }
            b[ii]=b[ii]/d0;
           }

          for(ii=0;ii<k-1;ii++)
           {
            for(jj=0;jj<k-1;jj++)
             {
              aa[ii+(k-1)*jj]=0;
              for(l=0;l<k;l++)aa[ii+(k-1)*jj]+=a[ii+(k-1)*l]*a[jj+(k-1)*l];
             }
           }

          MFNullVectorSmall(k,a,nrm,e);
          MFSolveFullSmall(k-1,aa,b,e);

          for(ii=0;ii<k;ii++)
           {
            org[ii]=0.;
            for(l=0;l<k-1;l++)org[ii]+=asave[l+(k-1)*ii]*b[l];
           }

#ifdef MFALLOWVERBOSE
          if(verbose)
           {
            printf("  Done finding edge equation org=(%lf",org[0]);fflush(stdout);
            for(ii=1;ii<k;ii++)printf(",%lf",org[ii]);
            printf(")\n");fflush(stdout);
            printf("                             nrm=(%lf",nrm[0]);fflush(stdout);
            for(ii=1;ii<k;ii++)printf(",%lf",nrm[ii]);
            printf(")\n");fflush(stdout);
    
            printf("Check: is org on all the faces?\n");fflush(stdout);
            d0=0.;for(l=0;l<k;l++)d0+=org[l]*nrm[l];
            printf("            nrm.org=%lf\n",d0);fflush(stdout);
            for(ii=0;ii<k-1;ii++)
             {
              oo=MFPolytopeFaceOrigin(P,iface[ii],e);
              nn=MFKV_CStar(MFPolytopeFaceNormal(P,iface[ii],e),e);
              d0=-oo;for(l=0;l<k;l++)d0+=nn[l]*org[l];
              printf("   face %d: org.n-o=%lf\n",ii,d0);fflush(stdout);
              d0=0.;for(l=0;l<k;l++)d0+=nn[l]*nrm[l];
              printf("            nrm.n=%lf\n",d0);fflush(stdout);
             }
           }
#endif

          r=R*R;for(l=0;l<k;l++)r-=org[l]*org[l];
          if(r>0)
           {
            r=sqrt(r);

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              d0=0.;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l])*(org[l]+r*nrm[l]);d0=sqrt(d0);
              printf("  This edge crosses the ball. |pt0|=%lf, R=%lf\n",d0,R);fflush(stdout);
              d0=0.;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l])*(org[l]+r*nrm[l]);d0=sqrt(d0);
              printf("                              |pt1|=%lf, R=%lf\n",d0,R);fflush(stdout);
             }
#endif

            good0=1;
            good1=1;
            for(ii=0;ii<MFPolytopeNumberOfFaces(P,e);ii++)
             {
              doit=1;for(l=0;l<k-1;l++)if(ii==iface[l])doit=0;
              if(doit)
               {
                oo=MFPolytopeFaceOrigin(P,ii,e);
                nn=MFKV_CStar(MFPolytopeFaceNormal(P,ii,e),e);
                d0=-oo;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l])*nn[l];
                if(d0>0)good0=0;
                d1=-oo;for(l=0;l<k;l++)d1+=(org[l]-r*nrm[l])*nn[l];
                if(d1>0)good1=0;
               }
             }

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              if(good0){printf("  The +r point is in the polytope.\n");fflush(stdout);}
                  else {printf("  The +r point is not in the polytope.\n");fflush(stdout);}
             }
#endif

            if(good0)
             {
              for(l=0;l<k;l++)f[l]=org[l]+r*nrm[l];
              MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
              IMFEvaluateFlow(F,vp,p0,pF,e);
              MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);
              d1=0.;for(l=0;l<k;l++)d1+=f[l]*nrm[l];
              t=-r/d1;
              if(t<0.&&(Omega==NULL||MFNRegionInterior(Omega,vp,e),e))
               {
                good0=1;
                for(ii=0;ii<k-1;ii++)
                 {
                  d0=0.;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l]+t*f[l])*asave[ii+(k-1)*l];
                  d1=0.;for(l=0;l<k;l++)d1+= org[l]                 *asave[ii+(k-1)*l];
                  if(d0<0 || d0>d1)good0=0;
                 }
                if(1||good0)
                 {
                  for(l=0;l<k;l++)f[l]=org[l]+r*nrm[l]+t*f[l];
                  MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,z,e);
                  result=IMFInterpolateOnDualFace(A,F,p0,z,chart,inter,e);
/*                result=IMFInterpolateOnDualFace2(A,F,p0,c,vf,z,chart,inter,e);*/
                  if(result!=NULL)
                   {
                    if(IMFExpansionNVGetT(result,e)<tmax&&(Omega==NULL||MFNRegionInterior(Omega,IMFExpansionNVGetSigma(result,e),e)))
                     {
#ifdef DOINTERPANIM
                      if(0){
                      for(j=0;j<MFNV_NC(vp,e);j++)fprintf(IMFInterpee," %lf",(MFNV_CStar(vp,e))[j]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);
                      IMFNInterpee++;
                      for(j=0;j<MFNV_NC(result,e);j++)fprintf(IMFInterpee," %lf",IMFExpansionU(IMFExpansionNVGetE(result,e),e)[j]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);
                      IMFNInterpee++;}
#endif

#ifdef MFALLOWVERBOSE
                      if(verbose){printf("Result is on chart %d\n",chart);fflush(stdout);}
#endif
                      MFChartResetChangedFlag(MFAtlasChart(A,chart,e),e);
                      goto FreeAndReturn;
                     }
                   }
                 }
               }
             }

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              if(good1){printf("  The -r point is in the polytope.\n");fflush(stdout);}
                  else {printf("  The +r point is not in the polytope.\n");fflush(stdout);}
             }
#endif

            if(good1)
             {
              for(l=0;l<k;l++)f[l]=org[l]-r*nrm[l];
              MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
              IMFEvaluateFlow(F,vp,p0,pF,e);
              MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);
              d1=0.;for(l=0;l<k;l++)d1+=f[l]*nrm[l];
              t=r/d1;
              if(t<0.&&(Omega==NULL||MFNRegionInterior(Omega,vp,e),e))
               {
                good1=1;
                for(ii=0;ii<k-1;ii++)
                 {
                  d0=0.;for(l=0;l<k;l++)d0+=(org[l]-r*nrm[l]+t*f[l])*asave[ii+(k-1)*l];
                  d1=0.;for(l=0;l<k;l++)d1+= org[l]                 *asave[ii+(k-1)*l];
                  if(d0<0 || d0>d1)good1=0;
                 }
                if(1||good1)
                 {
                  for(l=0;l<k;l++)f[l]=org[l]-r*nrm[l]+t*f[l];
                  MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,z,e);
                  result=IMFInterpolateOnDualFace(A,F,p0,z,chart,inter,e);
/*                result=IMFInterpolateOnDualFace2(A,F,A,p0,vf,z,chart,inter,e);*/
                  if(result!=NULL)
                   {
                    if(IMFExpansionNVGetT(result,e)<tmax&&(Omega==NULL||MFNRegionInterior(Omega,IMFExpansionNVGetSigma(result,e),e)))
                     {
#ifdef DOINTERPANIM
                    if(0){
                      for(j=0;j<MFNV_NC(vp,e);j++)fprintf(IMFInterpee," %lf",(MFNV_CStar(vp,e))[j]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);
                      IMFNInterpee++;
                      for(j=0;j<MFNV_NC(result,e);j++)fprintf(IMFInterpee," %lf",IMFExpansionU((IMFExpansionNVGetE(result,e)),e)[j],e);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);
                      IMFNInterpee++;}
#endif

#ifdef MFALLOWVERBOSE
                      if(verbose){printf("Result is on chart %d\n",chart);fflush(stdout);}
#endif

                      MFChartResetChangedFlag(MFAtlasChart(A,chart,e),e);
                      goto FreeAndReturn;
                     }
                   }
                 }
               }
             }
           }
         }
       }
     }
    MFChartResetChangedFlag(MFAtlasChart(A,chart,e),e);
   }

  result=NULL;
  printf("No such points on the boundary!!!\n");

FreeAndReturn:
  MFFreeKVector(vf,e);
  MFFreeNVector(z,e);
  MFFreeNVector(vp,e);
  MFFreeNVector(vF,e);

  return result;
 }

void MFNullVectorSmall(int k,double *A,double *phi, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNullVectorSmall"};

#ifdef HAVE_LAPACK

/*Find the null vector of a (k-1)xk matrix */

  double t0,t1;

  if(k==2)
   {
    t1=sqrt(A[0]*A[0]+A[1]*A[1]);
    t0    = A[1]/t1;
    phi[1]=-A[0]/t1;
    phi[0]=t0;
   }else if(k==3)
   {
    t0=sqrt(A[0]*A[0]+A[2]*A[2]+A[4]*A[4]);
    t1=sqrt(A[1]*A[1]+A[3]*A[3]+A[5]*A[5]);
    phi[0]=(A[2]*A[5]-A[3]*A[4])/t0/t1;
    phi[1]=(A[4]*A[1]-A[5]*A[0])/t0/t1;
    phi[2]=(A[0]*A[3]-A[1]*A[2])/t0/t1;
   }else{  
    int i,j;
    char jobvl='N';
    char jobvr='A';
    int lda;
    double *s;
    double *U;
    int ldU;
    double *V;
    int ldV;
    double *work;
    int lwork;
    int info=0;
    int m;
    double t;

    lda=k-1;
    ldU=k-1;
    U=NULL;

    ldV=k;
    V=(double*)malloc(k*(k-1)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(V==NULL)
   {
    sprintf(IMFInterpolationErrorMsg,"Out of memory, trying to allocate %d bytes",k*(k-1)*sizeof(double));
    MFSetError(e,12,RoutineName,IMFInterpolationErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

    s=(double*)malloc(k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(s==NULL)
   {
    sprintf(IMFInterpolationErrorMsg,"Out of memory, trying to allocate %d bytes",k*sizeof(double));
    MFSetError(e,12,RoutineName,IMFInterpolationErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

    i=-1;
    info=0;
    jobvl='N';
    jobvr='A';
    CALLDGESVD(&jobvl,&jobvr,&lda,&k,A,&lda,s,U,&ldU,V,&ldV,&t,&i,&info);
    lwork=round(t);

    work=(double*)malloc(lwork*sizeof(double));

#ifndef MFNOSAFETYNET
  if(work==NULL)
   {
    sprintf(IMFInterpolationErrorMsg,"Out of memory, trying to allocate %d bytes",lwork*sizeof(double));
    MFSetError(e,12,RoutineName,IMFInterpolationErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

    info=0;
    jobvl='N';
    jobvr='A';
    CALLDGESVD(&jobvl,&jobvr,&lda,&k,A,&lda,s,U,&ldU,V,&ldV,work,&lwork,&info);
  
    m=0;
    for(i=0;i<k;i++)
     {
      if(fabs(s[i])<1.e-7)
       {
        for(j=0;j<k;j++)phi[j]=V[i+ldV*j];
        m++;
       }
     }
    free(V);
    free(s);
    free(work);
   }
  return;
#else
  sprintf(IMFInterpolationErrorMsg,"IMFInterpolation requires dgesvd from Lapack");
  MFSetError(e,12,RoutineName,IMFInterpolationErrorMsg,__LINE__,__FILE__);
  return;
#endif
 }

void MFSolveFullSmall(int k,double *A,double *b, MFErrorHandler e)
 {
  static char RoutineName[]={"MFSolveFullSmall"};
/*solves the kxk system of linear eqs.;*/
  double t;

  if(k==1)
   {
    b[0]=b[0]/A[0];
   }else if(k==2)
   {
    t   =(A[3]*b[0]-A[2]*b[1])/(A[0]*A[3]-A[1]*A[2]);
    b[1]=(A[0]*b[1]-A[1]*b[0])/(A[0]*A[3]-A[1]*A[2]);
    b[0]=t;
   }else{  
    MFSolveFull(k,A,b,e);
   }

  return;
 }

MFNVector IMFInterpolateOnDualFace(MFAtlas A, IMFFlow F, MFKVector p0, MFNVector z, int chart, int *inter, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFInterpolateOnDualFace"};
  int i,k,n;
  int j;
  MFNVector c[kmax];
  double a[(kmax-1)*(kmax-1)];
  double t[kmax];
  double d;
  MFNSpace space;
  MFNVector diffi;
  MFNVector diffj;
  MFNVector diffz;
  MFNVector result;
  MFNVector sigma;
  int nsigma;
  double time;
  IMFExpansion U,uu;
  IMFExpansion V0,V1;
  MFKVector vs,s0;
  double *s;

  double *u;
  double *du;
  double *ds;
  double *dds;
  int good;
  int verbose=0;

  k=MFAtlasK(A,e);
  n=MFAtlasN(A,e);
  space=MFIMFNSpace(MFAtlasMF(A,e),e);

  c[0]=MFAtlasChartCenter(A,chart,e);
  diffi=MFCreateNVector(n,e);
  diffj=MFCreateNVector(n,e);
  diffz=MFCreateNVector(n,e);
  vs=MFCreateKVector(k,e);
  s=MFKV_CStar(vs,e);
  s0=MFCreateKVector(k,e);
  
  if(0){printf("IMFInterpolateOnDualFace, z=");MFPrintNVector(stdout,z,e);printf("\n");fflush(stdout);}
  nsigma=MFNV_NC(IMFExpansionNVGetSigma(c[0],e),e);

  good=0;
  for(i=0;i<k-1;i++)
   {
    j=MFAtlasHalfSpaceRightChart(A,inter[i],e);
    if(j==chart)j=MFAtlasHalfSpaceLeftChart(A,inter[i],e);
    c[i+1]=MFAtlasChartCenter(A,j,e);
    d=0.;for(j=0;j<nsigma;j++)d+=pow(MFNV_CStar(IMFExpansionNVGetSigma(c[0],e),e)[j]-
                                     MFNV_CStar(IMFExpansionNVGetSigma(c[i+1],e),e)[j],2.);
    if(sqrt(d)>1.e-14)good=1;
   }

  if(!good)
   {
    MFFreeNVector(diffi,e);
    MFFreeNVector(diffj,e);
    MFFreeNVector(diffz,e);
    MFFreeKVector(vs,e);
    MFFreeKVector(s0,e);

    return NULL;
   }

  MFNSpaceDirection(space,z,c[0],diffz,e);
  for(i=0;i<k-1;i++)
   {
    MFNSpaceDirection(space,c[i+1],c[0],diffi,e);
    a[i+(k-1)*i]=MFNSpaceInner(space,diffi,diffi,e);
    for(j=i+1;j<k-1;j++)
     {
      MFNSpaceDirection(space,c[j+1],c[0],diffj,e);
      a[i+(k-1)*j]=MFNSpaceInner(space,diffi,diffj,e);
      a[i+(k-1)*j]=a[j+(k-1)*i];
     }
    t[i+1]=MFNSpaceInner(space,diffz,diffi,e);
   }

  MFSolveFullSmall(k-1,a,t+1,e);

  t[0]=1.;for(i=1;i<k;i++)t[0]-=t[i];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Barycentric coordinates (%lf",t[0]);
    for(i=1;i<k;i++)printf(",%lf",t[i]);
    printf(")\n");
   }
#endif

  good=1;
  for(i=0;i<k;i++)if(t[i]<0.||t[i]>1.)good=0;
  for(i=0;i<k;i++)if(t[i]<.1||t[i]>.9)good=0;  /* EXPERIMENT */
  if(!good)
   {
    MFFreeNVector(diffi,e);
    MFFreeNVector(diffj,e);
    MFFreeNVector(diffz,e);
    MFFreeKVector(vs,e);
    MFFreeKVector(s0,e);
    return NULL;
   }

  if(k!=2)
   {
    fprintf(stderr,"Interpolation of more than 2d is not yet supported\n");
    fflush(stderr);
    abort();
   }

  u=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(u==NULL)
   {
    sprintf(IMFInterpolationErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFInterpolationErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  du=(double*)malloc(n*k*sizeof(double));

#ifndef MFNOSAFETYNET
  if(du==NULL)
   {
    sprintf(IMFInterpolationErrorMsg,"Out of memory, trying to allocate %d bytes",n*k*sizeof(double));
    MFSetError(e,12,RoutineName,IMFInterpolationErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  ds=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(ds==NULL)
   {
    sprintf(IMFInterpolationErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFInterpolationErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  dds=(double*)malloc(n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dds==NULL)
   {
    sprintf(IMFInterpolationErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFInterpolationErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif


/*
  position is straightforward, as are t and sigma.
  tangent and second derivatives get as scalars along direction from c[0] to c[1]
    then interpolate u_s and u_ss 
  then use the flow to build the full expansion. (same as initial curve)
*/

  sigma=MFCreateNVector(nsigma,e);
  V0=IMFExpansionNVGetE(c[0],e);
  V1=IMFExpansionNVGetE(c[1],e);
  for(j=0;j<n;j++)u[j]=t[0]*(IMFExpansionU(V0,e))[j]+t[1]*(IMFExpansionU(V1,e))[j];
  time=t[0]*IMFExpansionNVGetT(c[0],e)+t[1]*IMFExpansionNVGetT(c[1],e);
  for(j=0;j<nsigma;j++)MFNV_CStar(sigma,e)[j]=t[0]*MFNV_CStar(IMFExpansionNVGetSigma(c[0],e),e)[j];
                                           +t[1]*MFNV_CStar(IMFExpansionNVGetSigma(c[1],e),e)[j];

#ifdef DOINTERPANIM
  if(0){
  for(j=0;j<n;j++)fprintf(IMFInterpee," %lf",(IMFExpansionU(V0,e))[j]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);IMFNInterpee++;
  for(j=0;j<n;j++)fprintf(IMFInterpee," %lf",(IMFExpansionU(V1,e))[j]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);IMFNInterpee++;
  }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Barycentric coordinates (%lf",t[0]);
    for(i=1;i<k;i++)printf(",%lf",t[i]);
    printf(")\n");
   }
#endif

  MFKV_CStar(s0,e)[0]=0.;
  MFKV_CStar(s0,e)[1]=0.;
  MFNSpaceDirection(space,c[0],c[1],diffi,e);
  d=0.;for(j=0;j<n;j++)d+=MFNV_CStar(diffi,e)[j]*MFNV_CStar(diffi,e)[j];d=sqrt(d);

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Vector from 0 to 1 (%lf",MFNV_CStar(diffi,e)[0])/d;
    for(i=1;i<n;i++)printf(",%lf",MFNV_CStar(diffi,e)[i]/d);
    printf(")\n");
   }
#endif

  MFMVMulT(MFIMFNSpace(MFAtlasMF(A,e),e),MFAtlasChartTangentSpace(A,chart,e),diffi,vs,e);
  d=0.;for(j=0;j<k;j++)d+=s[j]*s[j];d=sqrt(d);
  for(j=0;j<k;j++)s[j]=s[j]/d;
  IMFEvaluateExpansionDirectionalDerivative(V0,s0,vs,diffj,e);
  for(j=0;j<n;j++)ds[j]=t[0]*MFNV_CStar(diffj,e)[j];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Tangent 0          (%lf",MFNV_CStar(diffj,e)[0]);
    for(i=1;i<n;i++)printf(",%lf",MFNV_CStar(diffj,e)[i]);
    printf(")\n");
   }
#endif

  IMFEvaluateExpansionSecondDirectionalDerivative(V0,s0,vs,vs,diffj,e);
  for(j=0;j<n;j++)dds[j]=t[0]*MFNV_CStar(diffj,e)[j];

  j=MFAtlasHalfSpaceRightChart(A,inter[0],e);
  if(j==chart)j=MFAtlasHalfSpaceLeftChart(A,inter[0],e);
  MFMVMulT(MFIMFNSpace(MFAtlasMF(A,e),e),MFAtlasChartTangentSpace(A,j,e),diffi,vs,e);
  d=0.;for(j=0;j<k;j++)d+=s[j]*s[j];d=sqrt(d);
  for(j=0;j<k;j++)s[j]=s[j]/d;
  IMFEvaluateExpansionDirectionalDerivative(V1,s0,vs,diffj,e);
  for(j=0;j<n;j++)ds[j]+=t[1]*MFNV_CStar(diffj,e)[j];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Tangent 1          (%lf",MFNV_CStar(diffj,e)[0]);
    for(i=1;i<n;i++)printf(",%lf",MFNV_CStar(diffj,e)[i]);
    printf(")\n");
    printf("Interpolated       (%lf",ds[0]);
    for(i=1;i<n;i++)printf(",%lf",ds[i]);
    printf(")\n");
   }
#endif

  IMFEvaluateExpansionSecondDirectionalDerivative(V1,s0,vs,vs,diffj,e);
  for(j=0;j<n;j++)dds[j]+=t[1]*MFNV_CStar(diffj,e)[j];

  uu=IMFCreateExpansion(n,1,e);
  IMFExpansionSetDerivatives(uu,u,ds,dds,NULL,e);
  U=IMFInflateExpansionWithFlow(uu,F,p0,e);

/*printf("Interpolated u, (%lf",u[0]);for(j=1;j<n;j++)printf(",%lf",u[j]);printf(")");fflush(stdout);
  printf(", time %lf",time);
  printf(", sigma");MFPrintNVector(stdout,sigma,e);printf("\n");fflush(stdout);*/

  result=IMFCreateExpansionNVector(U,time,sigma,chart,2,e);

  IMFExpansionNVSetChart0(result,chart,e);
  IMFExpansionNVSetS0(result,s0,e);

  MFFreeNVector(diffi,e);
  MFFreeNVector(diffj,e);
  MFFreeNVector(diffz,e);
  IMFFreeExpansion(U,e);
  IMFFreeExpansion(uu,e);
  MFFreeNVector(sigma,e);
  MFFreeKVector(vs,e);
  MFFreeKVector(s0,e);

  free(u);
  free(du);
  free(ds);
  free(dds);

  return result;
 }

MFNVector IMFGetInterpolationPointOnList(MFAtlas A,IMFFlow F, MFKVector p0, MFAtlas c, double tmax, MFNRegion Omega, int nList, int *list, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFGetInterpolationPointOnList"};
/*
    edges whose index is all >m give points where edge crosses sphere
    clip points against all faces with index >m not in index set of edge
    mark intersection point
    index of edge is ind1,ind2
*/
  int verbose=0;

  MFPolytope P;
  int iBndChart,nBndCharts;
  int chart;
  int doit;

  int i,j,k,l;
  int n,m;
  int n1,n2;
  double d0,d1,t;
  double r;
  int good0;
  int good1;
  int inter[kmax];
  MFNVector result;
  double R;

  int ii,jj;
  int noncube;
  double a[(kmax-1)*kmax],b[kmax-1],aa[(kmax-1)*(kmax-1)],asave[(kmax-1)*kmax];
  int iface[kmax-1];
  double d[kmax];
  double org[kmax],nrm[kmax];
  double oo;
  double *nn;
  double *f;
  MFKVector vf;
  MFNVector z;
  MFNVector vp;
  double *p;
  MFNVector vF;
  double *pF;
  int nv;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

#ifdef DOINTERPANIM
  if(IMFInterper!=NULL)
   {
    fclose(IMFInterper);
    IMFInterper=NULL;
   }
  if(IMFInterper==NULL)
   {
    IMFInterper=fopen("Interp1.dx","w");
    IMFNInterper=0;
   }
  if(IMFInterpee==NULL)
   {
    IMFInterpee=fopen("Interp2.dx","w");
    IMFNInterpee=0;
   }
  if(IMFCircle!=NULL)
   {
    fclose(IMFCircle);
    IMFCircle=NULL;
   }
  if(IMFCircle==NULL)
   {
    IMFCircle=fopen("IMFCircle.dx","w");
    IMFNCircle=0;
   }
#endif

  nBndCharts=MFAtlasNumberOfChartsWithBoundary(A,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("There are %d charts on the boundary list\n",nBndCharts);fflush(stdout);}
#endif

  k=MFAtlasK(A,e);
  if(k>kmax)
   {
    fprintf(stderr,"%s, I statically allocated some things here for k<10! Complain.\n");
    fflush(stderr);
    abort();
   }

  m=1;
  for(i=0;i<k;i++)m=m*2;
  m=m-1;

  vf=MFCreateKVector(k,e);
  f=MFKV_CStar(vf,e);
  z=MFCreateNVector(MFAtlasN(A,e),e);
  vp=MFCreateNVector(MFAtlasN(A,e),e);
  p=MFNV_CStar(vp,e);
  vF=MFCreateNVector(MFAtlasN(A,e),e);
  pF=MFNV_CStar(vF,e);

/* Clean boundary list */

  for(iBndChart=0;iBndChart<nBndCharts;iBndChart++)
   {
    chart=MFAtlasChartWithBoundary(A,iBndChart,e);
    P=MFChartPolytope(MFAtlasChart(A,chart,e),e);

    d0=0.;
    for(i=0;i<MFPolytopeNumberOfVertices(P,e);i++)
     {
      d1=MFPolytopeRadiusOfVertex(MFChartPolytope(MFAtlasChart(A,chart,e),e),i,e);
      if(d1>d0)d0=d1;
     }

    R=MFChartRadius(MFAtlasChart(A,chart,e),e);
    if(MFChartIsSingular(MFAtlasChart(A,chart,e),e)
           ||d0<R
           ||IMFExpansionNVGetT(MFChartCenter(MFAtlasChart(A,chart,e),e),e)>tmax
           ||(Omega!=NULL&&!MFNRegionInterior(Omega,IMFExpansionNVGetSigma(MFChartCenter(MFAtlasChart(A,chart,e),e),e),e)) )
     {
      MFAtlasRemoveChartFromBoundaryList(A,iBndChart,e);
      iBndChart--;
      nBndCharts--;
     }
   }

#ifdef DOINTERPANIM
  nBndCharts=MFAtlasNumberOfChartsWithBoundary(A,e);
  for(iBndChart=0;iBndChart<nBndCharts;iBndChart++)
   {
    chart=MFAtlasChartWithBoundary(A,iBndChart,e);
    P=MFChartPolytope(MFAtlasChart(A,chart,e),e);
    R=MFChartRadius(MFAtlasChart(A,chart,e),e);

    if(0)
     { /* Draw Polygon */
      for(i=0;i<MFPolytopeNumberOfVertices(P,e);i++)
       {
        n1=MFPolytopeVertexNumberOfIndices(P,i,e);
        for(j=i+1;j<MFPolytopeNumberOfVertices(P,e);j++)
         {
          n2=MFPolytopeVertexNumberOfIndices(P,j,e);
          m=n1;if(n2>m)m=n2;

          if(m>0)
           {
            n=MFPolytopeIntersectIndexSets(P,i,j,inter,e);
           }else{
            n=0;
           }

          if(n==k-1)
           {
            MFPolytopeVertex(P,i,vf,e);
            MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
            for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterpee," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);
            IMFNInterpee++;
            MFPolytopeVertex(P,j,vf,e);
            MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
            for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterpee," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);
            IMFNInterpee++;
           }
         }
       }
     }
    if(1)
     { /* Draw Ball */
      for(i=0;i<301;i++)
       {
        f[0]=R*cos(2*3.1415926*i/300.);
        f[1]=R*sin(2*3.1415926*i/300.);
        good0=1;
        for(ii=0;ii<MFPolytopeNumberOfFaces(P,e);ii++)
         {
          oo=MFPolytopeFaceOrigin(P,ii,e);
          nn=MFKV_CStar(MFPolytopeFaceNormal(P,ii,e),e);
          d0=-oo;for(l=0;l<k;l++)d0+=f[l]*nn[l];
          if(d0>0)good0=0;
         }

        MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);

        f[0]=R*cos(2*3.1415926*(i+1)/300.);
        f[1]=R*sin(2*3.1415926*(i+1)/300.);
        for(ii=0;ii<MFPolytopeNumberOfFaces(P,e);ii++)
         {
          oo=MFPolytopeFaceOrigin(P,ii,e);
          nn=MFKV_CStar(MFPolytopeFaceNormal(P,ii,e),e);
          d0=-oo;for(l=0;l<k;l++)d0+=f[l]*nn[l];
          if(d0>0)good0=0;
         }

        if(good0)
         {
          for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFCircle," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFCircle,"\n");fflush(IMFCircle);IMFNCircle++;
          MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
          for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFCircle," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFCircle,"\n");fflush(IMFCircle);
          IMFNCircle++;
         }
       }
     }
   }

  if(1)
   {
    for(iBndChart=0;iBndChart<nBndCharts;iBndChart++)
     {
      chart=MFAtlasChartWithBoundary(A,iBndChart,e);
      P=MFChartPolytope(MFAtlasChart(A,chart,e),e);

      R=MFChartRadius(MFAtlasChart(A,chart,e),e);
      for(i=0;i<MFPolytopeNumberOfVertices(P,e);i++)
       {
        n1=MFPolytopeVertexNumberOfIndices(P,i,e);
        for(j=i+1;j<MFPolytopeNumberOfVertices(P,e);j++)
         {
          n2=MFPolytopeVertexNumberOfIndices(P,j,e);
          m=n1;if(n2>m)m=n2;

          if(m>0)
           {
            n=MFPolytopeIntersectIndexSets(P,i,j,inter,e);
           }else{
            n=0;
           }

          noncube=1;for(l=0;l<n;l++)noncube=noncube&&inter[l]>m;

          if(n==k-1&&noncube)
           {
            for(ii=0;ii<k-1;ii++)
             {
              iface[ii]=-1;
              for(l=0;l<MFPolytopeNumberOfFaces(P,e);l++)
                if(MFPolytopeFaceIndex(P,l,e)==inter[ii])iface[ii]=l;
             }
  
            for(ii=0;ii<k-1;ii++)
             {
              b[ii]=MFPolytopeFaceOrigin(P,iface[ii],e);
              nn=MFKV_CStar(MFPolytopeFaceNormal(P,iface[ii],e),e);
              d0;for(l=0;l<k;l++)d0+=nn[l]*nn[l];
              for(l=0;l<k;l++)
               {
                a[ii+(k-1)*l]=nn[l]/d0;
                asave[ii+(k-1)*l]=a[ii+(k-1)*l];
               }
              b[ii]=b[ii]/d0;
             }

            for(ii=0;ii<k-1;ii++)
             {
              for(jj=0;jj<k-1;jj++)
               {
                aa[ii+(k-1)*jj]=0;
                for(l=0;l<k;l++)aa[ii+(k-1)*jj]+=a[ii+(k-1)*l]*a[jj+(k-1)*l];
               }
             }

            MFNullVectorSmall(k,a,nrm,e);
            MFSolveFullSmall(k-1,aa,b,e);

            for(ii=0;ii<k;ii++)
             {
              org[ii]=0.;
              for(l=0;l<k-1;l++)org[ii]+=asave[l+(k-1)*ii]*b[l];
             }

            r=R*R;for(l=0;l<k;l++)r-=org[l]*org[l];
            if(r>0)
             {
              r=sqrt(r);

              good0=1;
              good1=1;
              for(ii=0;ii<MFPolytopeNumberOfFaces(P,e);ii++)
               {
                doit=1;for(l=0;l<k-1;l++)if(ii==iface[l])doit=0;
                if(doit)
                 {
                  oo=MFPolytopeFaceOrigin(P,ii,e);
                  nn=MFKV_CStar(MFPolytopeFaceNormal(P,ii,e),e);
                  d0=-oo;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l])*nn[l];
                  if(d0>0)good0=0;
                  d1=-oo;for(l=0;l<k;l++)d1+=(org[l]-r*nrm[l])*nn[l];
                  if(d1>0)good1=0;
                 }
               }

              if(good0)
               {
                for(l=0;l<k;l++)f[l]=org[l]+r*nrm[l];
                MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
                for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterper," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterper,"\n");fflush(IMFInterper);
                IMFNInterper++;
                if(1||inter[0]<inter[1])
                 {
                  for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFCircle," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFCircle,"\n");fflush(IMFCircle);
                  IMFNCircle++;
                  IMFEvaluateFlow(F,vp,p0,pF,e);
                  MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);
                  d0=0.;for(l=0;l<k;l++)d0+=f[l]*f[l];
                  for(l=0;l<k;l++)f[l]=org[l]+r*nrm[l]+.2*f[l]/sqrt(d0);
                  MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
                  for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFCircle," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFCircle,"\n");fflush(IMFCircle);
                  IMFNCircle++;
                 }
                for(l=0;l<k;l++)f[l]=org[l]+r*nrm[l];
                MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);

                IMFEvaluateFlow(F,vp,p0,pF,e);
                MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);
                d1=0.;for(l=0;l<k;l++)d1+=f[l]*nrm[l];
                t=-r/d1;
                if(t<0.)
                 {
                  good0=1;
                  for(ii=0;ii<k-1;ii++)
                   {
                    d0=0.;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l]+t*f[l])*asave[ii+(k-1)*l];
                    d1=0.;for(l=0;l<k;l++)d1+= org[l]                 *asave[ii+(k-1)*l];
                    if(d0<0 || d0>d1)good0=0;
                   }
                  if(good0)
                   {
                    for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterper," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterper,"\n");fflush(IMFInterper);
                    IMFNInterper++;
                   }
                 }
               }

              if(good1)
               {
                for(l=0;l<k;l++)f[l]=org[l]-r*nrm[l];
                MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
                for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterper," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterper,"\n");fflush(IMFInterper);
                IMFNInterper++;
                if(inter[0]<inter[1])
                 {
                  for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFCircle," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFCircle,"\n");fflush(IMFCircle);
                  IMFNCircle++;
                  IMFEvaluateFlow(F,vp,p0,pF,e);
                  MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);
                  d0=0.;for(l=0;l<k;l++)d0+=f[l]*f[l];
                  for(l=0;l<k;l++)f[l]=org[l]-r*nrm[l]+.2*f[l]/sqrt(d0);
                  MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
                  for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFCircle," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFCircle,"\n");fflush(IMFCircle);
                  IMFNCircle++;
                 }
                for(l=0;l<k;l++)f[l]=org[l]-r*nrm[l];
                MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);

                IMFEvaluateFlow(F,vp,p0,pF,e);
                MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);
                d1=0.;for(l=0;l<k;l++)d1+=f[l]*nrm[l];
                t=r/d1;
                if(t<0.)
                 {
                  good1=1;
                  for(ii=0;ii<k-1;ii++)
                   {
                    d0=0.;for(l=0;l<k;l++)d0+=(org[l]-r*nrm[l]+t*f[l])*asave[ii+(k-1)*l];
                    d1=0.;for(l=0;l<k;l++)d1+= org[l]                 *asave[ii+(k-1)*l];
                    if(d0<0 || d0>d1)good1=0;
                   }
                  if(good1)
                   {
                    for(l=0;l<MFAtlasN(A,e);l++)fprintf(IMFInterper," %lf",(MFNV_CStar(vp,e))[l]);fprintf(IMFInterper,"\n");fflush(IMFInterper);
                    IMFNInterper++;
                   }
                 }
               }
             }
           }
         }
       }
     }
   }
#endif

  for(iBndChart=0;iBndChart<nList;iBndChart++)
   {
    chart=list[iBndChart];
    P=MFChartPolytope(MFAtlasChart(A,chart,e),e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("The %dth chart on the boundary list is chart %d\n",iBndChart,chart);fflush(stdout);}
#endif

    if(!MFChartHasBoundary(MFAtlasChart(A,chart,e),e))continue;
    nv=MFPolytopeNumberOfVertices(P,e);
/*  if(!MFChartChangedFlag(MFAtlasChart(A,chart,e),e))nv=0;*/

    R=MFChartRadius(MFAtlasChart(A,chart,e),e);
    for(i=0;i<nv;i++)
     {
      n1=MFPolytopeVertexNumberOfIndices(P,i,e);
      for(j=i+1;j<nv;j++)
       {
        n2=MFPolytopeVertexNumberOfIndices(P,j,e);
        m=n1;if(n2>m)m=n2;

        if(m>0)
         {
          n=MFPolytopeIntersectIndexSets(P,i,j,inter,e);
         }else{
          n=0;
         }

        noncube=1;for(l=0;l<n;l++)noncube=noncube&&inter[l]>m;

        if(n==k-1&&noncube)
         {

#ifdef MFALLOWVERBOSE
          if(verbose){printf("  Indices of this non-cube edge:\n");fflush(stdout);}
#endif
          for(ii=0;ii<k-1;ii++)
           {
            iface[ii]=-1;
            for(l=0;l<MFPolytopeNumberOfFaces(P,e);l++)
             {
              if(MFPolytopeFaceIndex(P,l,e)==inter[ii])iface[ii]=l;
             }

#ifdef MFALLOWVERBOSE
            if(verbose){printf("  %d (face %d)\n",inter[ii],iface[ii]);fflush(stdout);}
#endif
           }

          for(ii=0;ii<k-1;ii++)
           {
            b[ii]=MFPolytopeFaceOrigin(P,iface[ii],e);
            nn=MFKV_CStar(MFPolytopeFaceNormal(P,iface[ii],e),e);
            d0;for(l=0;l<k;l++)d0+=nn[l]*nn[l];
            for(l=0;l<k;l++)
             {
              a[ii+(k-1)*l]=nn[l]/d0;
              asave[ii+(k-1)*l]=a[ii+(k-1)*l];
             }
            b[ii]=b[ii]/d0;
           }

          for(ii=0;ii<k-1;ii++)
           {
            for(jj=0;jj<k-1;jj++)
             {
              aa[ii+(k-1)*jj]=0;
              for(l=0;l<k;l++)aa[ii+(k-1)*jj]+=a[ii+(k-1)*l]*a[jj+(k-1)*l];
             }
           }

          MFNullVectorSmall(k,a,nrm,e);
          MFSolveFullSmall(k-1,aa,b,e);

          for(ii=0;ii<k;ii++)
           {
            org[ii]=0.;
            for(l=0;l<k-1;l++)org[ii]+=asave[l+(k-1)*ii]*b[l];
           }

#ifdef MFALLOWVERBOSE
          if(verbose)
           {
            printf("  Done finding edge equation org=(%lf",org[0]);fflush(stdout);
            for(ii=1;ii<k;ii++)printf(",%lf",org[ii]);
            printf(")\n");fflush(stdout);
            printf("                             nrm=(%lf",nrm[0]);fflush(stdout);
            for(ii=1;ii<k;ii++)printf(",%lf",nrm[ii]);
            printf(")\n");fflush(stdout);
    
            printf("Check: is org on all the faces?\n");fflush(stdout);
            d0=0.;for(l=0;l<k;l++)d0+=org[l]*nrm[l];
            printf("            nrm.org=%lf\n",d0);fflush(stdout);
            for(ii=0;ii<k-1;ii++)
             {
              oo=MFPolytopeFaceOrigin(P,iface[ii],e);
              nn=MFKV_CStar(MFPolytopeFaceNormal(P,iface[ii],e),e);
              d0=-oo;for(l=0;l<k;l++)d0+=nn[l]*org[l];
              printf("   face %d: org.n-o=%lf\n",ii,d0);fflush(stdout);
              d0=0.;for(l=0;l<k;l++)d0+=nn[l]*nrm[l];
              printf("            nrm.n=%lf\n",d0);fflush(stdout);
             }
           }
#endif

          r=R*R;for(l=0;l<k;l++)r-=org[l]*org[l];
          if(r>0)
           {
            r=sqrt(r);

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              d0=0.;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l])*(org[l]+r*nrm[l]);d0=sqrt(d0);
              printf("  This edge crosses the ball. |pt0|=%lf, R=%lf\n",d0,R);fflush(stdout);
              d0=0.;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l])*(org[l]+r*nrm[l]);d0=sqrt(d0);
              printf("                              |pt1|=%lf, R=%lf\n",d0,R);fflush(stdout);
             }
#endif

            good0=1;
            good1=1;
            for(ii=0;ii<MFPolytopeNumberOfFaces(P,e);ii++)
             {
              doit=1;for(l=0;l<k-1;l++)if(ii==iface[l])doit=0;
              if(doit)
               {
                oo=MFPolytopeFaceOrigin(P,ii,e);
                nn=MFKV_CStar(MFPolytopeFaceNormal(P,ii,e),e);
                d0=-oo;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l])*nn[l];
                if(d0>0)good0=0;
                d1=-oo;for(l=0;l<k;l++)d1+=(org[l]-r*nrm[l])*nn[l];
                if(d1>0)good1=0;
               }
             }

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              if(good0){printf("  The +r point is in the polytope.\n");fflush(stdout);}
                  else {printf("  The +r point is not in the polytope.\n");fflush(stdout);}
             }
#endif

            if(good0)
             {
              for(l=0;l<k;l++)f[l]=org[l]+r*nrm[l];
              MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
              IMFEvaluateFlow(F,vp,p0,pF,e);
              MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);
              d1=0.;for(l=0;l<k;l++)d1+=f[l]*nrm[l];
              t=-r/d1;
              if(t<0.&&(Omega==NULL||MFNRegionInterior(Omega,vp,e),e))
               {
                good0=1;
                for(ii=0;ii<k-1;ii++)
                 {
                  d0=0.;for(l=0;l<k;l++)d0+=(org[l]+r*nrm[l]+t*f[l])*asave[ii+(k-1)*l];
                  d1=0.;for(l=0;l<k;l++)d1+= org[l]                 *asave[ii+(k-1)*l];
                  if(d0<0 || d0>d1)good0=0;
                 }
                if(1||good0)
                 {
                  for(l=0;l<k;l++)f[l]=org[l]+r*nrm[l]+t*f[l];
                  MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,z,e);
                  result=IMFInterpolateOnDualFace(A,F,p0,z,chart,inter,e);
                  if(result!=NULL)
                   {
                    if(IMFExpansionNVGetT(result,e)<tmax&&(Omega==NULL||MFNRegionInterior(Omega,IMFExpansionNVGetSigma(result,e),e)))
                     {
#ifdef DOINTERPANIM
                      for(j=0;j<MFNV_NC(vp,e);j++)fprintf(IMFInterpee," %lf",(MFNV_CStar(vp,e))[j]);fprintf(IMFInterpee,"\n");fflush(IMFInterper);
                      IMFNInterpee++;
                      for(j=0;j<MFNV_NC(result,e);j++)fprintf(IMFInterpee," %lf",IMFExpansionU((IMFExpansionNVGetE(result,e)),e)[j]);fprintf(IMFInterpee,"\n");fflush(IMFInterper);
                      IMFNInterpee++;
#endif

#ifdef MFALLOWVERBOSE
                      if(verbose){printf("Result is on chart %d\n",chart);fflush(stdout);}
#endif

                      MFChartResetChangedFlag(MFAtlasChart(A,chart,e),e);
                      goto FreeAndReturn;
                     }
                   }
                 }
               }
             }

#ifdef MFALLOWVERBOSE
            if(verbose)
             {
              if(good1){printf("  The -r point is in the polytope.\n");fflush(stdout);}
                  else {printf("  The +r point is not in the polytope.\n");fflush(stdout);}
             }
#endif

            if(good1)
             {
              for(l=0;l<k;l++)f[l]=org[l]-r*nrm[l];
              MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,vp,e);
              IMFEvaluateFlow(F,vp,p0,pF,e);
              MFChartProjectIntoTangentSpace(MFAtlasChart(A,chart,e),vF,vf,e);
              d1=0.;for(l=0;l<k;l++)d1+=f[l]*nrm[l];
              t=r/d1;
              if(t<0.&&(Omega==NULL||MFNRegionInterior(Omega,vp,e)))
               {
                good1=1;
                for(ii=0;ii<k-1;ii++)
                 {
                  d0=0.;for(l=0;l<k;l++)d0+=(org[l]-r*nrm[l]+t*f[l])*asave[ii+(k-1)*l];
                  d1=0.;for(l=0;l<k;l++)d1+= org[l]                 *asave[ii+(k-1)*l];
                  if(d0<0 || d0>d1)good1=0;
                 }
                if(1||good1)
                 {
                  for(l=0;l<k;l++)f[l]=org[l]-r*nrm[l]+t*f[l];
                  MFChartPointInTangentSpace(MFAtlasChart(A,chart,e),vf,z,e);
                  result=IMFInterpolateOnDualFace(A,F,p0,z,chart,inter,e);
                  if(result!=NULL)
                   {
                    if(IMFExpansionNVGetT(result,e)<tmax&&(Omega==NULL||MFNRegionInterior(Omega,IMFExpansionNVGetSigma(result,e),e)))
                     {
#ifdef DOINTERPANIM
                      for(j=0;j<MFNV_NC(vp,e);j++)fprintf(IMFInterpee," %lf",(MFNV_CStar(vp,e))[j]);fprintf(IMFInterpee,"\n");fflush(IMFInterper);
                      IMFNInterpee++;
                      for(j=0;j<MFNV_NC(result,e);j++)fprintf(IMFInterpee," %lf",IMFExpansionU((IMFExpansionNVGetE(result,e)),e)[j]);fprintf(IMFInterpee,"\n");fflush(IMFInterper);
                      IMFNInterpee++;
#endif

#ifdef MFALLOWVERBOSE
                      if(verbose){printf("Result is on chart %d\n",chart);fflush(stdout);}
#endif

                      MFChartResetChangedFlag(MFAtlasChart(A,chart,e),e);
                      goto FreeAndReturn;
                     }
                   }
                 }
               }
             }
           }
         }
       }
     }
    MFChartResetChangedFlag(MFAtlasChart(A,chart,e),e);
   }

  result=NULL;
  printf("No such points on the boundary!!!\n");

FreeAndReturn:
  MFFreeKVector(vf,e);
  MFFreeNVector(z,e);
  MFFreeNVector(vp,e);
  MFFreeNVector(vF,e);

  return result;
 }

#define nmax 10000

IMFFlow IMFTPBVPFlow=NULL;
MFChart IMFTPBVPChart=NULL;
MFAtlas IMFTPBVPC=NULL;
int IMFTPBVPCChart=-1;
double *IMFTPBVPs=NULL;

static void IMFTPBVPf (double,int,double*,int,double*,double*,double*,double*,MFErrorHandler);
static void IMFTPBVPfu(double,int,double*,int,double*,double*,double*,double*,MFErrorHandler);
static void IMFTPBVPfl(double,int,double*,int,double*,double*,double*,double*,MFErrorHandler);
static void IMFTPBVPa (int,int,double*,double*,int,double*,double*,double*,double*,double*,MFErrorHandler);
static void IMFTPBVPau(int,int,double*,double*,int,double*,double*,double*,double*,double*,MFErrorHandler);
static void IMFTPBVPal(int,int,double*,double*,int,double*,double*,double*,double*,double*,MFErrorHandler);

MFNVector IMFInterpolateOnDualFace2(MFAtlas A, IMFFlow F, MFKVector p0, MFAtlas c, MFKVector s, MFNVector z, int chart, int *inter, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFInterpolateOnDualFace2"};
  int i,j,l,n,m;
  int prev;
  MFNVector result;
  MFNVector u[nmax];
  int k;
  double a;
  MFImplicitMF M;
  MFKVector s0;
  MFNVector u0,ug;
  double lambda;
  IMFExpansion E;
  MFNKMatrix TS;
  double R,Rf,Rmax;
  double epsilon=.01;
  int nx,nu,np;
  int nE;
  int verbose=0;
  IMFFlow Fp;
  IMFExpansion E0,E1;
  MFNVector c0,c1;
  double d;
  int good;
  int nsigma;
  MFNVector diffi;
  MFNVector diffj;
  MFNVector diffz;
  double t[kmax];
  MFNVector C[kmax];
  MFNSpace space;
  double aa[(kmax-1)*(kmax-1)];
 
#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s-------------------------------------------\n",RoutineName);fflush(stdout);}
#endif

  k=MFAtlasK(A,e);
  n=MFAtlasN(A,e);

  if(k>kmax)
   {
    fprintf(stderr,"%s, I statically allocated some things here for k<10! Complain.\n");
    fflush(stderr);
    abort();
   }

  C[0]=MFAtlasChartCenter(A,chart,e);
  printf("Chart containing the edge is %d",chart);MFPrintNVector(stdout,C[0],e);printf("\n");fflush(stdout);
  nsigma=MFNV_NC(IMFExpansionNVGetSigma(C[0],e),e);
  good=0;
  for(i=0;i<k-1;i++)
   {
    j=MFAtlasHalfSpaceRightChart(A,inter[i],e);
    if(j==chart)j=MFAtlasHalfSpaceLeftChart(A,inter[i],e);
    C[i+1]=MFAtlasChartCenter(A,j,e);
    printf(" adj Chart %d is %d",i,j);fflush(stdout);MFPrintNVector(stdout,C[i+1],e);printf("\n");fflush(stdout);
    d=0.;for(j=0;j<nsigma;j++)d+=pow(MFNV_CStar(IMFExpansionNVGetSigma(C[0],e),e)[j]-
                                     MFNV_CStar(IMFExpansionNVGetSigma(C[i+1],e),e)[j],2.);
    if(sqrt(d)>1.e-14)good=1;
   }

  if(!good)return NULL;

  diffi=MFCreateNVector(n,e);
  diffj=MFCreateNVector(n,e);
  diffz=MFCreateNVector(n,e);
  space=MFIMFNSpace(MFAtlasMF(A,e),e);

  MFNSpaceDirection(space,z,C[0],diffz,e);
  for(i=0;i<k-1;i++)
   {
    MFNSpaceDirection(space,C[i+1],C[0],diffi,e);
    aa[i+(k-1)*i]=MFNSpaceInner(space,diffi,diffi,e);
    for(j=i+1;j<k-1;j++)
     {
      MFNSpaceDirection(space,C[j+1],C[0],diffj,e);
      aa[i+(k-1)*j]=MFNSpaceInner(space,diffi,diffj,e);
      aa[i+(k-1)*j]=aa[j+(k-1)*i];
     }
    t[i+1]=MFNSpaceInner(space,diffz,diffi,e);
   }

  MFSolveFullSmall(k-1,aa,t+1,e);

  t[0]=1.;for(i=1;i<k;i++)t[0]-=t[i];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Barycentric coordinates (%lf",t[0]);
    for(i=1;i<k;i++)printf(",%lf",t[i]);
    printf(")\n");
   }
#endif

  MFFreeNVector(diffi,e);
  MFFreeNVector(diffj,e);
  MFFreeNVector(diffz,e);

  good=1;
  for(i=0;i<k;i++)if(t[i]<0.||t[i]>1.)good=0;
  for(i=0;i<k;i++)if(t[i]<.1||t[i]>.9)good=0;  /* EXPERIMENT */
  if(!good)return NULL;

#ifdef DOINTERPANIM
  E0=IMFExpansionNVGetE(MFAtlasChartCenter(A,chart,e),e);
  for(j=0;j<3;j++)fprintf(IMFInterpee," %lf",(IMFExpansionU(E0,e))[j]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);IMFNInterpee++;
  for(i=0;i<2-1;i++)
   {
    m=MFAtlasHalfSpaceRightChart(A,inter[i],e);
    if(m==chart)m=MFAtlasHalfSpaceLeftChart(A,inter[i],e);
   }
  E1=IMFExpansionNVGetE(MFAtlasChartCenter(A,m,e),e);
  for(j=0;j<3;j++)fprintf(IMFInterpee," %lf",(IMFExpansionU(E1,e))[j]);fprintf(IMFInterpee,"\n");fflush(IMFInterpee);IMFNInterpee++;
#endif

/* Get trajectory from the initial curve to the center of chart */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Get trajectory from the initial curve to the center of chart\n");fflush(stdout);}
#endif

  n=0;
  prev=IMFExpansionNVGetPrev(MFAtlasChartCenter(A,chart,e),e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf(" n=%d, chart=%d, prev=%d\n",n,chart,prev);fflush(stdout);}
#endif

  n++;
  while(prev>0)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf(" n=%d, chart=%d",n,prev);}
#endif

    prev=IMFExpansionNVGetPrev(MFAtlasChartCenter(A,prev,e),e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf(", prev=%d\n",prev);fflush(stdout);}
#endif

    n++;
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("n=%d\n",n);fflush(stdout);}
#endif

  i=n-1;
  u[i]=MFCloneNVector(MFAtlasChartCenter(A,chart,e),e);
  prev=IMFExpansionNVGetPrev(u[i],e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("u[%d]=chart %d, prev=%d\n",i,chart,prev);fflush(stdout);}
#endif

  i--;
  while(prev>0)
   {

#ifdef MFALLOWVERBOSE
    if(verbose){printf("u[%d]=chart %d",i,prev);fflush(stdout);}
#endif

    u[i]=MFCloneNVector(MFAtlasChartCenter(A,prev,e),e);
    prev=IMFExpansionNVGetPrev(u[i],e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf(", prev=%d\n",prev);fflush(stdout);}
#endif

    i--;
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Done getting trajectory\n");fflush(stdout);}
#endif

/* Do the homotopy for the interpolated trajectory */

  k=MFAtlasK(c,e)+1;
  Fp=IMFCreateFatFlow(F,k,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("c(0x%8.8x) is an %d-manifold in %d space\n",c,MFAtlasK(c,e),MFAtlasN(c,e));fflush(stdout);}
#endif

  m=1;
  nx=m*(n-1)+2;
  M=MFIMFCreateTPBVP(1,nx,IMFFlowNU(Fp,e),k+1,IMFTPBVPf,IMFTPBVPfu,IMFTPBVPfl,
                            IMFFlowNU(Fp,e)+k,IMFTPBVPa,IMFTPBVPau,IMFTPBVPal,
                            0,NULL,NULL,NULL,
                            NULL,NULL,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("trajectory starts at chart %d on c, s(0x%8.8x)=",IMFExpansionNVGetChart0(u[0],e),IMFExpansionNVGetS0(u[0],e));fflush(stdout);}
  if(verbose){MFPrintKVector(stdout,IMFExpansionNVGetS0(u[0],e),e);printf("\n");fflush(stdout);}
#endif

  IMFTPBVPFlow=Fp;
  IMFTPBVPChart=MFAtlasChart(A,chart,e);
  IMFFlowP0=p0;
  k=MFKV_NC(IMFExpansionNVGetS0(u[0],e),e)+1;
  IMFTPBVPs=(double*)malloc((k-1)*sizeof(double));if(IMFTPBVPs==NULL)abort();

#ifndef MFNOSAFETYNET
  if(IMFTPBVPs==NULL)
   {
    sprintf(IMFInterpolationErrorMsg,"Out of memory, trying to allocate %d bytes",(k-1)*sizeof(double));
    MFSetError(e,12,RoutineName,IMFInterpolationErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<k-1;i++)IMFTPBVPs[i]=MFKV_C(IMFExpansionNVGetS0(u[0],e),i,e);
  IMFTPBVPC=c;
  IMFTPBVPCChart=IMFExpansionNVGetChart0(u[0],e);

  nx=MFTPBVPGetNX(M,e);
  nu=MFTPBVPGetNU(M,e);
  np=MFTPBVPGetNP(M,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("TPBVP: n=%d, nx=%d, nu=%d, np=%d\n",MFIMF_N(M,e),nx,nu,np);fflush(stdout);}
#endif

  TestEm(nu,np,MFTPBVPGetNBC(M,e),e);

  TS=MFCreateNKMatrixWithData(MFIMF_N(M,e),1,NULL,e);
  for(i=0;i<MFIMF_N(M,e);i++)MFNKMSetC(TS,i,0,0.,e);
  MFNKMSetC(TS,nu*nx+k-1,0,1.,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Construct initial guess\n");fflush(stdout);}
#endif

  u0=MFCreateNVector(nu*nx+np+nx+1,e);

  E0=IMFExpansionNVGetE(u[0],e);E1=IMFExpansionNVGetE(u[1],e);
  a=-.5/m;
  printf("E0 dataLn=%d\n",IMFExpansionDataLn(E0,e));fflush(stdout);
  printf("E1 dataLn=%d\n",IMFExpansionDataLn(E1,e));fflush(stdout);
  nE=IMFExpansionDataLn(E1,e);
  for(j=0;j<nE;j++)MFNVSetC(u0,j,(1-a)*IMFExpansionData(E0,e)[j]+a*IMFExpansionData(E1,e)[j],e);
  MFNVSetC(u0,nu*nx+np,(1-a)*IMFExpansionNVGetT(u[0],e)+a*IMFExpansionNVGetT(u[1],e),e);
  printf("left most point, t=%lf, u=(%lf,%lf,%lf,...)\n",MFNV_C(u0,nu*nx+np,e),MFNV_C(u0,0,e),MFNV_C(u0,1,e),MFNV_C(u0,2,e));
    fflush(stdout);
#ifdef DOTRAJ
  for(j=0;j<3;j++)fprintf(IMFInterpT," %lf",MFNV_C(u0,j,e));fprintf(IMFInterpT,"\n");fflush(IMFInterpT);IMFNInterpT++;
#endif

  for(i=1;i<n;i++)
   {
    E0=IMFExpansionNVGetE(u[i-1],e);
    E1=IMFExpansionNVGetE(u[i  ],e);
    printf("  u[%d] t=%lf, u=(%lf,%lf,%lf,...)\n",i-1,IMFExpansionNVGetT(u[i-1],e),IMFExpansionData(E0,e)[0],IMFExpansionData(E0,e)[1],IMFExpansionData(E0,e)[2]);
    fflush(stdout);
    for(l=0;l<m;l++)
     {
      a=(.5+l)/m;
      for(j=0;j<nE;j++)MFNVSetC(u0,j+nu*(1+l+m*(i-1)),(1-a)*IMFExpansionData(E0,e)[j]+a*IMFExpansionData(E1,e)[j],e);
#ifdef DOTRAJ
  for(j=0;j<3;j++)fprintf(IMFInterpT," %lf",MFNV_C(u0,j+nu*(l+m*(i-1)),e));fprintf(IMFInterpT,"\n");fflush(IMFInterpT);IMFNInterpT++;
#endif
      MFNVSetC(u0,nu*nx+np+1+l+m*(i-1),(1-a)*IMFExpansionNVGetT(u[i-1],e)+a*IMFExpansionNVGetT(u[i],e),e);
      printf("   l=%d, t=%lf, u=(%lf,%lf,%lf,...)\n",l,MFNV_C(u0,nu*nx+np+l+m*(i-1),e),MFNV_C(u0,0+nu*(l+m*(i-1)),e),MFNV_C(u0,1+nu*(l+m*(i-1)),e),MFNV_C(u0,2+nu*(l+m*(i-1)),e));
        fflush(stdout);
     }
   }
  E0=IMFExpansionNVGetE(u[n-1],e);
  printf("  u[%d] t=%lf, u=(%lf,%lf,%lf,...)\n",n-1,IMFExpansionNVGetT(u[n-1],e),IMFExpansionData(E0,e)[0],IMFExpansionData(E0,e)[1],IMFExpansionData(E0,e)[2]);

  a=(.5+m+1)/m;
  E0=IMFExpansionNVGetE(u[n-2],e);
  E1=IMFExpansionNVGetE(u[n-1],e);
  for(j=0;j<nE;j++)MFNVSetC(u0,j+nu*(1+m*(n-1)),(1-a)*IMFExpansionData(E0,e)[j]+a*IMFExpansionData(E1,e)[j],e);
  MFNVSetC(u0,nu*nx+np+m*(n-1)+1,(1-a)*IMFExpansionNVGetT(u[n-2],e)+a*IMFExpansionNVGetT(u[n-1],e),e);
  printf("right most point, t=%lf, u=(%lf,%lf,%lf,...)\n",MFNV_C(u0,nu*nx+np+1+m*(n-1),e),MFNV_C(u0,0+nu*(1+m*(n-1)),e),MFNV_C(u0,1+nu*(1+m*(n-1)),e),MFNV_C(u0,2+nu*(1+m*(n-1)),e));
        fflush(stdout);
#ifdef DOTRAJ
  for(j=0;j<3;j++)fprintf(IMFInterpT," %lf",MFNV_C(u0,j+nu*(1+m*(n-1)),e));fprintf(IMFInterpT,"\n");fflush(IMFInterpT);IMFNInterpT++;
#endif

  MFNVSetC(u0,nu*nx+0,1.,e); /* alpha */
  for(j=0;j<k-1;j++)
   {
    MFNVSetC(u0,nu*nx+j+1,0.,e); /* sigma_j */
   }
  MFNVSetC(u0,nu*nx+k,0.,e); /* lambda */

  ug=MFCreateNVector(MFIMF_N(M,e),e);

/* Simple homotopy */

  for(lambda=0.;lambda<=1.;lambda+=.1);
   {
    MFNVSetC(u0,nu*nx+k,lambda,e);

#ifdef MFALLOWVERBOSE
    if(verbose){printf("Project at lambda=%lf\n",lambda);fflush(stdout);}
#endif

    MFIMFProject(M,u0,TS,ug,e);
    for(i=0;i<MFIMF_N(M,e);i++)MFNVSetC(u0,i,MFNV_C(ug,i,e),e);
   }

/* Add the interpolated trajectory to the atlas */

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Add the interpolated trajectory to the atlas\n");fflush(stdout);}
#endif

  IMFExpansionNVSetChart0(u[0],IMFTPBVPCChart,e);
  s0=MFCreateKVector(k-1,e);
  for(i=0;i<k-1;i++)MFKVSetC(s0,i,IMFTPBVPs[i],e);
  IMFExpansionNVSetS0(u[0],s0,e);
  MFFreeKVector(s0,e);
  free(IMFTPBVPs);

  prev=-1;
  for(i=0;i<n;i++)
   {
    IMFExpansionNVSetPrev(u[i],prev,e);
    if(i>0)IMFExpansionNVSetType(u[i],3,e);
     else IMFExpansionNVSetType(u[i],1,e);
    E=IMFExpansionNVGetE(u[i],e);
    TS=IMFExpansionTS(E,e);
    R=IMFExpansionR(E,epsilon,e);
    Rf=IMFFlowR(F,epsilon,u[i],p0,TS,e);
    if(Rf<R)R=Rf;
    if(R>Rmax||R!=R)R=Rmax;
    prev=MFAtlasAddChartWithAll(A,u[i],TS,R,e);
   }
  
#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  IMFFreeFlow(Fp,e);
  return u[n-1];
 }

void IMFTPBVPf(double r, int nu, double *u, int np, double *p, double *u0, double *l0, double *f, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFTPBVPf"};
  int i;
  MFNVector vu;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("In IMFTPBVPf\n");fflush(stdout);
    printf("nu is %d, FlowN is %d\n",nu,IMFFlowNU(IMFTPBVPFlow,e),e);fflush(stdout);
   }
#endif

  vu=MFCreateWrappedNVector(nu,u,e);

  IMFEvaluateFlow(IMFTPBVPFlow,vu,IMFFlowP0,f,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("p[0]=%lf,f:(%lf,",p[0],f[0]);for(i=1;i<nu;i++)printf(",%lf",f[i]);printf(")\n");fflush(stdout);}
#endif

  for(i=0;i<nu;i++)f[i]=p[0]*f[i];

  MFFreeNVector(vu,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Done IMFTPBVPf\n");fflush(stdout);}
#endif

  return;
 }

void IMFTPBVPfu(double r, int nu, double *u, int np, double *p, double *u0, double *l0, double *fu, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFTPBVPfu"};
  int i;
  MFNVector vu;

  vu=MFCreateWrappedNVector(nu,u,e);

  IMFEvaluateDerivativeOfFlow(IMFTPBVPFlow,vu,IMFFlowP0,fu,e);
  for(i=0;i<nu*nu;i++)fu[i]=p[0]*fu[i];

  MFFreeNVector(vu,e);

  return;
 }

void IMFTPBVPfl(double r, int nu, double *u, int np, double *p, double *u0, double *l0, double *fl, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFTPBVPfl"};
  int i;
  MFNVector vu;

  vu=MFCreateWrappedNVector(nu,u,e);
  IMFEvaluateFlow(IMFTPBVPFlow,vu,IMFFlowP0,fl,e);
  for(i=0;i<nu;i++)fl[i+nu*0]=fl[i];
  for(i=0;i<nu*(np-1);i++)fl[i+nu]=0.;

  MFFreeNVector(vu,e);

  return;
 }

/* Parms:
      p[0] : alpha
      p[1],...,p[k-1] : sigma
      p[k] : lambda
*/

void IMFTPBVPa(int nbc, int nu, double *uL, double *uR, int np, double *p, double *u0L, double *u0R,double *l0, double *a, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFTPBVPa"};
  int i,k;
  MFKVector vsigma;
  MFKVector vs;
  MFNVector vuR;
  MFNVector vc;
  double *c;
  MFChart cchart;
  double nrm;
  double *sigma,*s;
  IMFExpansion E;
  static int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("In IMFTPBVPa\n");fflush(stdout);}
#endif

  k=MFAtlasK(IMFTPBVPC,e);
  vsigma=MFCreateKVector(k,e);
  sigma=MFKV_CStar(vsigma,e);
  vs=MFCreateKVector(k+1,e);
  s=MFKV_CStar(vs,e);
  vuR=MFCreateWrappedNVector(nu,uR,e);
  c=(double*)malloc(MFAtlasN(IMFTPBVPC,e)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(c==NULL)
   {
    sprintf(IMFInterpolationErrorMsg,"Out of memory, trying to allocate %d bytes",MFAtlasN(IMFTPBVPC,e)*sizeof(double));
    MFSetError(e,12,RoutineName,IMFInterpolationErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<MFAtlasN(IMFTPBVPC,e);i++)c[i]=0.;
  vc=MFCreateWrappedNVector(MFAtlasN(IMFTPBVPC,e),c,e);

  cchart=MFAtlasChart(IMFTPBVPC,IMFTPBVPCChart,e);
  nrm=0.;for(i=0;i<k-1;i++)nrm+=p[i]*p[i];nrm=sqrt(nrm);
  if(nrm>MFAtlasChartRadius(IMFTPBVPC,IMFTPBVPCChart,e))
   {
    abort();
/*  pivot to next, change cchart and s
    this is going to make it difficult to take derivatives.*/
   }

  for(i=0;i<k-1;i++)sigma[i]=p[i+1];

/* c is n+k, look up last k -> s on expansion, evaluate expansion and derivatives at this point to 
     get a vector of lenght nu */

  MFChartPointInTangentSpace(cchart,vsigma,vc,e);
  E=IMFSphereOnExpansionGetLocal(IMFTPBVPC,vc,e);

/* At t=0. */

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Conditions at t=0. (0..%d):\n",nu-1);
    printf(" uL=(%lf,%lf,%lf,...\n",uL[0],uL[1],uL[2]);fflush(stdout);
    printf(" s=(%lf), c(s)=(%lf,%lf,%lf), chart %d\n",c[4],c[0],c[1],c[2],cchart);fflush(stdout);
    printf(" E(c(s))=(%lf,%lf,%lf,...)\n",IMFExpansionData(E,e)[0],IMFExpansionData(E,e)[1],IMFExpansionData(E,e)[2]);fflush(stdout);
   }
#endif

  for(i=0;i<nu;i++)a[i]=uL[i];
  for(i=0;i<IMFExpansionDataLn(E,e);i++)a[i]-=IMFExpansionData(E,e)[i];

/* At t=1. */

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Conditions at t=1. (%d..%d):\n",nu,nbc);fflush(stdout);
   }
#endif

  MFChartProjectIntoTangentSpace(IMFTPBVPChart,vuR,vs,e);
  for(i=nu;i<nbc;i++)a[i]=s[i-nu]-p[k]*IMFTPBVPs[i-nu];
  
  MFFreeKVector(vsigma,e);
  MFFreeKVector(vs,e);
  MFFreeNVector(vuR,e);
  free(c);
  MFFreeNVector(vc,e);
  IMFFreeExpansion(E,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("Done IMFTPBVPa\n");fflush(stdout);}
#endif

  return;
 }

void IMFTPBVPau(int nbc, int nu, double *uL, double *uR, int np, double *p, double *u0L, double *u0R,double *l0, double *au, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFTPBVPau"};
  int i,j,k;
  MFNKMatrix cTS;
  double *TS;
  MFKVector vsigma;
  MFNVector vc;
  double *c;
  double *sigma;
  IMFExpansion E;
  MFChart cchart;
  int n1;

  k=MFAtlasK(IMFTPBVPC,e);
  for(i=0;i<2*nbc*nu;i++)au[i]=0.;

  cTS=MFChartTangentSpace(IMFTPBVPChart,e);
  TS=MFNKM_CStar(cTS,e);

  vsigma=MFCreateKVector(k,e);
  sigma=MFKV_CStar(vsigma,e);
  for(i=0;i<k-1;i++)sigma[i]=p[i+1];

  c=(double*)malloc(MFAtlasN(IMFTPBVPC,e)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(c==NULL)
   {
    sprintf(IMFInterpolationErrorMsg,"Out of memory, trying to allocate %d bytes",MFAtlasN(IMFTPBVPC,e)*sizeof(double));
    MFSetError(e,12,RoutineName,IMFInterpolationErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<MFAtlasN(IMFTPBVPC,e);i++)c[i]=0.;
  vc=MFCreateWrappedNVector(MFAtlasN(IMFTPBVPC,e),c,e);

  cchart=MFAtlasChart(IMFTPBVPC,IMFTPBVPCChart,e);
  MFChartPointInTangentSpace(cchart,vsigma,vc,e);
  E=IMFSphereOnExpansionGetLocal(IMFTPBVPC,vc,e);
  n1=MFChartN(IMFTPBVPChart,e);

  for(i=0;i<nu;i++)
      au[i+nbc*(i+0*nu)]=1.;

  for(i=nu;i<nbc;i++)
   {
    for(j=0;j<n1;j++)
      au[i+nbc*(j+nu)]=TS[j+n1*(i-nu)];
    for(j=n1;j<nu;j++)
      au[i+nbc*(j+nu)]=0.;
   }

  MFFreeKVector(vsigma,e);
  free(c);
  MFFreeNVector(vc,e);
  
  return;
 }

void IMFTPBVPal(int nbc, int nu, double *uL, double *uR, int np, double *p, double *u0L, double *u0R,double *l0, double *al, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFTPBVPal"};
  int i,j,k;
  MFNKMatrix cTS;
  double *TS;
  MFKVector vsigma;
  MFNVector vc;
  double *c;
  double *sigma;
  double *du;
  IMFExpansion E;
  int nc,n1;

  for(i=0;i<nbc*np;i++)al[i]=0.;

  k=MFAtlasK(IMFTPBVPC,e);
  nc=MFAtlasN(IMFTPBVPC,e);
  cTS=MFAtlasChartTangentSpace(IMFTPBVPC,IMFTPBVPCChart,e);
  TS=MFNKM_CStar(cTS,e);

  vsigma=MFCreateKVector(k,e);
  sigma=MFKV_CStar(vsigma,e);
  for(i=0;i<k-1;i++)sigma[i]=p[i+1];

  c=(double*)malloc(MFAtlasN(IMFTPBVPC,e)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(c==NULL)
   {
    sprintf(IMFInterpolationErrorMsg,"Out of memory, trying to allocate %d bytes",MFAtlasN(IMFTPBVPC,e)*sizeof(double));
    MFSetError(e,12,RoutineName,IMFInterpolationErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<MFAtlasN(IMFTPBVPC,e);i++)c[i]=0.;
  vc=MFCreateWrappedNVector(MFAtlasN(IMFTPBVPC,e),c,e);

  du=(double*)malloc(nu*sizeof(double));

#ifndef MFNOSAFETYNET
  if(du==NULL)
   {
    sprintf(IMFInterpolationErrorMsg,"Out of memory, trying to allocate %d bytes",nu*sizeof(double));
    MFSetError(e,12,RoutineName,IMFInterpolationErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif


  E=IMFSphereOnExpansionGetLocal(IMFTPBVPC,vc,e);
  n1=IMFExpansionDataLn(E,e);

  for(j=0;j<k;j++)
   {
    IMFEvaluateExpansionDirectionalDerivativeE(E,c+nc-k-1,TS+nc-k-1+nc*j,du,e);
    for(i=0;i<n1;i++)
      al[i+nbc*(j+1)]+=du[i];
   }

  for(i=nu;i<nbc;i++)
   al[i+nbc*k]=-IMFTPBVPs[i-nu];

  MFFreeKVector(vsigma,e);
  free(c);
  free(du);
  MFFreeNVector(vc,e);
  
  return;
 }

#define nbcMAX 100
#define nuMAX 100
#define npMAX 100

void TestA(int nbc, int nu, double *uL, double *uR, int np, double *p, double *u0L, double *u0R,double *l0,
                               void (*a)(int,int,double*,double*,int,double*,double*,double*,double*,double*,MFErrorHandler),
                               void (*au)(int,int,double*,double*,int,double*,double*,double*,double*,double*,MFErrorHandler),
                               void (*al)(int,int,double*,double*,int,double*,double*,double*,double*,double*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"TestA"};
  int i;
  int j;
  double Au[nbcMAX*2*nuMAX];
  double Al[nbcMAX*npMAX];
  double Ap[nbcMAX];
  double Am[nbcMAX];
  double save;
  double eps=1.e-4;
  double err;

  printf("Test derivatives of boundary conditions\n");fflush(stdout);

  (*au)(nbc,nu,uL,uR,np,p,u0L,u0R,l0,Au,e);
  (*al)(nbc,nu,uL,uR,np,p,u0L,u0R,l0,Al,e);
  for(i=0;i<nu;i++)
   {
    save=uL[i];
    uL[i]=save+eps;
    (*a)(nbc,nu,uL,uR,np,p,u0L,u0R,l0,Ap,e);
    uL[i]=save-eps;
    (*a)(nbc,nu,uL,uR,np,p,u0L,u0R,l0,Am,e);
    uL[i]=save;
    for(j=0;j<nbc;j++)
     {
      err=fabs(Au[j+nbc*i]-.5*(Ap[j]-Am[j])/eps);
      if(err>eps)
        printf(" da[%d]/duL[%d] = %le (exact) %le (diffed) error %le\n",j,i,Au[j+nbc*i],.5*(Ap[j]-Am[j])/eps,err);
     }
   }
  for(i=0;i<nu;i++)
   {
    save=uR[i];
    uR[i]=save+eps;
    (*a)(nbc,nu,uL,uR,np,p,u0L,u0R,l0,Ap,e);
    uR[i]=save-eps;
    (*a)(nbc,nu,uL,uR,np,p,u0L,u0R,l0,Am,e);
    uR[i]=save;
    for(j=0;j<nbc;j++)
     {
      err=fabs(Au[j+nbc*(i+nu)]-.5*(Ap[j]-Am[j])/eps);
      if(err>eps)
        printf(" da[%d]/duR[%d] = %le (exact) %le (diffed) error %le\n",j,i,Au[j+nbc*(i+nu)],.5*(Ap[j]-Am[j])/eps,err);
     }
   }
  for(i=0;i<np;i++)
   {
    save=p[i];
    p[i]=save+eps;
    (*a)(nbc,nu,uL,uR,np,p,u0L,u0R,l0,Ap,e);
    p[i]=save-eps;
    (*a)(nbc,nu,uL,uR,np,p,u0L,u0R,l0,Am,e);
    p[i]=save;
    for(j=0;j<nbc;j++)
     {
      err=fabs(Al[j+nbc*i]-.5*(Ap[j]-Am[j])/eps);
      if(err>eps)
        printf(" da[%d]/dp[%d] = %le (exact) %le (diffed) error %le\n",j,i,Al[j+nbc*i],.5*(Ap[j]-Am[j])/eps,err);
     }
   }
 }


void TestF(double r, int nu, double *u, int np, double *p, double *u0,double *l0,
                               void (*f)(double,int,double*,int,double*,double*,double*,double*,MFErrorHandler),
                               void (*fu)(double,int,double*,int,double*,double*,double*,double*,MFErrorHandler),
                               void (*fl)(double,int,double*,int,double*,double*,double*,double*,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"TestF"};
  int i;
  int j;
  double Fu[nuMAX*nuMAX];
  double Fl[nuMAX*npMAX];
  double Fp[nuMAX];
  double Fm[nuMAX];
  double save;
  double eps=1.e-4;
  double err;

  printf("Test derivatives of flow\n");fflush(stdout);

  (*fu)(r,nu,u,np,p,u0,l0,Fu,e);
  (*fl)(r,nu,u,np,p,u0,l0,Fl,e);
  for(i=0;i<nu;i++)
   {
    save=u[i];
    u[i]=save+eps;
    (*f)(r,nu,u,np,p,u0,l0,Fp,e);
    u[i]=save-eps;
    (*f)(r,nu,u,np,p,u0,l0,Fm,e);
    u[i]=save;
    for(j=0;j<nu;j++)
     {
      err=fabs(Fu[j+nu*i]-.5*(Fp[j]-Fm[j])/eps);
      if(err>eps)
        printf(" df[%d]/du[%d] = %le (exact) %le (diffed) error %le\n",j,i,Fu[j+nu*i],.5*(Fp[j]-Fm[j])/eps,err);
     }
   }
  for(i=0;i<np;i++)
   {
    save=p[i];
    p[i]=save+eps;
    (*f)(r,nu,u,np,p,u0,l0,Fp,e);
    p[i]=save-eps;
    (*f)(r,nu,u,np,p,u0,l0,Fm,e);
    p[i]=save;
    for(j=0;j<nu;j++)
     {
      eps=fabs(Fl[j+nu*i]-.5*(Fp[j]-Fm[j])/eps);
      if(err>eps)
        printf(" df[%d]/dp[%d] = %le (exact) %le (diffed) error %le\n",j,i,Fl[j+nu*i],.5*(Fp[j]-Fm[j])/eps,err);
     }
   }
 }

void TestEm(int nu, int np, int nbc, MFErrorHandler e)
 {
  static char RoutineName[]={"TestEm"};
  int i;
  double r;
  double u[nuMAX];
  double u0[nuMAX];
  double uL[nuMAX];
  double u0L[nuMAX];
  double uR[nuMAX];
  double u0R[nuMAX];
  double p[npMAX];
  double l0[npMAX];

  r=1.;
  for(i=0;i<nu;i++)
   {
    u[i]=1.;
    uL[i]=1.;
    uR[i]=1.;
    u0[i]=1.;
    u0L[i]=1.;
    u0R[i]=1.;
   }
  for(i=0;i<np;i++)
   {
    p[i]=1.;
    l0[i]=1.;
   }

  printf("Test derivatives\n");fflush(stdout);
  TestF(r,nu,u,np,p,u0,l0,IMFTPBVPf,IMFTPBVPfu,IMFTPBVPfl,e);
  TestA(nbc,nu,uL,uR,np,p,u0L,u0R,l0,IMFTPBVPa,IMFTPBVPau,IMFTPBVPal,e);
 }

#ifdef __cplusplus
}
#endif
