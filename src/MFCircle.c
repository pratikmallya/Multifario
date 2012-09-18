/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   August 26, 1999
 */

static char *id="@(#) $Id: MFCircle.c,v 1.4 2011/07/21 17:42:46 mhender Exp $";

#include <MFImplicitMF.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

extern double MFEpsilon;

static char MFCircleMFErrorHandlerMsg[256]="";

static void MFFreeCircleData(void*,MFErrorHandler);
static int MFProjectCircle(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentCircle(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static double MFScaleCircle(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static void MFWriteCircleData(FILE*,void*,MFErrorHandler);
static MFImplicitMF MFReadCircle(FILE*,MFErrorHandler);

static int MFCircleProjectToSave(MFNVector,double*,void*,MFErrorHandler);
static int MFCircleProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
static int MFCircleProjectForBB(MFNVector,double*,void*,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

/*! \fn MFImplicitMF MFIMFCreateCircle(double x,double y,double R, MFErrorHandler e);
 *  \brief Creates an implicitly defined manifold for a circle embedded in the plane, centered at (x,y) with
 *         radius R.
 *
 *  \param x The x coordinate of the center of the circle
 *  \param y The y coordinate of the center of the circle
 *  \param R The radius.
 *  \param e An error handler.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateCircle(double x,double y,double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateCircle"};
  MFImplicitMF circle;
  MFNSpace space;
  double *data;

  circle=MFIMFCreateBaseClass(2,1,"Circle",e);

  space=MFCreateNSpace(2,e);
  MFIMFSetSpace(circle,space,e);
  MFFreeNSpace(space,e);

  data=(double*)malloc(3*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFCircleMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",3*sizeof(double));
    MFSetError(e,12,RoutineName,MFCircleMFErrorHandlerMsg,__LINE__,__FILE__);
    free(circle);
    return NULL;
   }
#endif

  data[0]=x;
  data[1]=y;
  data[2]=R;
  MFIMFSetData(circle,data,e);
  MFIMFSetFreeData(circle,MFFreeCircleData,e);
  MFIMFSetProject(circle,MFProjectCircle,e);
  MFIMFSetTangent(circle,MFTangentCircle,e);
  MFIMFSetScale(circle,MFScaleCircle,e);
  MFIMFSetWriteData(circle,MFWriteCircleData,e);
  MFIMFSetProjectForSave(circle,MFCircleProjectToSave,e);
  MFIMFSetProjectForDraw(circle,MFCircleProjectToDraw,e);
  MFIMFSetProjectForBB(circle,MFCircleProjectForBB,e);

  MFIMFSetVectorFactory(circle,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(circle,MFNKMatrixFactory,e);

  return circle;
 }

void MFFreeCircleData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeCircleData"};
  free(d);
  return;
 }

int MFProjectCircle(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectCircle"};
  double t,A,B,C,R,nx,ny;
  double u0x,u0y;
  double v0x,v0y;
  double err;

  double *u, *Phi, *u0;

  u0=MFNV_CStar(vu0,e);
  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  v0x=Phi[0];
  v0y=Phi[1];
  nx= v0y;
  ny=-v0x;
  u0x=u0[0]-((double*)d)[0];
  u0y=u0[1]-((double*)d)[1];
  R=((double*)d)[2];

  A=nx*nx+ny*ny;
  B=nx*u0x+ny*u0y;
  C=u0x*u0x+u0y*u0y-R*R;
  if(A>0.)
   {
    if(B>0)
      t=(-B+sqrt(B*B-A*C))/A;
     else
      t=(-B-sqrt(B*B-A*C))/A;
   }else if(A==0.)
     t=-.5*C/B;
    else return 0;

  u[0]=u0[0]+t*nx;
  u[1]=u0[1]+t*ny;

  *index=0;

  return 1;
 }

int MFTangentCircle(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentCircle"};
  double s;
  double u0x,u0y;
  double v0x,v0y;
  double v1x,v1y;
  double err;
  double *u, *Phi;

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  u0x=(u[0]-((double*)d)[0])/((double*)d)[2];
  u0y=(u[1]-((double*)d)[1])/((double*)d)[2];

  v0x=-u0y;
  v0y= u0x;

  Phi[0]=-u0y;
  Phi[1]= u0x;

  return 1;
 }

double MFScaleCircle(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScaleCircle"};
  return sqrt(.5*MFEpsilon*((double*)d)[2]);
 }

void MFWriteCircleData(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteCircleData"};
  static double *data;

  data=(double*)d;

  fprintf(fid,"%lf %lf %lf\n",data[0],data[1],data[2]);
  fflush(fid);
  return;
 }

MFImplicitMF MFReadCircle(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadCircle"};
  MFImplicitMF circle;
  double x=0.;
  double y=0.;
  double R=0.;

  fscanf(fid,"%lf %lf %lf\n",&x,&y,&R);
  circle=MFIMFCreateCircle(x,y,R,e);

  return circle;
 }

int MFCircleProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  if(x==NULL)return 2;

  x[0]=MFNV_C(u,0,e);
  x[1]=MFNV_C(u,1,e);

  return 0;
 }

int MFCircleProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  if(x==NULL)return 2;

  x[0]=MFNV_C(u,0,e);
  x[1]=MFNV_C(u,1,e);

  return 0;
 }

int MFCircleProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  if(x==NULL)return 2;

  x[0]=MFNV_C(u,0,e);
  x[1]=MFNV_C(u,1,e);

  return 0;
 }

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif
