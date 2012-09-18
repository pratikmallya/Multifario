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

static char *id="@(#) $Id: MFSphere.c,v 1.3 2007/02/13 01:22:34 mhender Exp $";

static char MFSphereMFErrorHandlerMsg[256]="";

#include <MFImplicitMF.h>
#include <MFPrint.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

static void MFFreeSphereData(void*,MFErrorHandler);
static int MFProjectSphere(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentSphere(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static double MFScaleSphere(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static void MFWriteSphereData(FILE*,void*,MFErrorHandler);
static MFImplicitMF MFReadSphere(FILE*,MFErrorHandler);

static int MFSphereProjectToSave(MFNVector,double*,void*,MFErrorHandler);
static int MFSphereProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
static int MFSphereProjectForBB(MFNVector,double*,void*,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

/*! \fn MFImplicitMF MFIMFCreateSphere(double x,double y,double z,double R, MFErrorHandler e);
 *  \brief Creates an implicitly defined manifold for a 2-sphere embedded in 3-space, centered at (x,y,z) with
 *         radius R.
 *
 *  \param x The x coordinate of the center of the sphere
 *  \param y The y coordinate of the center of the sphere
 *  \param z The z coordinate of the center of the sphere
 *  \param R The radius.
 *  \param e An error handler.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateSphere(double x,double y,double z,double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateSphere"};
  MFImplicitMF sphere;
  MFNSpace space;
  double *data;

  sphere=MFIMFCreateBaseClass(3,2,"Sphere",e);

  space=MFCreateNSpace(3,e);
  MFIMFSetSpace(sphere,space,e);
  MFFreeNSpace(space,e);

  data=(double*)malloc(4*sizeof(double));

#ifndef MFNPOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFSphereMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",4*sizeof(double));
    MFSetError(e,12,RoutineName,MFSphereMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(sphere);
    return NULL;
   }
#endif

  data[0]=x;
  data[1]=y;
  data[2]=z;
  data[3]=R;
  MFIMFSetData(sphere,(void*)data,e);
  MFIMFSetFreeData(sphere,MFFreeSphereData,e);
  MFIMFSetProject(sphere,MFProjectSphere,e);
  MFIMFSetTangent(sphere,MFTangentSphere,e);
  MFIMFSetR(sphere,-1.,e);
  MFIMFSetScale(sphere,MFScaleSphere,e);
  MFIMFSetWriteData(sphere,MFWriteSphereData,e);
  MFIMFSetProjectForSave(sphere,MFSphereProjectToSave,e);
  MFIMFSetProjectForDraw(sphere,MFSphereProjectToDraw,e);
  MFIMFSetProjectForBB(sphere,MFSphereProjectForBB,e);

  MFIMFSetVectorFactory(sphere,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(sphere,MFNKMatrixFactory,e);

  return sphere;
 }

void MFFreeSphereData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeSphereData"};
  free(d);
  return;
 }

int MFProjectSphere(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectSphere"};
  double t,A,B,C,R,nx,ny,nz;
  double u0x,u0y,u0z;
  double v0x,v0y,v0z;
  double v1x,v1y,v1z;
  double err;
  double *u,*u0,*Phi;
  int verbose=0;

  u0=MFNV_CStar(vu0,e);
  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  v0x=Phi[0];
  v0y=Phi[1];
  v0z=Phi[2];
  v1x=Phi[0+3];
  v1y=Phi[1+3];
  v1z=Phi[2+3];
  nx=v0y*v1z-v0z*v1y;
  ny=v0z*v1x-v0x*v1z;
  nz=v0x*v1y-v0y*v1x;
  u0x=u0[0]-((double*)d)[0];
  u0y=u0[1]-((double*)d)[1];
  u0z=u0[2]-((double*)d)[2];
  R=((double*)d)[3];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("In Project sphere\n");
    printf(" R %lf o=(%lf,%lf,%lf)\n",R,((double*)d)[0],((double*)d)[1],((double*)d)[2]);fflush(stdout);
    printf(" u0 (%lf,%lf,%lf) : ",u0x,u0y,u0z);MFPrintNVector(stdout,vu0,e);printf("\n");
    printf(" t1 (%lf,%lf,%lf) t2 (%lf,%lf,%lf) :\n",v0x,v0y,v0z,v1x,v1y,v1z);MFPrintNKMatrix(stdout,mPhi,e);printf("\n");
   }
#endif

  A=nx*nx+ny*ny+nz*nz;
  B=nx*u0x+ny*u0y+nz*u0z;
  C=u0x*u0x+u0y*u0y+u0z*u0z-R*R;

#ifdef MFALLOWVERBOSE
  if(verbose)printf(" quadratic  (%lf)t^2+(%lf)t+(%lf)=0\n",A,B,C);
#endif

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
  u[2]=u0[2]+t*nz;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf(" n  (%lf,%lf,%lf)\n",nz,ny,nz);
    printf(" u  (%lf,%lf,%lf) : ",u[0],u[1],u[2]);MFPrintNVector(stdout,vu,e);printf("\n");

/* Test */

    err=fabs((u[0]-u0[0])*Phi[3]+(u[1]-u0[1])*Phi[4]+(u[2]-u0[2])*Phi[5]);
    err+=fabs((u[0]-u0[0])*Phi[0]+(u[1]-u0[1])*Phi[1]+(u[2]-u0[2])*Phi[2]);
    err+=fabs(sqrt( pow(u[0]-((double*)d)[0],2)+pow(u[1]-((double*)d)[1],2)+pow(u[2]-((double*)d)[2],2))-((double*)d)[3]);
    printf("Error %le\n",err);
   }
#endif

  *index=0;
  return 1;
 }

int MFTangentSphere(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentSphere"};
  double s;
  double u0x,u0y,u0z;
  double v0x,v0y,v0z;
  double v1x,v1y,v1z;
  double err;
  double *u,*Phi;
  int verbose=0;

#ifdef MFALLOWVERBOSE
  if(verbose)printf("In %s\n",RoutineName);
#endif

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  u0x=(u[0]-((double*)d)[0])/((double*)d)[3];
  u0y=(u[1]-((double*)d)[1])/((double*)d)[3];
  u0z=(u[2]-((double*)d)[2])/((double*)d)[3];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf(" R %lf o=(%lf,%lf,%lf)\n",((double*)d)[3],((double*)d)[0],((double*)d)[1],((double*)d)[2]);fflush(stdout);
    printf(" u0 (%lf,%lf,%lf) : ",u0x,u0y,u0z);MFPrintNVector(stdout,vu,e);printf("\n");
   }
#endif

  if(u0x!=1.)
   {
    v0x=1.-u0x*u0x;
    v0y=  -u0x*u0y;
    v0z=  -u0x*u0z;
   }else if(u0y!=1.)
   {
    v0x=  -u0y*u0x;
    v0y=1.-u0y*u0y;
    v0z=  -u0y*u0z;
   }else{
    v0x=  -u0z*u0x;
    v0y=  -u0z*u0y;
    v0z=1.-u0z*u0z;
   }
  s=1./sqrt(v0x*v0x+v0y*v0y+v0z*v0z);
  v0x=s*v0x;
  v0y=s*v0y;
  v0z=s*v0z;

  v1x=-v0y*u0z+v0z*u0y;
  v1y=-v0z*u0x+v0x*u0z;
  v1z=-v0x*u0y+v0y*u0x;

  Phi[0]=v0x;
  Phi[1]=v0y;
  Phi[2]=v0z;
  Phi[0+3]=v1x;
  Phi[1+3]=v1y;
  Phi[2+3]=v1z;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf(" t1 (%lf,%lf,%lf) t2 (%lf,%lf,%lf) :\n",v0x,v0y,v0z,v1x,v1y,v1z);MFPrintNKMatrix(stdout,mPhi,e);printf("\n");

/* Test */
    err=fabs(Phi[0]*Phi[0]+Phi[1]*Phi[1]+Phi[2]*Phi[2]-1.);
    err+=fabs(Phi[3]*Phi[0]+Phi[4]*Phi[1]+Phi[5]*Phi[2]);
    err+=fabs(Phi[3]*Phi[3]+Phi[4]*Phi[4]+Phi[5]*Phi[5]-1.);
    err+=fabs(Phi[0]*u0x+Phi[1]*u0y+Phi[2]*u0z);
    err+=fabs(Phi[3]*u0x+Phi[4]*u0y+Phi[5]*u0z);

    printf("Error %le\n",err);
   }
#endif

  return 1;
 }

extern double MFEpsilon;

double MFScaleSphere(int n,int k,MFNVector u,MFNKMatrix Phi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScaleSphere"};
/*return sqrt(.5*((double*)d)[3]/.1);*/

  return .95*sqrt(2*MFEpsilon)*((double*)d)[3];
 }

void MFWriteSphereData(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteSphereData"};
  double *data;

  data=(double*)d;

  fprintf(fid,"%lf %lf %lf %lf\n",data[0],data[1],data[2],data[3]);
  fflush(fid);
  return;
 }

MFImplicitMF MFReadSphere(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadSphere"};
  MFImplicitMF sphere;
  double x=0.;
  double y=0.;
  double z=0.;
  double R=0.;

  fscanf(fid,"%lf %lf %lf %lf\n",&x,&y,&z,&R);
  sphere=MFIMFCreateSphere(x,y,z,R,e);

  return sphere;
 }

int MFSphereProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  if(x==NULL)return 3;

  x[0]=MFNV_C(u,0,e);
  x[1]=MFNV_C(u,1,e);
  x[2]=MFNV_C(u,2,e);

  return 0;
 }

int MFSphereProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  if(x==NULL)return 3;

  x[0]=MFNV_C(u,0,e);
  x[1]=MFNV_C(u,1,e);
  x[2]=MFNV_C(u,2,e);

  return 0;
 }

int MFSphereProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  if(x==NULL)return 3;

  x[0]=MFNV_C(u,0,e);
  x[1]=MFNV_C(u,1,e);
  x[2]=MFNV_C(u,2,e);

  return 0;
 }

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif
