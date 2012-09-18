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

static char *id="@(#) $Id: MFTorus.c,v 1.3 2007/02/13 01:22:34 mhender Exp $";

static char MFTorusMFErrorHandlerMsg[256]="";

extern double MFEpsilon;

#include <MFImplicitMF.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

static void MFFreeTorusData(void*,MFErrorHandler);
static int MFProjectTorus(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentTorus(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static double MFScaleTorus(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static void MFWriteTorusData(FILE*,void*,MFErrorHandler);
static MFImplicitMF MFReadTorus(FILE*,MFErrorHandler);

static int MFTorusProjectToSave(MFNVector,double*,void*,MFErrorHandler);
static int MFTorusProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
static int MFTorusProjectForBB(MFNVector,double*,void*,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

#define MFTWOPI 6.2831853071795862320

/*! \fn MFImplicitMF MFIMFCreateTorus(double x,double y,double z,double RI,double RO);
 *  \brief Creates an implicitly defined manifold for a 2-torus embedded in 3-space, centered at (x,y,z) with
 *         inner radius RI and outer radius RO. If RO were 0. (don't! the eqs are singular then), the torus
 *         would be a circle of radius RI in the xy-plane. You are free to set RI<RO, but it won't look like a torus.
 *
 *  \param x The x coordinate of the center of the torus
 *  \param y The y coordinate of the center of the torus
 *  \param z The z coordinate of the center of the torus
 *  \param RI The inner radius.
 *  \param RO The outer radius.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateTorus(double x,double y,double z,double RI, double RO, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateTorus"};
  MFImplicitMF torus;
  MFNSpace space;
  double *data;

  torus=MFIMFCreateBaseClass(3,2,"Torus",e);

  space=MFCreateNSpace(3,e);
  MFIMFSetSpace(torus,space,e);
  MFFreeNSpace(space,e);

  data=(double*)malloc(5*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFTorusMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",5*sizeof(double));
    MFSetError(e,12,RoutineName,MFTorusMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(torus);
    return NULL;
   }
#endif

  data[0]=x;
  data[1]=y;
  data[2]=z;
  data[3]=RI;
  data[4]=RO;
  MFIMFSetData(torus,data,e);
  MFIMFSetFreeData(torus,MFFreeTorusData,e);
  MFIMFSetProject(torus,MFProjectTorus,e);
  MFIMFSetTangent(torus,MFTangentTorus,e);
  MFIMFSetScale(torus,MFScaleTorus,e);
  MFIMFSetWriteData(torus,MFWriteTorusData,e);

  MFIMFSetProjectForSave(torus,MFTorusProjectToSave,e);
  MFIMFSetProjectForDraw(torus,MFTorusProjectToDraw,e);
  MFIMFSetProjectForBB(torus,MFTorusProjectForBB,e);

  MFIMFSetVectorFactory(torus,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(torus,MFNKMatrixFactory,e);

  return torus;
 }

void MFFreeTorusData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeTorusData"};
  free(d);
  return;
 }

int MFProjectTorus(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectTorus"};
  double R,R1,R2;
  double nx,ny,nz;
  double x0,y0,z0;
  double x,y,z;
  double l0x,l0y,l0z;
  double l1x,l1y,l1z;
  double s,t,dt,F,dF;
  double disc,ddisc;
  double err;
  double *u0,*Phi,*u;

  u0=MFNV_CStar(vu0,e);
  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

/* (sqrt(x^2+y^2)-R1)^2+z^2 -R2^2 = 0; */

  l0x=Phi[0];
  l0y=Phi[1];
  l0z=Phi[2];
  l1x=Phi[0+3];
  l1y=Phi[1+3];
  l1z=Phi[2+3];
  nx=l0y*l1z-l0z*l1y;
  ny=l0z*l1x-l0x*l1z;
  nz=l0x*l1y-l0y*l1x;
  x0=u0[0]-((double*)d)[0];
  y0=u0[1]-((double*)d)[1];
  z0=u0[2]-((double*)d)[2];
  R1=((double*)d)[3];
  R2=((double*)d)[4];

  n=0;
  t=0.;
  while(TRUE)
   {
    disc=sqrt(pow(x0+t*nx,2)+pow(y0+t*ny,2));
    F=pow(disc-R1,2)+pow(z0+t*nz,2)-pow(R2,2);
    if(fabs(F)<1.e-7)
     {
      x=x0+t*nx;
      y=y0+t*ny;
      z=z0+t*nz;

      u[0]=((double*)d)[0]+x0+t*nx;
      u[1]=((double*)d)[1]+y0+t*ny;
      u[2]=((double*)d)[2]+z0+t*nz;
      if(FALSE)printf("\n");fflush(stdout);

      *index=0;
      return 1;
     }
    ddisc=(nx*(x0+t*nx)+ny*(y0+t*ny))/disc;
    dF=2*(disc-R1)*ddisc+2*(z0+t*nz)*nz;
    dt=-F/dF;
    if(dt>.01)dt=.01;
    if(dt<-.01)dt=-.01;

    if(n>20)return 0;

    t=t+dt;
    n++;
   }
  *index=0;
  return 1;
 }

int MFTangentTorus(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentTorus"};
  double D,s;
  double x,y,z;
  double R1,R2;
  double *u,*Phi;

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

/* (sqrt(x^2+y^2)-R1)^2+z^2 -R2^2 = 0; */

  x=u[0]-((double*)d)[0];
  y=u[1]-((double*)d)[1];
  z=u[2]-((double*)d)[2];
  R1=((double*)d)[3];
  R2=((double*)d)[4];

/*printf("Torus Tangent (%lf,%lf,%lf)\n",x,y,z);fflush(stdout);*/

/* (1-R1/sqrt(x^2+y^2)) (2x x' + 2y y') +2z z' = 0 */

  D=1-R1/sqrt(x*x+y*y);
  Phi[0]= y*D;
  Phi[1]=-x*D;
  Phi[2]=0.;
  s=1./sqrt(Phi[0]*Phi[0]+Phi[1]*Phi[1]);
  Phi[0]=Phi[0]*s;
  Phi[1]=Phi[1]*s;

  if(fabs(y)>.1)
   {
    Phi[4]=-y*z;
    Phi[5]=D*(x*x+y*y);
    Phi[3]=x/y*Phi[4];
   }else{
    Phi[3]=-x*z;
    Phi[5]=D*(x*x+y*y);
    Phi[4]=y/x*Phi[3];
   }

  s=1./sqrt(Phi[3]*Phi[3]+Phi[4]*Phi[4]+Phi[5]*Phi[5]);
  Phi[3]=Phi[3]*s;
  Phi[4]=Phi[4]*s;
  Phi[5]=Phi[5]*s;

  return 1;
 }

double MFScaleTorus(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScaleTorus"};

/* (sqrt(x^2+y^2)-R1)^2+z^2 -R2^2 = 0; */
/* (1-R1/sqrt(x^2+y^2)) (x*dx + y*dy) +z dz = 0 */
/* (1-R1/sqrt(x^2+y^2)) (x*dDx +y*dDy)+z*dDz+(1-R1/sqrt(x^2+y^2)) (dxDx+dyDy)+dzDz+R1/sqrt(x^2+y^2)(xDx+yDy)(xdx+ydy) +dzDz = 0 */

  double dg,ddg;
  double sqrtR;
  double x,y,z,R1,R2;
  double nx,ny,nz;
  double l,r;
  double dx,dy,dz;
  double Dx,Dy,Dz;
  double d00,d01,d11;
  double *u,*Phi;

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  x=u[0]-((double*)d)[0];
  y=u[1]-((double*)d)[1];
  z=u[2]-((double*)d)[2];
  R1=((double*)d)[3];
  R2=((double*)d)[4];
  sqrtR=sqrt(x*x+y*y);

  nx=Phi[1]*Phi[5]-Phi[2]*Phi[4];
  ny=Phi[2]*Phi[3]-Phi[0]*Phi[5];
  nz=Phi[0]*Phi[4]-Phi[1]*Phi[3];

  dx=Phi[0];
  dy=Phi[1];
  dz=Phi[2];
  Dx=Phi[0];
  Dy=Phi[1];
  Dz=Phi[2];
  ddg=(1-R1/sqrtR)*(dx*Dx+dy*Dy)+dz*Dz
      +R1/sqrtR/(x*x+y*y)*(x*Dx+y*Dy)*(x*dx+y*dy)+dz*Dz;
  dg=(1-R1/sqrtR)*(x*nx+y*ny)+z*nz;
  d00=-ddg/dg;

  dx=Phi[0];
  dy=Phi[1];
  dz=Phi[2];
  Dx=Phi[3];
  Dy=Phi[4];
  Dz=Phi[5];
  ddg=(1-R1/sqrtR)*(dx*Dx+dy*Dy)+dz*Dz
      +R1/sqrtR/(x*x+y*y)*(x*Dx+y*Dy)*(x*dx+y*dy)+dz*Dz;
  dg=(1-R1/sqrtR)*(x*nx+y*ny)+z*nz;
  d01=-ddg/dg;

  dx=Phi[3];
  dy=Phi[4];
  dz=Phi[5];
  Dx=Phi[3];
  Dy=Phi[4];
  Dz=Phi[5];
  ddg=(1-R1/sqrtR)*(dx*Dx+dy*Dy)+dz*Dz
      +R1/sqrtR/(x*x+y*y)*(x*Dx+y*Dy)*(x*dx+y*dy)+dz*Dz;
  dg=(1-R1/sqrtR)*(x*nx+y*ny)+z*nz;
  d11=-ddg/dg;

/* Ev sat (d00-l)*(d11-l)+d01*d01=0 */

  if(d00+d11>0)
    l=.5*(d00+d11)+.5*sqrt((d00-d11)*(d00-d11)+4*d01*d01);
   else
    l=-.5*(d00+d11)+.5*sqrt((d00-d11)*(d00-d11)+4*d01*d01);

  r=sqrt(2*MFEpsilon/l);

/*
  printf("  MFScaleTorus, (%lf,%lf,%lf), R1=%lf, R2=%lf\n",x,y,z,R1,R2);fflush(stdout);
  printf("                c0=%lf c1=%lf, lev=%lf, r=%lf\n",d00,d11,l,r);fflush(stdout);
*/

  return r;
 }

void MFWriteTorusData(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteTorusData"};
  double *data;

  data=(double*)d;

  fprintf(fid,"%lf %lf %lf %lf %lf\n",data[0],data[1],data[2],data[3],data[4]);
  fflush(fid);
  return;
 }

MFImplicitMF MFReadTorus(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadTorus"};
  MFImplicitMF torus;
  double x=0.;
  double y=0.;
  double z=0.;
  double RI=0.;
  double RO=0.;

  fscanf(fid,"%lf %lf %lf %lf %lf\n",&x,&y,&z,&RI,&RO);
  torus=MFIMFCreateTorus(x,y,z,RI,RO,e);

  return torus;
 }

int MFTorusProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  if(x==NULL)return 3;

  x[0]=MFNV_C(u,0,e);
  x[1]=MFNV_C(u,1,e);
  x[2]=MFNV_C(u,2,e);

  return 0;
 }

int MFTorusProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  if(x==NULL)return 3;

  x[0]=MFNV_C(u,0,e);
  x[1]=MFNV_C(u,1,e);
  x[2]=MFNV_C(u,2,e);

  return 0;
 }

int MFTorusProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
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
