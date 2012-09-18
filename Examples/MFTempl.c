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

static char *id="@(#) $Id: MFTempl.c,v 1.3 2007/02/13 01:22:34 mhender Exp $";

static char MFTemplMFErrorHandlerMsg[256]="";

extern double MFEpsilon;

#include <MFImplicitMF.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

static void MFFreeTemplData(void*,MFErrorHandler);
static int MFProjectTempl(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentTempl(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static double MFScaleTempl(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static void MFWriteTemplData(FILE*,void*,MFErrorHandler);
static MFImplicitMF MFReadTempl(FILE*,MFErrorHandler);

static int MFTemplProjectToSave(MFNVector,double*,void*,MFErrorHandler);
static int MFTemplProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
static int MFTemplProjectForBB(MFNVector,double*,void*,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

#define MFTWOPI 6.2831853071795862320

/*! \fn MFImplicitMF MFIMFCreateTempl(double x,double y,double z,double RI,double RO);
 *  \brief Creates an inv manifold for a 1-d stiff system embedded in 2-space
 *  the system is x_dot = -x, y_dot = -stiffness(y-x^2)
 *  \param stiffness only
 *  \returns The manifold.
 */
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


MFImplicitMF MFIMFCreateTempl(double stiff, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateTempl"};
  MFImplicitMF templ;
  MFNSpace space;
  double *data;

  templ=MFIMFCreateBaseClass(3,2,"Templ",e);

  space=MFCreateNSpace(3,e);
  MFIMFSetSpace(templ,space,e);
  MFFreeNSpace(space,e);

  data=(double*)malloc(1*sizeof(double)); /*done*/

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFTemplMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",5*sizeof(double));
    MFSetError(e,12,RoutineName,MFTemplMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(templ);
    return NULL;
   }
#endif

  data[0]=stiff;
  MFIMFSetData(templ,data,e);
  MFIMFSetFreeData(templ,MFFreeTemplData,e);
  MFIMFSetProject(templ,MFProjectTempl,e);
  MFIMFSetTangent(templ,MFTangentTempl,e);
  MFIMFSetScale(templ,MFScaleTempl,e);
  MFIMFSetWriteData(templ,MFWriteTemplData,e);

  MFIMFSetProjectForSave(templ,MFTemplProjectToSave,e);
  MFIMFSetProjectForDraw(templ,MFTemplProjectToDraw,e);
  MFIMFSetProjectForBB(templ,MFTemplProjectForBB,e);

  MFIMFSetVectorFactory(templ,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(templ,MFNKMatrixFactory,e);

  return templ;
 }

void MFFreeTemplData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeTemplData"};
  free(d);
  return;
 }

int MFProjectTempl(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectTempl"};
  double nx,ny;
  double x0,y0,z0,zz0;
  double x,y,z,zz;
  double s,t,dt,F,dF;
  double disc,ddisc;
  double err;
  double *u0,*Phi,*u;
  double stifff;
  MFNVector tang;

  printf("Num col why doesnt this get outputted? %lf \n",MFNKMatrixK(mPhi,e));fflush(stdout);
  printf("Num row why doesnt this get outputted? %lf \n",MFNKMatrixN(mPhi,e));fflush(stdout);

  tang=MFMColumn(mPhi,0,e);

  nx=-MFNV_C(tang,1,e);
  ny= MFNV_C(tang,0,e);
  
  x0=MFNV_C(vu0,0,e);
  y0=MFNV_C(vu0,1,e);
  z0=MFNV_C(vu0,2,e);
  stifff = ((double*)d)[0];

  n=0;
  t=0.;
  while(TRUE)
   {
    F=(y0+t*ny)-stifff/(stifff-2.)*pow(x0+t*nx,2);
    if(fabs(F)<1.e-7)
     {
      x=x0+t*nx;
      y=y0+t*ny;
      z=z0;

      MFNVSetC(vu,0,x,e);
      MFNVSetC(vu,1,y,e);
      MFNVSetC(vu,2,z,e);

      *index=0;
      return 1;
     }
    dF=ny-stifff/(stifff-2.)*2*(x0+t*nx)*nx;
    dt=-F/dF;
    
    if(n>20)return 0;

    t=t+dt;
    n++;
   }
  *index=0;
  MFFreeNVector(tang,e);
  return 1;
 }

int MFTangentTempl(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentTempl"};
  double D,s;
  double x,y,z,zz;
  double R1,R2;
  double *u,*Phi;
  double stifff;
  double phi0,phi1;
  MFNVector tang1,tang2;


  x=MFNV_C(vu,0,e);
  y=MFNV_C(vu,1,e);
  z=MFNV_C(vu,2,e);
  stifff = ((double*)d)[0];

  phi0=1;
  phi1=2*x*(stifff)/(stifff-2);
  s=1./sqrt(phi0*phi0+phi1*phi1);
  phi0=phi0*s;
  phi1=phi1*s;

  tang1=MFMColumn(mPhi,0,e);
  tang2=MFMColumn(mPhi,1,e);

  MFNVSetC(tang1,0,phi0,e);
  MFNVSetC(tang1,1,phi1,e);
  MFNVSetC(tang1,2,0,e);
  MFNVSetC(tang2,0,0,e);
  MFNVSetC(tang2,1,0,e);
  MFNVSetC(tang2,2,1,e);

  MFMSetColumn(mPhi,0,tang1,e);     /* Was 1 */
  MFMSetColumn(mPhi,1,tang2,e);     /* Was 2 */

  MFFreeNVector(tang1,e);
  MFFreeNVector(tang2,e);
  return 1;
 }

double MFScaleTempl(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScaleTempl"};

  double dg,ddg;
  double sqrtR;
  double x,y,z,R1,R2;
  double nx,ny,nz;
  double l,r;
  double dx,dy,dz;
  double Dx,Dy,Dz;
  double d00,d01,d11;
  double *u,*Phi;
  double stifff;

  /*Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);*/

  /*x=u[0];
  y=u[1];*/
  x=MFNV_C(vu,0,e);
  y=MFNV_C(vu,1,e);
  stifff=((double*)d)[0];

  l=2*stifff/(stifff-2.)/pow(sqrt(1+pow(2*stifff/(stifff-2.)*x,2)),3);
  r=sqrt(2*MFEpsilon/l);

  return r;
 }

void MFWriteTemplData(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteTemplData"};
  double *data;

  data=(double*)d;

  /*fprintf(fid,"%lf %lf %lf %lf %lf\n",data[0],data[1],data[2],data[3],data[4]);*/
  fprintf(fid,"%lf\n",data[0]);
  fflush(fid);
  return;
 }

MFImplicitMF MFReadTempl(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadTempl"};
  MFImplicitMF templ;
  /*double x=0.;
  double y=0.;
  double z=0.;
  double RI=0.;
  double RO=0.;*/
  double stifff=0.;

  /*fscanf(fid,"%lf %lf %lf %lf %lf\n",&x,&y,&z,&RI,&RO);*/
  fscanf(fid,"%lf\n",&stifff);
  templ=MFIMFCreateTempl(stifff,e);

  return templ;
 }

int MFTemplProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  if(x==NULL)return 3;

  x[0]=MFNV_C(u,0,e);
  x[1]=MFNV_C(u,1,e);
  x[2]=MFNV_C(u,2,e);
  /*added this line*/
  /*x[3]=MFNV_C(u,3,e);*/

  return 0;
 }

int MFTemplProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  if(x==NULL)return 3;

  x[0]=MFNV_C(u,0,e);
  x[1]=MFNV_C(u,1,e);
  x[2]=MFNV_C(u,2,e);
  /*added this line*/
  /*x[3]=MFNV_C(u,3,e);*/

  return 0;
 }

int MFTemplProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  if(x==NULL)return 3;

  x[0]=MFNV_C(u,0,e);
  x[1]=MFNV_C(u,1,e);
  x[2]=MFNV_C(u,2,e);
  /*added this line*/
  /*x[3]=MFNV_C(u,3,e);*/

  return 0;
 }

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif
