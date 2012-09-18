/* 
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   March 1, 1999
 */

static char *id="@(#) $Id: MFPlane.c,v 1.2 2007/02/13 01:22:34 mhender Exp $";

static char MFPlaneMFErrorHandlerMsg[256]="";

#include <MFImplicitMF.h>
#include <MFErrorHandler.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

static int MFProjectPlane(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentPlane(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static double MFScalePlane(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static MFImplicitMF MFReadPlane(FILE*,MFErrorHandler);
static int MFPlaneProjectToSave(MFNVector,double*,void*,MFErrorHandler);
static int MFPlaneProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
static int MFPlaneProjectForBB(MFNVector,double*,void*,MFErrorHandler);
MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

/*! \fn MFImplicitMF MFIMFCreatePlane(void);
 *  \brief Creates the plane. (not a very exciting manifold I'll grant you).
 *
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreatePlane(MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreatePlane"};
  MFImplicitMF plane;
  MFNSpace space;

  plane=MFIMFCreateBaseClass(2,2,"Plane",e);

  space=MFCreateNSpace(2,e);
  MFIMFSetSpace(plane,space,e);
  MFFreeNSpace(space,e);

  MFIMFSetProject(plane,MFProjectPlane,e);
  MFIMFSetTangent(plane,MFTangentPlane,e);
  MFIMFSetScale(plane,MFScalePlane,e);
  MFIMFSetProjectForSave(plane,MFPlaneProjectToSave,e);
  MFIMFSetProjectForDraw(plane,MFPlaneProjectToDraw,e);
  MFIMFSetProjectForBB(plane,MFPlaneProjectForBB,e);
  MFIMFSetVectorFactory(plane,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(plane,MFNKMatrixFactory,e);

  return plane;
 }

int MFProjectPlane(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectPlane"};
  double *u,*u0,*Phi;

  u0=MFNV_CStar(vu0,e);
  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  u[0]=u0[0];
  u[1]=u0[1];
  *index=0;
  return 1;
 }

int MFTangentPlane(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentPlane"};
  double *u,*Phi;

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  Phi[0]=1.;
  Phi[1]=0.;
  Phi[2]=0.;
  Phi[3]=1.;

  return 1;
 }

double MFScalePlane(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScalePlane"};
  return .1;
 }

MFImplicitMF MFReadPlane(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadPlane"};
  MFImplicitMF plane;

  plane=MFIMFCreatePlane(e);

  return plane;
 }

int MFPlaneProjectToSave(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  if(x==NULL)return 2;

  x[0]=MFNV_C(u,0,e);
  x[1]=MFNV_C(u,1,e);

  return 0;
 }

int MFPlaneProjectToDraw(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  if(x==NULL)return 2;

  x[0]=MFNV_C(u,0,e);
  x[1]=MFNV_C(u,1,e);

  return 0;
 }

int MFPlaneProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
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
