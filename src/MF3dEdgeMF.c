/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   Sept. 27, 2002  modified MF3dPolygon
 */

static char *id="@(#) $Id: MF3dEdgeMF.c,v 1.4 2011/07/21 17:42:46 mhender Exp $";

static char MFEdgeMFErrorHandlerMsg[256]="";

#include <MFImplicitMF.h>
#include <MFPrint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

static void MFFreeEdgeIn3SpaceData(void*,MFErrorHandler);
static int MFProjectEdgeIn3Space(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentEdgeIn3Space(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static double MFScaleEdgeIn3Space(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static void MFWriteEdgeIn3SpaceData(FILE*,void*,MFErrorHandler);
static MFImplicitMF MFReadEdgeIn3Space(FILE*,MFErrorHandler);
static int MFEdgeIn3SpaceProjectForBB(MFNVector,double*,void*,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

struct MFEdgeIn3SpaceData
 {
  double r;
  double ox,oy,oz;
  double ax,ay,az;
 };

/*! \fn MFImplicitMF MFIMFCreateEdgeIn3SpaceWithRadius(double *o,double *d,double R, MFErrorHandler e);
 *  \brief Creates a Euclidean line which contains the segment between the points l and r (left and right).
 *
 *  \param o An array of length at least 3, with the coordinates of the left end point of the interval.
 *  \param d An array of length at least 3, with the coordinates of the right end point of the interval.
 *  \param R  A radius to use for charts on the manifold.
 *  \param e An error handler.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreateEdgeIn3SpaceWithRadius(double *v0, double *v1,double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreateEdgeIn3SpaceWithRadius"};
  MFImplicitMF edge;
  MFNSpace space;
  struct MFEdgeIn3SpaceData *p;
  double d;
  int verbose=0;
  int i;

#ifdef MFALLOWVERBOSE
  if(verbose)printf("%s\n",RoutineName);
#endif

  edge=MFIMFCreateBaseClass(3,1,"EdgeIn3Space",e);

  space=MFCreateNSpace(3,e);
  MFIMFSetSpace(edge,space,e);
  MFFreeNSpace(space,e);

  p=(struct MFEdgeIn3SpaceData*)malloc(sizeof(struct MFEdgeIn3SpaceData));

#ifndef MFNOSAFETYNET
  if(p==NULL)
   {
    sprintf(MFEdgeMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFEdgeIn3SpaceData));
    MFSetError(e,12,RoutineName,MFEdgeMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  p->r=R;
  p->ox=v0[0];
  p->oy=v0[1];
  p->oz=v0[2];
  p->ax=v1[0]-v0[0];
  p->ay=v1[1]-v0[1];
  p->az=v1[2]-v0[2];
  d=1./sqrt(p->ax*p->ax+p->ay*p->ay+p->az*p->az);
  p->ax=p->ax*d;
  p->ay=p->ay*d;
  p->az=p->az*d;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, data=0x%8.8x, o=(%lf,%lf,%lf), n=(%lf,%lf,%lf)\n",RoutineName,p,p->ox,p->oy,p->oz,p->ax,p->ay,p->az);fflush(stdout);}
#endif

  MFIMFSetData(edge,(void*)p,e);
  MFIMFSetFreeData(edge,MFFreeEdgeIn3SpaceData,e);
  MFIMFSetProject(edge,MFProjectEdgeIn3Space,e);
  MFIMFSetProjectForBB(edge,MFEdgeIn3SpaceProjectForBB,e);
  MFIMFSetProjectForDraw(edge,MFEdgeIn3SpaceProjectForBB,e);
  MFIMFSetTangent(edge,MFTangentEdgeIn3Space,e);
  MFIMFSetScale(edge,MFScaleEdgeIn3Space,e);
  MFIMFSetR(edge,R,e);
  MFIMFSetWriteData(edge,MFWriteEdgeIn3SpaceData,e);

  MFIMFSetVectorFactory(edge,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(edge,MFNKMatrixFactory,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return edge;
 }

void MFFreeEdgeIn3SpaceData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeEdgeIn3SpaceData"};
  struct MFEdgeIn3SpaceData *p;

  p=(struct MFEdgeIn3SpaceData*)d;
  free(p);
  return;
 }

int MFProjectEdgeIn3Space(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectEdgeIn3Space"};
  static int i;
  double a,b,nm;
  struct MFEdgeIn3SpaceData *p;
  double *u0,*Phi,*u;
  double s;
  int verbose=0;

  u0=MFNV_CStar(vu0,e);
  u=MFNV_CStar(vu,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  p=(struct MFEdgeIn3SpaceData*)d;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s, data=0x%8.8x, o=(%lf,%lf,%lf), n=(%lf,%lf,%lf)\n",RoutineName,p,p->ox,p->oy,p->oz,p->ax,p->ay,p->az);fflush(stdout);}
#endif

  s=(u0[0]-p->ox)*p->ax+(u0[1]-p->oy)*p->ay+(u0[2]-p->oz)*p->az;

#ifdef MFALLOWVERBOSE
  if(verbose){ printf("    u=");MFPrintNVector(stdout,vu,e);fflush(stdout);
  printf("    u=(%lf,%lf,%lf), s=%lf\n",u0[0],u0[1],u0[2],s);fflush(stdout);}
#endif

  u[0]=p->ox+s*p->ax;
  u[1]=p->oy+s*p->ay;
  u[2]=p->oz+s*p->az;

  *index=0;
  return 1;
 }

int MFTangentEdgeIn3Space(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentEdgeIn3Space"};
  static int i;
  struct MFEdgeIn3SpaceData *p;
  double *Phi,*u;

  Phi=MFNKM_CStar(mPhi,e);

  p=(struct MFEdgeIn3SpaceData*)d;
  Phi[0]=p->ax;
  Phi[1]=p->ay;
  Phi[2]=p->az;

  return 0;
 }

double MFScaleEdgeIn3Space(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScaleEdgeIn3Space"};
  struct MFEdgeIn3SpaceData *p;

  p=(struct MFEdgeIn3SpaceData*)d;
  return p->r;
 }

void MFWriteEdgeIn3SpaceData(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteEdgeIn3SpaceData"};
  struct MFEdgeIn3SpaceData *data;

  data=(struct MFEdgeIn3SpaceData*)d;
  fprintf(fid,"%lf %lf %lf %lf %lf %lf %lf\n",data->r,data->ox,data->oy,data->oz,data->ax,data->ay,data->az);
  fflush(fid);
  return;
 }

MFImplicitMF MFReadEdgeIn3Space(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadEdgeIn3Space"};
  MFImplicitMF edge;
  struct MFEdgeIn3SpaceData *p;
  double d;
  int verbose=0;
  int i;
  double r=0.;
  double ox=0.;
  double oy=0.;
  double oz=0.;
  double ax=0.;
  double ay=0.;
  double az=0.;
  double v0[3];
  double v1[3];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s\n",RoutineName);
    fflush(stdout);
   }
#endif

  fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf\n",&r,&ox,&oy,&oz,&ax,&ay,&az);

  v0[0]=ox;
  v0[1]=oy;
  v0[2]=oz;
  v1[0]=ox+ax;
  v1[1]=oy+ay;
  v1[2]=oz+az;

  edge=MFIMFCreateEdgeIn3SpaceWithRadius(v0,v1,r,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return edge;
 }

int MFEdgeIn3SpaceProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
 {
  int i;

  if(x==NULL)return 3;

  for(i=0;i<3;i++)x[i]=MFNV_C(u,i,e);

  return 0;
 }

/*! @} */

/*! @} */

#ifdef __cplusplus
}
#endif
