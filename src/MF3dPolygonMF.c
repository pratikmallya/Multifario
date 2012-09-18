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

static char *id="@(#) $Id: MF3dPolygonMF.c,v 1.4 2011/07/21 17:42:46 mhender Exp $";

static char MFPolygonMFErrorHandlerMsg[256]="";

#include <MFImplicitMF.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

static void MFFreePolygonIn3SpaceData(void*,MFErrorHandler);
static int MFProjectPolygonIn3Space(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentPolygonIn3Space(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static double MFScalePolygonIn3Space(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static void MFWritePolygonIn3SpaceData(FILE*,void*,MFErrorHandler);
static MFImplicitMF MFReadPolygonIn3Space(FILE*,MFErrorHandler);
static int MFPolygonIn3SpaceProjectForBB(MFNVector,double*,void*,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

struct MFPolygonIn3SpaceData
 {
  double r;
  double ox,oy,oz;
  double ax,ay,az;
  double bx,by,bz;
  double nx,ny,nz;
 };

/*! \fn MFImplicitMF MFIMFCreatePolygonIn3SpaceWithRadius(int nv,double *v,double R, MFErrorHandler e);
 *  \brief Creates a Euclidean plane which contains the given polygon (which is assumed to be flat).
 *
 *  \param nv The number of vertices.
 *  \param v  The vertices stored as (v[0+3*iv],v[1+3*iv],v[2+3*iv]).
 *  \param R  A radius to use for charts on the manifold.
 *  \param e  An error handler.
 *  \returns An implicitly defined manifold.
 */
MFImplicitMF MFIMFCreatePolygonIn3SpaceWithRadius(int n,double *v,double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFIMFCreatePolygonIn3SpaceWithRadius"};
  MFImplicitMF polygon;
  MFNSpace space;
  struct MFPolygonIn3SpaceData *p;
  double d;
  int verbose=0;
  int i;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s\n",RoutineName);
    printf("%d\n",n);
    for(i=0;i<n;i++)
     printf("   %d (%lf,%lf,%lf)\n",i,v[3*i],v[3*i+1],v[3*i+2]);
    fflush(stdout);
   }
#endif

  polygon=MFIMFCreateBaseClass(3,2,"PolygonIn3Space",e);

  space=MFCreateNSpace(3,e);
  MFIMFSetSpace(polygon,space,e);
  MFFreeNSpace(space,e);

  p=(struct MFPolygonIn3SpaceData*)malloc(sizeof(struct MFPolygonIn3SpaceData));

#ifndef MFNOSAFETYNET
  if(p==NULL)
   {
    sprintf(MFPolygonMFErrorHandlerMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFPolygonIn3SpaceData));
    MFSetError(e,12,RoutineName,MFPolygonMFErrorHandlerMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  p->r=R;
  p->ox=v[0];
  p->oy=v[1];
  p->oz=v[2];

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  o=(%lf,%lf,%lf)\n",p->ox,p->oy,p->oz);fflush(stdout);}
#endif

  p->ax=v[3]-v[0];
  p->ay=v[4]-v[1];
  p->az=v[5]-v[2];
  d=1./sqrt(p->ax*p->ax+p->ay*p->ay+p->az*p->az);
  p->ax=p->ax*d;
  p->ay=p->ay*d;
  p->az=p->az*d;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  a=(%lf,%lf,%lf)\n",p->ax,p->ay,p->az);fflush(stdout);}
#endif

  p->bx=v[6]-v[0];
  p->by=v[7]-v[1];
  p->bz=v[8]-v[2];
  d=p->ax*p->bx+p->ay*p->by+p->az*p->bz;
  p->bx-=p->ax*d;
  p->by-=p->ay*d;
  p->bz-=p->az*d;
  d=1./sqrt(p->bx*p->bx+p->by*p->by+p->bz*p->bz);
  p->bx=p->bx*d;
  p->by=p->by*d;
  p->bz=p->bz*d;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  b=(%lf,%lf,%lf)\n",p->bx,p->by,p->bz);fflush(stdout);}
#endif

  p->nx=p->ay*p->bz-p->az*p->by;
  p->ny=p->az*p->bx-p->ax*p->bz;
  p->nz=p->ax*p->by-p->ay*p->bx;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("  n=(%lf,%lf,%lf)\n",p->nx,p->ny,p->nz);fflush(stdout);}
#endif

  MFIMFSetData(polygon,(void*)p,e);
  MFIMFSetFreeData(polygon,MFFreePolygonIn3SpaceData,e);
  MFIMFSetProject(polygon,MFProjectPolygonIn3Space,e);
  MFIMFSetTangent(polygon,MFTangentPolygonIn3Space,e);
  MFIMFSetScale(polygon,MFScalePolygonIn3Space,e);
  MFIMFSetR(polygon,R,e);
  MFIMFSetProjectForBB(polygon,MFPolygonIn3SpaceProjectForBB,e);
  MFIMFSetProjectForDraw(polygon,MFPolygonIn3SpaceProjectForBB,e);
  MFIMFSetWriteData(polygon,MFWritePolygonIn3SpaceData,e);

  MFIMFSetVectorFactory(polygon,MFNVectorFactory,e);
  MFIMFSetMatrixFactory(polygon,MFNKMatrixFactory,e);

/* Test */

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("Tangent [ %lf %lf %lf ]\n",
             (p->ax*p->ax+p->ay*p->ay+p->az*p->az),
             (p->ax*p->bx+p->ay*p->by+p->az*p->bz),
             (p->ax*p->nx+p->ay*p->ny+p->az*p->nz));
    printf("        [ %lf %lf %lf ]\n",
             (p->bx*p->ax+p->by*p->ay+p->bz*p->az),
             (p->bx*p->bx+p->by*p->by+p->bz*p->bz),
             (p->bx*p->nx+p->by*p->ny+p->bz*p->nz));
    printf("        [ %lf %lf %lf ]\n",
             (p->nx*p->ax+p->ny*p->ay+p->nz*p->az),
             (p->nx*p->bx+p->ny*p->by+p->nz*p->bz),
             (p->nx*p->nx+p->ny*p->ny+p->nz*p->nz));
    for(i=1;i<n;i++)
     {
      d=(p->nx*(v[3*i]-v[0])+p->ny*(v[3*i+1]-v[1])+p->nz*(v[3*i+2]-v[2]));
      printf("  (x[%d]-x[0]).n=%lf\n",i,d);
     }
    printf("done %s\n",RoutineName);fflush(stdout);
   }
#endif

  return polygon;
 }

void MFFreePolygonIn3SpaceData(void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreePolygonIn3SpaceData"};
  struct MFPolygonIn3SpaceData *p;

  p=(struct MFPolygonIn3SpaceData*)d;
  free(p);
  return;
 }

int MFProjectPolygonIn3Space(int n,int k,MFNVector vu0,MFNKMatrix mPhi,MFNVector vu,void *d,int *index, MFErrorHandler e)
 {
  static char RoutineName[]={"MFProjectPolygonIn3Space"};
  static int i;
  double a,b,nm;
  struct MFPolygonIn3SpaceData *p;
  double *u0,*Phi,*u;
  int verbose=0;

  u0=MFNV_CStar(vu0,e);
  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  p=(struct MFPolygonIn3SpaceData*)d;

  a=(u0[0]-p->ox)*p->ax+(u0[1]-p->oy)*p->ay+(u0[2]-p->oz)*p->az;
  b=(u0[0]-p->ox)*p->bx+(u0[1]-p->oy)*p->by+(u0[2]-p->oz)*p->bz;
  nm=(u0[0]-p->ox)*p->nx+(u0[1]-p->oy)*p->ny+(u0[2]-p->oz)*p->nz;

  u[0]=p->ox+a*p->ax+b*p->bx;
  u[1]=p->oy+a*p->ay+b*p->by;
  u[2]=p->oz+a*p->az+b*p->bz;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("  project (%lf,%lf,%lf)=(%lf,%lf,%lf)\n"
           "                       +%lf*(%lf,%lf,%lf)\n"
           "                       +%lf*(%lf,%lf,%lf)\n"
           "                       +%lf*(%lf,%lf,%lf)\n",u0[0],u0[1],u0[2],p->ox,p->oy,p->oz,a,p->ax,p->ay,p->az,b,p->bx,p->by,p->bz,nm,p->nx,p->ny,p->nz);
    printf("        =(%lf,%lf,%lf)\n",p->ox+a*p->ax+b*p->bx+nm*p->nx,p->oy+a*p->ay+b*p->by+nm*p->ny,p->oz+a*p->az+b*p->bz+nm*p->nz);
    printf("        =(%lf,%lf,%lf)\n",u[0],u[1],u[2]);
    printf("done %s\n",RoutineName);
    fflush(stdout);
   }
#endif

  *index=0;
  return 1;
 }

int MFTangentPolygonIn3Space(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFTangentPolygonIn3Space"};
  static int i;
  struct MFPolygonIn3SpaceData *p;
  double *Phi,*u;

  Phi=MFNKM_CStar(mPhi,e);
  u=MFNV_CStar(vu,e);

  p=(struct MFPolygonIn3SpaceData*)d;
  Phi[0]=p->ax;
  Phi[1]=p->ay;
  Phi[2]=p->az;
  Phi[3]=p->bx;
  Phi[4]=p->by;
  Phi[5]=p->bz;
  if(0)
   {
    printf("%s, \n",RoutineName);
    printf("    [%lf,%lf]\n",Phi[0],Phi[3]);
    printf("    [%lf,%lf]\n",Phi[1],Phi[4]);
    printf("    [%lf,%lf]\n",Phi[2],Phi[5]);
   }

  return 1;
 }

double MFScalePolygonIn3Space(int n,int k,MFNVector vu,MFNKMatrix mPhi,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFScalePolygonIn3Space"};
  struct MFPolygonIn3SpaceData *p;

  p=(struct MFPolygonIn3SpaceData*)d;
  return p->r;
 }

void MFWritePolygonIn3SpaceData(FILE *fid,void *d, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWritePolygonIn3SpaceData"};
  struct MFPolygonIn3SpaceData *data;

  data=(struct MFPolygonIn3SpaceData*)d;
  fprintf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",data->r,data->ox,data->oy,data->oz,data->ax,data->ay,data->az,data->bx,data->by,data->bz,data->nx,data->ny,data->nz);
  fflush(fid);
  return;
 }

MFImplicitMF MFReadPolygonIn3Space(FILE *fid, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadPolygonIn3Space"};
  MFImplicitMF polygon;
  struct MFPolygonIn3SpaceData *p;
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
  double bx=0.;
  double by=0.;
  double bz=0.;
  double nx=0.;
  double ny=0.;
  double nz=0.;
  double v[9];

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("%s\n",RoutineName);
    fflush(stdout);
   }
#endif

  fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&r,&ox,&oy,&oz,&ax,&ay,&az,&bx,&by,&bz,&nx,&ny,&nz);

  v[0]=ox;
  v[1]=oy;
  v[2]=oz;
  v[3]=ox+ax;
  v[4]=oy+ay;
  v[5]=oz+az;
  v[6]=ox+bx;
  v[7]=oy+by;
  v[8]=oz+bz;

  polygon= MFIMFCreatePolygonIn3SpaceWithRadius(3,v,r,e);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return polygon;
 }

int MFPolygonIn3SpaceProjectForBB(MFNVector u, double *x, void *d, MFErrorHandler e)
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
