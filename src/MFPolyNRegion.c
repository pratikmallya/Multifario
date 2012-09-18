/* 
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   November 11, 1997
 *              February 2, 1999   Ported to C
 *              October 6, 2000    Added Polyhedron
 */

static char *id="@(#) $Id: MFPolyNRegion.c,v 1.9 2011/07/21 17:42:46 mhender Exp $";

static char MFNRegionErrorMsg[256]="";

#include <MFNVector.h>
#include <MFNRegion.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#ifdef __cplusplus
 extern "C" {
#endif

static int CubeTest(MFNVector,void*,MFErrorHandler);
static void CubeFree(void*,MFErrorHandler);
static int RectangleTest(MFNVector,void*,MFErrorHandler);
static void RectangleFree(void*,MFErrorHandler);
static int HyperCubeTest(MFNVector,void*,MFErrorHandler);
static void HyperCubeFree(void*,MFErrorHandler);
static int HyperCubeTestByCorners(MFNVector,void*,MFErrorHandler);
static int DodIcosTest(MFNVector,void*,MFErrorHandler);
static void DodIcosFree(void*,MFErrorHandler);
static int MF3dPolygonTest(MFNVector,void*,MFErrorHandler);
static void MF3dPolygonFree(void*,MFErrorHandler);
static int MF3dEdgeTest(MFNVector,void*,MFErrorHandler);
static void MF3dEdgeFree(void*,MFErrorHandler);
static int MFPolyhedronTest(MFNVector,void*,MFErrorHandler);
static void MFPolyhedronFree(void*,MFErrorHandler);

MFNRegion MFNRegionCreateCube(double x0,double y0,double z0,double x1,double y1,double z1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreateCube"};
  MFNRegion cube;
  double *data;

  cube=MFNRegionCreateBaseClass("Cube",e);
  MFNRegionSetTest(cube,CubeTest,e);
  MFNRegionSetFreeData(cube,CubeFree,e);

  data=(double*)malloc(6*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",6*sizeof(double));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(cube);
    return NULL;
   }
#endif

  data[0]=x0;
  data[1]=y0;
  data[2]=z0;
  data[3]=x1;
  data[4]=y1;
  data[5]=z1;
  MFNRegionSetData(cube,(void*)data,e);

  return cube;
 }

int CubeTest(MFNVector v,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"CubeTest"};
  double *d;
  int rc;
  int n;
  double *x;

  n=MFNV_NC(v,e);
  x=MFNV_CStar(v,e);

  d=(double*)data;
  if(x[0]<d[0])rc=0;
   else if(x[1]<d[1])rc=0;
   else if(x[2]<d[2])rc=0;
   else if(x[0]>d[3])rc=0;
   else if(x[1]>d[4])rc=0;
   else if(x[2]>d[5])rc=0;
   else rc=1;

  return rc;
 }

void CubeFree(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"CubeFree"};

  free(data);
  return;
 }

MFNRegion MFNRegionCreateRectangle(double x0,double y0,double x1,double y1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreateRectangle"};
  MFNRegion rectangle;
  double *data;

  rectangle=MFNRegionCreateBaseClass("Rectangle",e);
  MFNRegionSetTest(rectangle,RectangleTest,e);
  MFNRegionSetFreeData(rectangle,RectangleFree,e);

  data=(double*)malloc(4*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",4*sizeof(double));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(rectangle);
    return NULL;
   }
#endif

  data[0]=x0;
  data[1]=y0;
  data[2]=x1;
  data[3]=y1;
  MFNRegionSetData(rectangle,data,e);

  return rectangle;
 }

int RectangleTest(MFNVector v,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"RectangleTest"};
  double *d;
  int rc;
  int n;
  double *x;

  n=MFNV_NC(v,e);
  x=MFNV_CStar(v,e);

  d=(double*)data;
  if(x[0]<d[0])rc=0;
   else if(x[1]<d[1])rc=0;
   else if(x[0]>d[2])rc=0;
   else if(x[1]>d[3])rc=0;
   else rc=1;

  return rc;
 }

void RectangleFree(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"RectangleFree"};

  free(data);
  return;
 }

MFNRegion MFNRegionCreateHyperCube(int n, double R, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreateHyperCube"};
  MFNRegion hypercube;
  double *data;

  hypercube=MFNRegionCreateBaseClass("HyperCube",e);
  MFNRegionSetTest(hypercube,HyperCubeTest,e);
  MFNRegionSetFreeData(hypercube,HyperCubeFree,e);

  data=(double*)malloc(sizeof(double));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(double));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(hypercube);
    return NULL;
   }
#endif

  data[0]=R;

  MFNRegionSetData(hypercube,data,e);

  return hypercube;
 }

MFNRegion MFNRegionCreateHyperCubeByCorners(int n, MFNVector LL, MFNVector UR, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreateHyperCubeByCorners"};
  MFNRegion hypercube;
  int i;
  double l,r;
  double *data;

  hypercube=MFNRegionCreateBaseClass("HyperCubeByCorners",e);
  MFNRegionSetTest(hypercube,HyperCubeTestByCorners,e);
  MFNRegionSetFreeData(hypercube,HyperCubeFree,e);

  data=(double*)malloc((2*n+1)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(data==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(double));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    free(hypercube);
    return NULL;
   }
#endif

  data[0]=n;
  for(i=0;i<n;i++)
   {
    l=MFNV_C(LL,i,e);
    r=MFNV_C(UR,i,e);
    if(l>r)
     {
      l=MFNV_C(UR,i,e);
      r=MFNV_C(LL,i,e);
     }
    data[2*i+1]=l;
    data[2*i+2]=r;
/*  printf("Direction %d range [%lf,%lf]\n",i,l,r);fflush(stdout);*/
   }

  MFNRegionSetData(hypercube,data,e);

  return hypercube;
 }

int HyperCubeTest(MFNVector v,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"HyperCubeTest"};
  double R;
  int i;
  int rc;
  int n;
  double *x;

  R=((double*)data)[0];
  n=MFNV_NC(v,e);
  if(!strcmp(MFNVGetId(v,e),"DENSE"))
   {
    x=MFNV_CStar(v,e);
    rc=1;
    for(i=0;i<n;i++)
     {
      if(x[i]>R)rc=0;
      if(x[i]<-R)rc=0;
     }
   }else{
    rc=1;
    for(i=0;i<n&&rc;i++)
     {
      if(MFNV_C(v,i,e)>R)rc=0;
      if(MFNV_C(v,i,e)<-R)rc=0;
     }
   }

  return rc;
 }

int HyperCubeTestByCorners(MFNVector v,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"HyperCubeTestByCorners"};
  double *c;
  int i;
  int rc;
  int n;
  double *x;

/*printf("%s\n",RoutineName);fflush(stdout);*/
  c=(double*)data;
  n=round(c[0]);
  if(!strcmp(MFNVGetId(v,e),"DENSE"))
   {
    x=MFNV_CStar(v,e);
    rc=1;
    for(i=0;i<n&&rc;i++)
     {
      if(x[i]<c[2*i+1])rc=0;
      if(x[i]>c[2*i+2])rc=0;
     }
   }else{
    rc=1;
    for(i=0;i<n&&rc;i++)
     {
/*printf("%d   %lf < %lf < %lf\n",i,c[2*i+2],MFNV_C(v,i,e),c[2*i+1]);fflush(stdout);*/
      if(MFNV_C(v,i,e)<c[2*i+1])rc=0;
      if(MFNV_C(v,i,e)>c[2*i+2])rc=0;
     }
   }

  return rc;
 }

void HyperCubeFree(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"HyperCubeFree"};

  free(data);
  return;
 }

MFNRegion MFNRegionCreateDodecahedronMinusIcosahedron(MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreateDodecahedronMinusIcosahedron"};
  MFNRegion diff;

  diff=MFNRegionCreateBaseClass("DodecahedronMinusIcosahedron",e);
  MFNRegionSetTest(diff,DodIcosTest,e);
  MFNRegionSetFreeData(diff,DodIcosFree,e);

  return diff;
 }

int DodIcosTest(MFNVector v,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"DodIcosTest"};
#include <Dodecahedron.h>
#include <Icosahedron.h>
  int i;
  double nx,ny,nz;
  double ox,oy,oz;
  double ax,ay,az;
  double bx,by,bz;
  double t;
  int n;
  double *x;

  n=MFNV_NC(v,e);
  x=MFNV_CStar(v,e);

  for(i=0;i<DODnF;i++)
   {
    ox=DODx[3*(DODf[DODnFV*i]-1)  ];
    oy=DODx[3*(DODf[DODnFV*i]-1)+1];
    oz=DODx[3*(DODf[DODnFV*i]-1)+2];
    ax=DODx[3*(DODf[DODnFV*i+1]-1)  ]-ox;
    ay=DODx[3*(DODf[DODnFV*i+1]-1)+1]-oy;
    az=DODx[3*(DODf[DODnFV*i+1]-1)+2]-oz;
    bx=DODx[3*(DODf[DODnFV*i+3]-1)  ]-ox;
    by=DODx[3*(DODf[DODnFV*i+3]-1)+1]-oy;
    bz=DODx[3*(DODf[DODnFV*i+3]-1)+2]-oz;
    nx=ay*bz-az*by;
    ny=az*bx-ax*bz;
    nz=ax*by-ay*bx;
    t=-ox*nx-oy*ny-oz*nz;
    if(t*((x[0]-ox)*nx+(x[1]-oy)*ny+(x[2]-oz)*nz)>0)
     {
      return 0;
     }
   }
  for(i=0;i<ICOSnF;i++)
   {
    ox=ICOSx[3*(ICOSf[ICOSnFV*i]-1)  ];
    oy=ICOSx[3*(ICOSf[ICOSnFV*i]-1)+1];
    oz=ICOSx[3*(ICOSf[ICOSnFV*i]-1)+2];
    ax=ICOSx[3*(ICOSf[ICOSnFV*i+1]-1)  ]-ox;
    ay=ICOSx[3*(ICOSf[ICOSnFV*i+1]-1)+1]-oy;
    az=ICOSx[3*(ICOSf[ICOSnFV*i+1]-1)+2]-oz;
    bx=ICOSx[3*(ICOSf[ICOSnFV*i+3]-1)  ]-ox;
    by=ICOSx[3*(ICOSf[ICOSnFV*i+3]-1)+1]-oy;
    bz=ICOSx[3*(ICOSf[ICOSnFV*i+3]-1)+2]-oz;
    nx=ay*bz-az*by;
    ny=az*bx-ax*bz;
    nz=ax*by-ay*bx;
    t=-ox*nx-oy*ny-oz*nz;
    if(t*((x[0]-ox)*nx+(x[1]-oy)*ny+(x[2]-oz)*nz)<0)
     {
      return 0;
     }
   }
  return 1;
 }

void DodIcosFree(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"DodIcosFree"};
  return;
 }

struct MF3dPolygon
  {
   int n;
   double *v;
   double x0,y0,z0;
   double a0,b0,c0;
   double a1,b1,c1;
  };

MFNRegion MFNRegionCreatePolygonal3dRegion(int n, double *v, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreatePolygonal3dRegion"};
  MFNRegion polygon;
  struct MF3dPolygon *p;
  int i;
  double d;

  polygon=MFNRegionCreateBaseClass("3dPolygon",e);
  MFNRegionSetTest(polygon,MF3dPolygonTest,e);
  MFNRegionSetFreeData(polygon,MF3dPolygonFree,e);

  p=(struct MF3dPolygon*)malloc(sizeof(struct MF3dPolygon));

#ifndef MFNOSAFETYNET
  if(p==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MF3dPolygon));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  p->n=n;
  p->v=(double*)malloc(3*(n+1)*sizeof(double));

#ifndef MFNOSAFETYNET
  if(p->v==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",3*(n+1)*sizeof(double));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  for(i=0;i<3*n;i++)(p->v)[i]=v[i];
  (p->v)[3*n]=0.;
  (p->v)[3*n+1]=0.;
  (p->v)[3*n+2]=0.;
  p->x0=0.;
  p->y0=0.;
  p->z0=0.;
  for(i=0;i<n;i++)
   {
    (p->v)[3*n]+=v[3*i];
    (p->v)[3*n+1]+=v[3*i+1];
    (p->v)[3*n+2]+=v[3*i+2];
    p->x0+=v[3*i];
    p->y0+=v[3*i+1];
    p->z0+=v[3*i+2];
   }
  (p->v)[3*n]/=n;
  (p->v)[3*n+1]/=n;
  (p->v)[3*n+2]/=n;
  p->x0/=n;
  p->y0/=n;
  p->z0/=n;

  p->a0=v[0]-p->x0;
  p->b0=v[1]-p->y0;
  p->c0=v[2]-p->z0;
  d=1./sqrt(p->a0*p->a0+p->b0*p->b0+p->c0*p->c0);
  p->a0*=d;
  p->b0*=d;
  p->c0*=d;

  p->a1=v[3]-p->x0;
  p->b1=v[4]-p->y0;
  p->c1=v[5]-p->z0;

  d=p->a1*p->a0+p->b1*p->b0+p->c1*p->c0;
  p->a1-=d*p->a0;
  p->b1-=d*p->b0;
  p->c1-=d*p->c0;
  d=1/sqrt(p->a1*p->a1+p->b1*p->b1+p->c1*p->c1);
  p->a1*=d;
  p->b1*=d;
  p->c1*=d;

  MFNRegionSetData(polygon,(void*)p,e);

  return polygon;
 }

int MF3dPolygonTest(MFNVector v,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MF3dPolygonTest"};
  int i,j;
  double ax,ay,az;
  double bx,by,bz;
  double sum;
  double rx,ry,rz;
  struct MF3dPolygon *p;
  double a,b;
  double alpha,beta;
  double n1,n2;
  double aBar,bBar,cBar;
  double integralTmp,q,fluxX,fluxY;
  double x1,y1,x2,y2,xp,yp;
  double exact;
  int n;
  double *x;
  int verbose=0;

  n=MFNV_NC(v,e);
  x=MFNV_CStar(v,e);

  p=(struct MF3dPolygon*)data;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s (%lf,%lf,%lf)\n",RoutineName,x[0],x[1],x[2]);fflush(stdout);}
#endif

  xp=(x[0]-p->x0)*p->a0+(x[1]-p->y0)*p->b0+(x[2]-p->z0)*p->c0;
  yp=(x[0]-p->x0)*p->a1+(x[1]-p->y0)*p->b1+(x[2]-p->z0)*p->c1;

#ifdef MFALLOWVERBOSE
  if(verbose)
   {
    printf("         (%lf,%lf)\n",xp,yp);fflush(stdout);
    printf("         n=%d\n",p->n);fflush(stdout);
   }
#endif

  sum=0.;
  for(i=0;i<p->n;i++)
   {
    ax=(p->v)[3*i  ];
    ay=(p->v)[3*i+1];
    az=(p->v)[3*i+2];
    j=i+1;if(j>p->n-1)j=0;
    bx=(p->v)[3*j  ];
    by=(p->v)[3*j+1];
    bz=(p->v)[3*j+2];

    x1=(ax-p->x0)*p->a0+(ay-p->y0)*p->b0+(az-p->z0)*p->c0;
    y1=(ax-p->x0)*p->a1+(ay-p->y0)*p->b1+(az-p->z0)*p->c1;
    x2=(bx-p->x0)*p->a0+(by-p->y0)*p->b0+(bz-p->z0)*p->c0;
    y2=(bx-p->x0)*p->a1+(by-p->y0)*p->b1+(bz-p->z0)*p->c1;

    a=(x1-xp);
    b=(y1-yp);
    alpha=(x2-x1);
    beta =(y2-y1);
    n1=beta;
    n2=-alpha;

    aBar=a*a+b*b;
    bBar=2.*(a*alpha+b*beta);
    cBar=alpha*alpha+beta*beta;
    q=4.*aBar*cBar-bBar*bBar;
    if(b*alpha==a*beta && bBar>=0 && bBar<=aBar)return 1;

#ifdef MFALLOWVERBOSE
    if(verbose){printf("    edge %d (%lf,%lf)<->%d (%lf,%lf), q=%le\n",i,x1,y1,j,x2,y2,q);fflush(stdout);}
#endif

    if(q!=0.)
     {
      integralTmp=2./sqrt(q)*(atan2(2.*cBar+bBar,sqrt(q))
                             -atan2(bBar,sqrt(q)));
     }else{
         integralTmp=2.*(1./bBar-1./(2*cBar+bBar));
     }
    exact=a*integralTmp+.5*alpha/cBar*
          (log(fabs(aBar+bBar+cBar))-log(aBar))
        -.5*alpha*bBar/cBar*integralTmp;
    fluxX=n1*exact;
    exact=b*integralTmp+.5*beta/cBar*
          (log(fabs(aBar+bBar+cBar))-log(aBar))
        -.5*beta*bBar/cBar*integralTmp;
    fluxY=n2*exact;
    sum+=fluxX+fluxY;
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("         sum=%le, result=%d\n",sum,fabs(sum)>1.5*3.1415926);fflush(stdout);}
#endif

  return fabs(sum)>1.5*3.1415926;
 }

void MF3dPolygonFree(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MF3dPolygonFree"};
  struct MF3dPolygon *p;

  p=(struct MF3dPolygon*)data;

  free(p->v);
  free(p);
  return;
 }

#define TRIANGULATED
#ifdef TRIANGULATED

struct MFPolyhedron
  {
   int nv;
   double *v;
   int nt;
   int *t;
   double xinf,yinf,zinf;
  };

int MFCreateTriangulatedPolygon(int nv,double *v,      int nfv, int *fv, int **tmp, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateTriangulatedPolygon"};
  int i;

/* This assumes the face is convex */

/*   Could do better by identifying the "fattest" convex triangle */

  *tmp=(int*)realloc(*tmp,3*(nfv-2)*sizeof(int));
  for(i=0;i<nfv-2;i++)
   {
    (*tmp)[3*i+0]=fv[0];
    (*tmp)[3*i+1]=fv[i+1];
    (*tmp)[3*i+2]=fv[i+2];
   }

  return nfv-2;
 }

MFNRegion MFNRegionCreatePolyhedral3dRegion(int nv, double *v, int nf, int *nfv, int **fv, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreatePolyhedral3dRegion"};
  MFNRegion polyhedron;
  struct MFPolyhedron *p;
  int i,n;
  int j;
  double d;
  int *tmp;

  polyhedron=MFNRegionCreateBaseClass("Polyhedron",e);
  MFNRegionSetTest(polyhedron,MFPolyhedronTest,e);
  MFNRegionSetFreeData(polyhedron,MFPolyhedronFree,e);

  p=(struct MFPolyhedron*)malloc(sizeof(struct MFPolyhedron));

#ifndef MFNOSAFETYNET
  if(p==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MFPolyhedron));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  p->nv=nv;
  p->v=(double*)malloc(3*p->nv*sizeof(double));

#ifndef MFNOSAFETYNET
  if(p->v==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",3*n*sizeof(double));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  p->xinf=v[0];
  p->yinf=v[0];
  p->zinf=v[0];
  for(i=0;i<nv;i++)
   {
    (p->v)[3*i+0]=v[3*i+0];
    (p->v)[3*i+1]=v[3*i+1];
    (p->v)[3*i+2]=v[3*i+2];
    if(v[3*i+0]>p->xinf)p->xinf=v[3*i+0];
    if(v[3*i+1]>p->yinf)p->yinf=v[3*i+1];
    if(v[3*i+2]>p->zinf)p->zinf=v[3*i+2];
   }
/*
  p->xinf=p->xinf+20;
  p->yinf=p->yinf+20;
  p->zinf=p->zinf+20;
 */
/* Break faces up into triangles */

  tmp=NULL;
  p->t=NULL;
  p->nt=0;
  for(i=0;i<nf;i++)
   {
    n=MFCreateTriangulatedPolygon(nv,v,nfv[i],fv[i],&tmp,e);

    p->t=realloc((void*)(p->t),3*(p->nt+n)*sizeof(int));

#ifndef MFNOSAFETYNET
    if(p->t==NULL)
     {
      sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",3*(p->nt+n)*sizeof(double));
      MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return NULL;
     }
#endif

    for(j=0;j<3*n;j++)p->t[3*p->nt+j]=tmp[j];
    p->nt+=n;
   }
  free(tmp);

  MFNRegionSetData(polyhedron,(void*)p,e);

  return polyhedron;
 }
#endif

int MFPolyhedronTest(MFNVector v,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolyhedronTest"};
  int i,j;
  double sum,flux,vol;
  double r1,r1x,r1y,r1z;
  double r2,r2x,r2y,r2z;
  double r3,r3x,r3y,r3z;
  struct MFPolyhedron *p;
  double a,b;
  double area;
  int n;
  double *x;
  int verbose=0;

  n=MFNV_NC(v,e);
  x=MFNV_CStar(v,e);

  p=(struct MFPolyhedron*)data;

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s (%lf,%lf,%lf)  ",RoutineName,x[0],x[1],x[2]);fflush(stdout);}
#endif
#if 0
  sum=0.;
  for(i=0;i<p->nt;i++)
   {
    r1x=p->v[3*(p->t[3*i  ])  ] -x[0];
    r1y=p->v[3*(p->t[3*i  ])+1] -x[1];
    r1z=p->v[3*(p->t[3*i  ])+2] -x[2];
    r1=sqrt(r1x*r1x+r1y*r1y+r1z*r1z);
    r2x=p->v[3*(p->t[3*i+1])  ] -x[0];
    r2y=p->v[3*(p->t[3*i+1])+1] -x[1];
    r2z=p->v[3*(p->t[3*i+1])+2] -x[2];
    r2=sqrt(r2x*r2x+r2y*r2y+r2z*r2z);
    r3x=p->v[3*(p->t[3*i+2])  ] -x[0];
    r3y=p->v[3*(p->t[3*i+2])+1] -x[1];
    r3z=p->v[3*(p->t[3*i+2])+2] -x[2];
    r3=sqrt(r3x*r3x+r3y*r3y+r3z*r3z);
    flux=(r1*r2*r3+r1*(r2x*r3x+r2y*r3y+r2z*r3z)
                  +r2*(r3x*r1x+r3y*r1y+r3z*r1z)
                  +r3*(r1x*r2x+r1y*r2y+r1z*r2z))/
         sqrt(2*(r2*r3+r2x*r3x+r2y*r3y+r2z*r3z)
               *(r3*r1+r3x*r1x+r3y*r1y+r3z*r1z)
               *(r1*r2+r1x*r2x+r1y*r2y+r1z*r2z));
    vol=r1x*(r2y*r3z-r2z*r3y)
       +r1y*(r2z*r3x-r2x*r3z)
       +r1z*(r2x*r3y-r2y*r3x);
    flux=2*acos(flux);
    if(vol<0.)flux=-flux;
  if(verbose){printf("   flux %d=%lf\n",i,flux);fflush(stdout);}
    sum+=flux;
   }

#ifdef MFALLOWVERBOSE
  if(verbose){printf("         sum=%le, result=%d\n",sum,fabs(sum)>1.5*3.1415926);fflush(stdout);}
#endif
  if( fabs(sum)>2.*3.1415926 )
   {
    printf("%s (%14.7lf,%14.7lf,%14.7lf)  sum=%21.7le>?%21.7lf, result=%d\n",RoutineName,x[0],x[1],x[2],fabs(sum),1.5*3.1415926,fabs(sum)>1.5*3.1415926);fflush(stdout);
   }

  return fabs(sum)>2.*3.1415926;
#endif

{
  double dx,dy,dz;
  double nx,ny,nz;
  double v0x,v0y,v0z;
  double v1x,v1y,v1z;
  double v2x,v2y,v2z;
  double t;
  double a11,a12,a21,a22,b1,b2,x1,x2;
  int result;
  int onBnd;
  double eps=.01;
  int try;
  int verbose=0;

/*if(0.<=x[0] && x[0]<=1. && 0.<=x[1] && x[1]<=1. && 0.<=x[2] && x[2]<=1. )verbose=1;*/

  onBnd=1;
  try=0;
  while(onBnd)
   {
  onBnd=0;

  dx=p->xinf-x[0]+40.*rand()/RAND_MAX;
  dy=p->yinf-x[1]+40.*rand()/RAND_MAX;
  dz=p->zinf-x[2]+40.*rand()/RAND_MAX;

  sum=0.;

  if(verbose)
   {
    printf("%s try %d x=(%lf,%lf,%lf)  xinf=(%lf,%lf,%lf)\n",RoutineName,try,x[0],x[1],x[2],x[0]+dx,x[1]+dy,x[2]+dz);fflush(stdout);

    t=-x[2]/dz;    printf("    z=0 at (%11.7lf,%11.7lf,%11.7lf) (%11.7lf,%11.7lf,%11.7lf) \n",x[0]+t*dx,x[1]+t*dy,x[2]+t*dz,x[0]+t*dx,x[1]+t*dy,x[2]+t*dz);fflush(stdout);
    t=(1-x[2])/dz; printf("    z=1 at (%11.7lf,%11.7lf,%11.7lf) (%11.7lf,%11.7lf,%11.7lf) \n",x[0]+t*dx,x[1]+t*dy,x[2]+t*dz,x[0]+t*dx,x[1]+t*dy,x[2]+t*dz);fflush(stdout);
    
    t=-x[1]/dy;    printf("    y=0 at (%11.7lf,%11.7lf,%11.7lf) (%11.7lf,%11.7lf,%11.7lf) \n",x[0]+t*dx,x[1]+t*dy,x[2]+t*dz,x[0]+t*dx,x[1]+t*dy,x[2]+t*dz);fflush(stdout);
    t=(1-x[1])/dy; printf("    y=1 at (%11.7lf,%11.7lf,%11.7lf) (%11.7lf,%11.7lf,%11.7lf) \n",x[0]+t*dx,x[1]+t*dy,x[2]+t*dz,x[0]+t*dx,x[1]+t*dy,x[2]+t*dz);fflush(stdout);
    
    t=-x[0]/dx;    printf("    x=0 at (%11.7lf,%11.7lf,%11.7lf) (%11.7lf,%11.7lf,%11.7lf) \n",x[0]+t*dx,x[1]+t*dy,x[2]+t*dz,x[0]+t*dx,x[1]+t*dy,x[2]+t*dz);fflush(stdout);
    t=(1-x[0])/dx; printf("    x=1 at (%11.7lf,%11.7lf,%11.7lf) (%11.7lf,%11.7lf,%11.7lf) \n",x[0]+t*dx,x[1]+t*dy,x[2]+t*dz,x[0]+t*dx,x[1]+t*dy,x[2]+t*dz);fflush(stdout);
   }

try++;
  for(i=0;i<p->nt;i++)
   {
    v0x=p->v[3*(p->t[3*i+0])  ];
    v0y=p->v[3*(p->t[3*i+0])+1];
    v0z=p->v[3*(p->t[3*i+0])+2];

    v1x=p->v[3*(p->t[3*i+1])  ];
    v1y=p->v[3*(p->t[3*i+1])+1];
    v1z=p->v[3*(p->t[3*i+1])+2];

    v2x=p->v[3*(p->t[3*i+2])  ];
    v2y=p->v[3*(p->t[3*i+2])+1];
    v2z=p->v[3*(p->t[3*i+2])+2];

    nx=(v1z-v0z)*(v2y-v0y)-(v1y-v0y)*(v2z-v0z);
    ny=(v1x-v0x)*(v2z-v0z)-(v1z-v0z)*(v2x-v0x);
    nz=(v1y-v0y)*(v2x-v0x)-(v1x-v0x)*(v2y-v0y);

    t=-( (x[0]-v0x)*nx+(x[1]-v0y)*ny+(x[2]-v0z)*nz )/( dx*nx+dy*ny+dz*nz );

    a11=(v1x-v0x)*(v1x-v0x)+(v1y-v0y)*(v1y-v0y)+(v1z-v0z)*(v1z-v0z);
    a12=(v1x-v0x)*(v2x-v0x)+(v1y-v0y)*(v2y-v0y)+(v2z-v0z)*(v1z-v0z);
    a21=(v2x-v0x)*(v1x-v0x)+(v2y-v0y)*(v1y-v0y)+(v1z-v0z)*(v2z-v0z);
    a22=(v2x-v0x)*(v2x-v0x)+(v2y-v0y)*(v2y-v0y)+(v2z-v0z)*(v2z-v0z);

    b1= (x[0]+t*dx-v0x)*(v1x-v0x) + (x[1]+t*dy-v0y)*(v1y-v0y) + (x[2]+t*dz-v0z)*(v1z-v0z);
    b2= (x[0]+t*dx-v0x)*(v2x-v0x) + (x[1]+t*dy-v0y)*(v2y-v0y) + (x[2]+t*dz-v0z)*(v2z-v0z);

    if( fabs(a11*a22-a12*a21)<1.e-7)
     {
      x1=-1.;
      x2=-1.;
     }else{
      x1=( a22*b1-a12*b2)/(a11*a22-a12*a21);
      x2=(-a21*b1+a11*b2)/(a11*a22-a12*a21);
     }

    if(fabs(a11*x1+a12*x2-b1)+fabs(a21*x1+a22*x2-b2)>1.e-7)
     {printf("*****A.x-b = (%11.7lf,%11.7lf)\n",a11*x1+a12*x2-b1,a21*x1+a22*x2-b2);fflush(stdout);}

    if(fabs(a11*a22-a12*a21)<1.e-7)
     {printf("*****determinant zero\n");fflush(stdout);
      printf("    [ %11.7lf %11.7lf ] [ x1 ] = [ %11.7lf ]\n",a11,a12,b1);fflush(stdout);
      printf("    [ %11.7lf %11.7lf ] [ x2 ] = [ %11.7lf ]\n",a21,a22,b2);fflush(stdout);exit(12);}

    if(x1!=x1 || x2!=x2)
     {printf("*****solution was nan\n");fflush(stdout);
      printf("    [ %11.7lf %11.7lf ] [ x1 ] = [ %11.7lf ]\n",a11,a12,b1);fflush(stdout);
      printf("    [ %11.7lf %11.7lf ] [ x2 ] = [ %11.7lf ]\n",a21,a22,b2);fflush(stdout);exit(12);}
    if(fabs(v0x+x1*(v1x-v0x)+x2*(v2x-v0x)-(x[0]+t*dx))
      +fabs(v0y+x1*(v1y-v0y)+x2*(v2y-v0y)-(x[1]+t*dy))
      +fabs(v0z+x1*(v1z-v0z)+x2*(v2z-v0z)-(x[2]+t*dz))>1.e-7)
     {
      printf(" (%11.7lf,%11.7lf,%11.7lf)!=(%11.7lf,%11.7lf,%11.7lf)\n",
                    v0x+x1*(v1x-v0x)+x2*(v2x-v0x),
                    v0y+x1*(v1y-v0y)+x2*(v2y-v0y),
                    v0z+x1*(v1z-v0z)+x2*(v2z-v0z),
                    x[0]+t*dx,
                    x[1]+t*dy,
                    x[2]+t*dz);fflush(stdout);
      exit(12);
     }

    if(0.<=t && t<=1.)
     {
      if(x1+x2<1.&&x1>0.&&x2>0.)sum=sum+1;
      if(fabs(x1)<eps && x2>-eps && x2<1+eps )onBnd=1;
      if(fabs(x2)<eps && x1>-eps && x1<1+eps )onBnd=1;
      if(fabs(x1+x2-1.)<eps && x1>-eps && x2>-eps ) onBnd=1;
     }

    if(verbose)
     {
      printf(" tri %2d, (x1,x2)=(%11.7lf,%11.7lf) x1+x2=%11.7lf sum=%11.7lf onBnd=%d",i,
                          x1,x2,x1+x2,
                          sum,onBnd);fflush(stdout);
      printf("                t=%11.7lf  x=(%11.7lf,%11.7lf,%11.7llf)\n",t,
                    v0x+x1*(v1x-v0x)+x2*(v2x-v0x),
                    v0y+x1*(v1y-v0y)+x2*(v2y-v0y),
                    v0z+x1*(v1z-v0z)+x2*(v2z-v0z));
     }
   }

  }
  if(sum-2*floor(sum/2.)<.5)result=0;
   else                     result=1;

/*if(result!=-1){printf("done (%14.7lf,%14.7lf,%14.7lf) trys=%d sum=%21.7le result=%d\n",x[0],x[1],x[2],try,sum,result);fflush(stdout);}
  if(verbose){printf("***********\n");fflush(stdout);}
 */
  return result;
}
 }

void MFPolyhedronFree(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPolyhedronFree"};
  struct MFPolyhedron *p;

  p=(struct MFPolyhedron*)data;

  free(p->v);
  free(p);
  return;
 }

struct MF3dEdge
  {
   double x0,y0,z0;
   double a,b,c;
   double d;
  };

MFNRegion MFNRegionCreateEdge3dRegion(double *v0, double *v1, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNRegionCreateEdge3dRegion"};
  MFNRegion edge;
  struct MF3dEdge *p;
  int i;
  double d;

  edge=MFNRegionCreateBaseClass("3dEdge",e);
  MFNRegionSetTest(edge,MF3dEdgeTest,e);
  MFNRegionSetFreeData(edge,MF3dEdgeFree,e);

  p=(struct MF3dEdge*)malloc(sizeof(struct MF3dEdge));
  if(p==NULL)
   {
    sprintf(MFNRegionErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct MF3dEdge));
    MFSetError(e,12,RoutineName,MFNRegionErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
  p->x0=v0[0];
  p->y0=v0[1];
  p->z0=v0[2];
  p->a=v1[0]-v0[0];
  p->b=v1[1]-v0[1];
  p->c=v1[2]-v0[2];
  p->d=p->a*p->a+p->b*p->b+p->c*p->c;

  MFNRegionSetData(edge,(void*)p,e);

  return edge;
 }

int MF3dEdgeTest(MFNVector v,void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MF3dEdgeTest"};
  int i,j;
  double s;
  int verbose=0;
  int n;
  double *x;
  struct MF3dEdge *p;

  n=MFNV_NC(v,e);
  x=MFNV_CStar(v,e);

  p=(struct MF3dEdge*)data;

  s=((x[0]-p->x0)*p->a+(x[1]-p->y0)*p->b+(x[2]-p->z0)*p->c)/p->d;

  return s*(1-s)>-1.e-7;
 }

void MF3dEdgeFree(void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"MF3dEdgeFree"};
  struct MF3dEdge *p;

  p=(struct MF3dEdge*)data;

  free(p);
  return;
 }

#ifdef __cplusplus
}
#endif
