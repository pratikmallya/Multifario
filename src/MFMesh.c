/*      author: Mike Henderson mhender@watson.ibm.com */

#include <MFMesh.h>

#include <MFAtlas.h>
#include <MFAtlasFriends.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <MFEnumPolytope.h>
#include <MFEnumDualPolytope.h>
#include <MFMultifariosMethod.h>

#include <math.h>
#include <MFDraw.h>
#include <sh.h>

/* nV,v are all of the vertices */
/* nE,e are all of the edges */
/* nF,nFe, fe are the faces in the bounding polygon */

static checkSize(double **,int,int*,MFErrorHandler);
static char MFErrorMsg[1024];

int MFAtlasLeftPolytope(MFAtlas A,int i, MFErrorHandler e);
int MFAtlasRightPolytope(MFAtlas A,int i, MFErrorHandler e);
int MFPolytopeVertexSmallestIndex(MFPolytope,int,MFErrorHandler);
int MFPolytopeVertexNumberOfIndices(MFPolytope,int,MFErrorHandler);
int *MFPolytopeVertexIndexSet(MFPolytope,int,MFErrorHandler);

#define ANIMATE

MFEnumDualPolytope MFGenerateMeshInPolyhedron(int nV,double **v,int nE,int **ev, int **evs, int nF, int *nFe, int **fe,int **fes, double R,MFErrorHandler e)
 {
  static char RoutineName[]={"MFGenerateMeshInPolyhedron"};
  MFImplicitMF M;
  MFAtlas A;
  MFNVector u;
  MFNKMatrix Phi;
  double o[3]={0.,0.,0.};
  double b[9]={1.,0.,0.,
               0.,1.,0.,
               0.,0.,1.};
  MFEnumDualPolytope mesh;
  MFChart chart;
  int nMeshPts;
  int mMeshPts;
  double *meshPts;
  double d,dt,t;
  int i,j;
  int *ePts;
  int *nEPts;

  int ne;
  int nv;
  double *V;
  int n;
  int nf;
  int *nfv;
  int **fv;
  int k;
  int npt;
  MFNRegion Omega;
  double delta;

  int nPt0;
  double *Pt0;
  double *R0;
  double *Pt;
  double *RPt;
  int    nPt;

  char filename[1024];
  int full=255;
  int zero=  0;
  float x0,y0,z0,x1,y1,z1;

  MFContinuationMethod H;

printf("Enter %s\n",RoutineName);fflush(stdout);

  nMeshPts=0;
  mMeshPts=0;
  meshPts=(double*)NULL;
  checkSize(&meshPts,3*nMeshPts,&mMeshPts,e);

  ePts=(int*)malloc(nE*sizeof(int));
  if(ePts==(int*)NULL)
   {
    sprintf(MFErrorMsg,"out of memory allocating an array of size %d.",nE*sizeof(int));
    MFSetError(e,12,RoutineName,MFErrorMsg,__LINE__,__FILE__);
   }

  nEPts=(int*)malloc(nE*sizeof(int));
  if(nEPts==(int*)NULL)
   {
    sprintf(MFErrorMsg,"out of memory allocating an array of size %d.",nE*sizeof(int));
    MFSetError(e,12,RoutineName,MFErrorMsg,__LINE__,__FILE__);
   }

/* Place a ball at each vertex */
printf("Add points at vertices, nV=%d\n",nV);fflush(stdout);

  for(i=0;i<nV;i++)
   {
    checkSize(&meshPts,3*(nMeshPts+1),&mMeshPts,e);
    meshPts[3*nMeshPts  ]=v[i][0];
    meshPts[3*nMeshPts+1]=v[i][1];
    meshPts[3*nMeshPts+2]=v[i][2];
    nMeshPts++;
   }

/* Distribute points along each edge */
printf("Add points on edges, nE=%d\n",nE);fflush(stdout);

  for(i=0;i<nE;i++)
   {
    ePts[i]=nMeshPts;
    d=sqrt( (v[ev[i][1]][0]-v[ev[i][0]][0])*(v[ev[i][1]][0]-v[ev[i][0]][0])
           +(v[ev[i][1]][1]-v[ev[i][0]][1])*(v[ev[i][1]][1]-v[ev[i][0]][1])
           +(v[ev[i][1]][2]-v[ev[i][0]][2])*(v[ev[i][1]][2]-v[ev[i][0]][2]));
    n=floor(d/R);
    dt=1./n;
    nEPts[i]=0;
    for(j=1;j<n;j++)
     {
      t=j*dt;
      checkSize(&meshPts,3*(nMeshPts+1),&mMeshPts,e);
      meshPts[3*nMeshPts  ]=v[ev[i][0]][0] + t*(v[ev[i][1]][0]-v[ev[i][0]][0]);
      meshPts[3*nMeshPts+1]=v[ev[i][0]][1] + t*(v[ev[i][1]][1]-v[ev[i][0]][1]);
      meshPts[3*nMeshPts+2]=v[ev[i][0]][2] + t*(v[ev[i][1]][2]-v[ev[i][0]][2]);
      nMeshPts++;
      nEPts[i]++;
     }
   }

/* Distribute points on each face */
printf("Add points on faces, nF=%d\n",nF);fflush(stdout);

  for(i=0;i<nF;i++)
   {
    nv=nFe[i];
    V=(double*)malloc(3*nv*sizeof(double));

    nPt0=nv;
    for(j=0;j<nFe[i];j++)
     {
      nPt0+=nEPts[fe[i][j]];
     }
    Pt0=(double*)malloc(3*nPt0*sizeof(double));
    R0=(double*)malloc(nPt0*sizeof(double));
    nPt0=0;
    for(j=0;j<nFe[i];j++)
     {
      if(fes[i][j]>0)
       {
        Pt0[3*nPt0  ]=v[ev[fe[i][j]][0]][0];
        Pt0[3*nPt0+1]=v[ev[fe[i][j]][0]][1];
        Pt0[3*nPt0+2]=v[ev[fe[i][j]][0]][2];
       }else{
        Pt0[3*nPt0  ]=v[ev[fe[i][j]][1]][0];
        Pt0[3*nPt0+1]=v[ev[fe[i][j]][1]][1];
        Pt0[3*nPt0+2]=v[ev[fe[i][j]][1]][2];
       }
      R0[nPt0]=R;
      V[3*j  ]=Pt0[3*nPt0  ];
      V[3*j+1]=Pt0[3*nPt0+1];
      V[3*j+2]=Pt0[3*nPt0+2];
      nPt0++;
     }
    for(j=0;j<nFe[i];j++)
     {
      for(k=0;k<nEPts[fe[i][j]];k++)
       {
        Pt0[3*nPt0  ]=meshPts[3*(ePts[fe[i][j]]+k)  ];
        Pt0[3*nPt0+1]=meshPts[3*(ePts[fe[i][j]]+k)+1];
        Pt0[3*nPt0+2]=meshPts[3*(ePts[fe[i][j]]+k)+2];
        R0[nPt0]=R;
        nPt0++;
       }
     }

    Pt=(double*)NULL;
    RPt=(double*)NULL;
    nPt=MFComputePointsOn3dPolygon(R,nv,V,nPt0,Pt0,R0,&Pt,&RPt,e);

    checkSize(&meshPts,3*(nMeshPts+nPt),&mMeshPts,e);
    for(k=0;k<nPt;k++)
     {
      meshPts[3*nMeshPts  ]=Pt[3*k  ];
      meshPts[3*nMeshPts+1]=Pt[3*k+1];
      meshPts[3*nMeshPts+2]=Pt[3*k+2];
      nMeshPts++;
     }
    free(V);
    free(Pt0);
    free(Pt);
    free(RPt);
   }

/* Distribute points on the interior */

printf("Add points on interior\n");fflush(stdout);

/* Create the polyhedron */

  nv=nV;
  V=(double*)malloc(3*nv*sizeof(double));
  for(i=0;i<nv;i++)
   {
    V[3*i  ]=meshPts[3*i  ];
    V[3*i+1]=meshPts[3*i+1];
    V[3*i+2]=meshPts[3*i+2];
   }

  nf=nF;
  nfv=(int*)malloc(nf*sizeof(int));
  fv=(int**)malloc(nf*sizeof(int));
  for(i=0;i<nf;i++)
   {
    nfv[i]=nFe[i];
    fv[i]=(int*)malloc(nfv[i]*sizeof(int));
    for(j=0;j<nFe[i];j++)
     {
      if(fes[i][j]>0)
       {
        fv[i][j]=ev[fe[i][j]][0];
       }else{
        fv[i][j]=ev[fe[i][j]][1];
       }
     }
   }

  Omega=MFNRegionCreatePolyhedral3dRegion(nv, V, nf, nfv, fv, e);

  free(V);
  for(i=0;i<nf;i++)
    free(fv[i]);
  free(fv);
  free(nfv);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);   /* Write polyhedra to a plotfile */
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e); /* Write points to a file */
  MFMultifarioSetFilename(H,"CubeMesh",e);

  M=MFIMFCreateFlat(3,3,o,b,e);
  A=MFCreateAtlas(M,e);

  u=MFCreateNVector(3,e);

#ifdef ANIMATE
  shSetOutputFormat("tiff");
  shSetOutputResolution(2048,2048);
  MFDrawInitializeCube(0.,89.8,-0.1,1.1,-0.1,1.1,-0.1,1.1,e);
#endif

  for(i=0;i<nMeshPts;i++)
   {
    u=MFCreateNVector(3,e);
    MFNVSetC(u,0,meshPts[3*i  ],e);
    MFNVSetC(u,1,meshPts[3*i+1],e);
    MFNVSetC(u,2,meshPts[3*i+2],e);
    Phi=MFIMFTangentSpace(M,u,e);

    chart=MFCreateChartWithCubeSize(M,u,Phi,R,10.,e);
    MFAtlasAddChart(A,chart,e); 
#ifdef ANIMATE
    MFDrawClear(e);
    shlinc(&full,&zero,&zero);
    for(j=0;j<nE;j++)
     {
      x0=v[ev[j][0]][0];
      y0=v[ev[j][0]][1];
      z0=v[ev[j][0]][2];
      x1=v[ev[j][1]][0];
      y1=v[ev[j][1]][1];
      z1=v[ev[j][1]][2];
      shline( &x0,&y0,&z0,&x1,&y1,&z1);
     }
    MFDrawAtlasOnce(A,e);
    sprintf(filename,"%s%4.4d","MeshPolyhedron",i);shSetOutputFilename(filename);
    MFDrawDisplay(e);
#endif

    MFFreeNVector(u,e);
    MFFreeNKMatrix(Phi,e);
   }

  u=MFCreateNVector(3,e);

  j=0;
  while((i=MFAtlasPointOnBoundaryInsideRegion(A,Omega,u,&Phi,&delta,e))>-1)
   {
    chart=MFCreateChartWithCubeSize(M,u,Phi,R,10.,e);
    MFAtlasAddChart(A,chart,e); 
j++;

#ifdef ANIMATE
    MFDrawClear(e);
    MFDrawAtlasOnce(A,e);
    sprintf(filename,"%s%4.4d","MeshPolyhedron",MFAtlasNumberOfCharts(A,e));shSetOutputFilename(filename);
    MFDrawDisplay(e);
#endif

    MFFreeNVector(u,e);
    MFFreeNKMatrix(Phi,e);
    u=MFCreateNVector(3,e);
   }
#ifdef ANIMATE
  MFDrawClose(e);
#endif

 {
  int kkk;
  MFPolytope P;
  int index;
  int l,r;
  int nindices;
  int *indices;
  int i0,i1,i2,i3;
  int valid;
  FILE *fid;
  FILE *did;

  fid=fopen("junk.dx","w");

  did=fopen("junkVerts.dx","w");
  for(i=0;i<MFAtlasNumberOfCharts(A,e);i++)
   {
    fprintf(did,"        %lf %lf %lf\n",MFNV_C(MFAtlasCenterOfChart(A,i,e),0,e),
                                   MFNV_C(MFAtlasCenterOfChart(A,i,e),1,e),
                                   MFNV_C(MFAtlasCenterOfChart(A,i,e),2,e) );fflush(stdout);
   }
  fclose(did);

  fprintf(fid,"object \"Vertices\" class array type float rank 1 shape 3 items %d data file junkVerts.dx\n",MFAtlasNumberOfCharts(A,e));fflush(stdout);
  fprintf(fid,"\n");fflush(stdout);

  did=fopen("junkEdges.dx","w");
  kkk=0;
  for(i=0;i<MFAtlasNumberOfCharts(A,e);i++)
   {
    P=MFChartPolytope(MFAtlasChart(A,i,e),e);
    for(j=0;j<MFPolytopeNumberOfFaces(P,e);j++)
     {
      index=MFPolytopeFaceIndex(P,j,e);
      l= MFAtlasLeftPolytope(A,index,e);
      r= MFAtlasRightPolytope(A,index,e);
      if(i<=l && i<=r && l!=r)
       {
        if(l<r){fprintf(did,"     %d %d\n",l,r);fflush(stdout);}
         else
               {fprintf(did,"     %d %d\n",r,l);fflush(stdout);}
        kkk++;
       }
     }
   }
  fclose(did);

  fprintf(fid,"object \"Edges\" class array type int rank 1 shape 2 items %d data file junkEdges.dx\n",kkk);fflush(stdout);

  fprintf(fid,"attribute \"ref\" string \"positions\"\n");fflush(stdout);
  fprintf(fid,"attribute \"element type\" string \"lines\"\n");fflush(stdout);
  fprintf(fid,"\n");fflush(stdout);

  did=fopen("junkTets.dx","w");
  kkk=0;
  for(i=0;i<MFAtlasNumberOfCharts(A,e);i++)
   {
    valid=1;
    P=MFChartPolytope(MFAtlasChart(A,i,e),e);
    for(j=0;j<MFPolytopeNumberOfVertices(P,e);j++)
     {
      nindices=MFPolytopeVertexNumberOfIndices(P,j,e);
       indices=MFPolytopeVertexIndexSet(P,j,e);

      for(k=0;k<nindices;k++)
       {
        l= MFAtlasLeftPolytope(A,indices[k],e);
        r= MFAtlasRightPolytope(A,indices[k],e);
        if(l==r)valid=0;
       }
      if(valid)
       {
        i0=-1;
        i1=-1;
        i2=-1;
        i3=-1;
        for(k=0;k<nindices;k++)
         {
           {
            l= MFAtlasLeftPolytope(A,indices[k],e);
            r= MFAtlasRightPolytope(A,indices[k],e);
            if(i0==-1)
             {
              if(l<r){i0=l;i1=r;}
               else  {i0=r;i1=l;}
             }
            if(i2==-1)
             {
              if(l!=i0 && l!=i1)i2=l;
              if(r!=i0 && r!=i1)i2=r;

              if(i2!=-1 && i2<i1){l=i1;i1=i2;i2=l;}
             }
            if(i3==-1)
             {
              if(l!=i0 && l!=i1&&l!=i2)i3=l;
              if(r!=i0 && r!=i1&&r!=i2)i3=r;
              if(i3!=-1 && i3<i1){l=i1;i1=i3;i3=l;}
              if(i3!=-1 && i3<i2){l=i2;i2=i3;i3=l;}
             }
           }
         }
if(i3==-1)
 {
      printf("Bad Tet [%3d,%3d,%3d,%3d]\n",i0,i1,i2,i3);fflush(stdout);
      for(k=0;k<nindices;k++)
       {
        l= MFAtlasLeftPolytope(A,indices[k],e);
        r= MFAtlasRightPolytope(A,indices[k],e);
        printf("   Edge [%d,%d]\n",l,r);fflush(stdout);
       }
 }else{
        fprintf(did,"     %d %d %d %d\n",i0,i1,i2,i3);fflush(fid);
        kkk++;
 }
       }
     }
   }
  fclose(did);

fprintf(fid,"object \"Tets\" class array type int rank 1 shape 4 items %d data file junkTets.dx\n",kkk);fflush(stdout);
fprintf(fid,"attribute \"ref\" string \"positions\"\n");fflush(stdout);
fprintf(fid,"attribute \"element type\" string \"tetrahedra\"\n");fflush(stdout);
fprintf(fid,"\n");fflush(stdout);
fprintf(fid,"object \"f0\" class field\n");fflush(stdout);
fprintf(fid,"component \"positions\" value \"Vertices\"\n");fflush(stdout);
fprintf(fid,"component \"connections\" value \"Edges\"\n");fflush(stdout);
fprintf(fid,"\n");fflush(stdout);
fprintf(fid,"object \"f1\" class field\n");fflush(stdout);
fprintf(fid,"component \"positions\" value \"Vertices\"\n");fflush(stdout);
fprintf(fid,"component \"connections\" value \"Tets\"\n");fflush(stdout);
fprintf(fid,"\n");fflush(stdout);
fprintf(fid,"object \"default\" class group\n");fflush(stdout);
fprintf(fid," member 0 \"f0\"\n");fflush(stdout);
fprintf(fid," member 1 \"f1\"\n");fflush(stdout);
fprintf(fid,"end\n");fflush(stdout);
fclose(fid);
 }

mesh=NULL;
printf("Enumerate dual\n");fflush(stdout);

  mesh=MFEnumDualOfAtlas(A,e);

  MFDualPolytopeToDXFile("CubeMeshDual",mesh,e);
  MFCloseAtlas(H,A,e);

  MFFreeAtlas(A,e);
  MFFreeImplicitMF(M,e);

printf("Exit %s\n",RoutineName);fflush(stdout);

  return mesh;
 }

static checkSize(double **array,int n, int *m, MFErrorHandler e)
 {
  static char RoutineName[]={"MFGenerateMeshInPolyhedron"};
  if(n>=*m)
   {
    (*m)+=1000;
    *array=(double*)realloc(*array,(*m)*sizeof(double));
    if(*array==(double*)NULL)
     {
      sprintf(MFErrorMsg,"out of memory allocating an array of size %d.",(*m)*sizeof(double));
      MFSetError(e,12,RoutineName,MFErrorMsg,__LINE__,__FILE__);
     }
   }
  return;
 }
