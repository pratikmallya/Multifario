#include <stdio.h>
#include <MFAtlas.h>
#include <MFDraw.h>
#include <MFErrorHandler.h>
#include <sh.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
 extern "C" {
#endif

void MFPlotfileToDX(FILE*,char*,MFErrorHandler);
void PFPolytopes(FILE*,int,int,int*,int*,int*,double**,int**,int*,int*,int**,int**,int*,int*,int**,int***,int**,MFErrorHandler);
void PFDualPolytopes(FILE*,int,int,int*,int*,int*,double**,int**,int*,int*,int**,int**,int*,int*,int**,int***,int**,MFErrorHandler);
void MFMeshToDXFile2(char*,int,int,double*,int*,int,int*,int*,int,int*,int**,int*,MFErrorHandler);

void MFPlotfileToDX(FILE *fid,char *name, MFErrorHandler e)
 {
  static char RoutineName[]={"PlotfileToDX"};
  int mV,iV,iV0;
  double *v=(double*)NULL;
  int *vt=(int*)NULL;
  int mE,iE;
  int *ev=(int*)NULL;
  int *et=(int*)NULL;
  int mF,iF;
  int *nfv=(int*)NULL;
  int **fv=(int**)NULL;
  int *ft=(int*)NULL;
  int i,j,l,n,m,k;
  int verbose=0;

  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}

  iV=0;
  iV0=0;
  mV=0;

  iE=0;
  mE=0;

  iF=0;
  mF=0;

  fscanf(fid,"Dimension of vertices, %d\n",&m);
  fscanf(fid,"Dimension of manifold, %d\n",&k);

  if(verbose){printf("Add in Polytopes\n");fflush(stdout);}
  PFPolytopes(fid,m,k,&iV,&iV0,&mV,&v,&vt,&iE,&mE,&ev,&et,&iF,&mF,&nfv,&fv,&ft,e);
  if(verbose){printf("  There are now %d pts, %d edges %d faces\n",iV,iE,iF);fflush(stdout);}
  if(verbose){printf("Dump it to %s\n",name);fflush(stdout);}
  MFMeshToDXFile2(name,iV,m,v,vt,iE,ev,et,iF,nfv,fv,ft,e);

  if(v!=(double*)NULL)free(v);
  if(nfv!=(int*)NULL)free(nfv);
  if(fv!=(int**)NULL)
   {
    for(i=0;i<iF;i++)
     {
      if(fv[i]!=(int*)NULL)free(fv[i]);
     }
    free(fv);
   }
  return;
 }

void PFPolytopes(FILE *fid,int m, int k, int *piV,int *piV0,int *pmV,double **pv,int **pvt,int *piE,int *pmE,int **pev,int **pet, int *piF,int *pmF,int **pnfv,int ***pfv, int **pft, MFErrorHandler e)
 {
  int mV,iV,iV0;
  double *v;
  int *vt;
  int mE,iE;
  int *ev;
  int *et;
  int mF,iF;
  int *nfv;
  int **fv;
  int *ft;
  int i,j,l,n;
  int ipoly,nv,ne,nf;
  int verbose=0;
  int iv,ie,iE0,je,np,v0,v1;
  int trip;
  int sing,bnd;
  char c;
  double R;

  iV=*piV;
  iV0=*piV0;
  mV=*pmV;
  v=*pv;
  vt=*pvt;
  iE=*piE;
  mE=*pmE;
  ev=*pev;
  et=*pet;
  iF=*piF;
  mF=*pmF;
  nfv=*pnfv;
  fv=*pfv;
  ft=*pft;

  while(!feof(fid))
   {
    fscanf(fid,"Polyhedron %d, ",&ipoly);
    c=fgetc(fid);
    ungetc(c,fid);
    if(c=='R')
     {
      fscanf(fid,"R=%lf, %d vertices, %d edges, %d faces,  boundary %d, singular %d\n",&R,&nv,&ne,&nf,&bnd,&sing);
     }else{
      fscanf(fid,"%d vertices, %d edges, %d faces,  boundary %d, singular %d\n",&nv,&ne,&nf,&bnd,&sing);
     }

    if(verbose){printf("Polyhedron %d, %d vertices, %d edges, %d faces\n",ipoly,nv,ne,nf);fflush(stdout);}

    while(iV+nv>=mV)
     {
      mV+=100;
      v=(double*)realloc((void*)v,m*mV*sizeof(double));
      vt=(int*)realloc((void*)vt,mV*sizeof(int));
     }
    iV0=iV;
    for(iv=0;iv<nv;iv++)
     {
      fscanf(fid,"Vertex %*d (%lf",&(v[0+m*iV]));
      for(i=1;i<m;i++)fscanf(fid,",%lf",&(v[i+m*iV]));
      fscanf(fid,"), %*d [%*[ 0-9+-.,]]\n");
      vt[iV]=0;
      if(verbose)
       {
        printf("  Vertex %d (%lf",iV,v[0+m*iV]);
        for(i=1;i<m;i++)printf(",%lf",v[i+m*iV]);
        printf(")\n");fflush(stdout);
       }
      iV++;
     }

    while(iE+ne>=mE)
     {
      mE+=100;
      ev=(int*)realloc((void*)ev,2*mE*sizeof(int));
      et=(int*)realloc((void*)et,mE*sizeof(int));
     }
    iE0=iE;
    for(ie=0;ie<ne;ie++)
     {
      if(k>1)
        fscanf(fid,"Edge %*d (%d,%d), %*d [%*[ 0-9+-.,]]\n",&v0,&v1);
      else if(k==1)
        fscanf(fid,"Edge %*d (%d,%d), %*d []\n",&v0,&v1);
       else
        fscanf(fid,"Edge %*d (%d,%d), %*d [%*[ 0-9+-.,]]\n",&v0,&v1);
      ev[0+2*iE]=iV0+v0;
      ev[1+2*iE]=iV0+v1;
      et[iE]=0;
      if(verbose){printf("  Edge %d (%d,%d)\n",iE,ev[0+2*iE],ev[1+2*iE]);fflush(stdout);}
      iE++;
     }

    for(i=0;i<nf;i++)fscanf(fid,"Face %*d neighbor %*d\n");

    if(k==2&&ne>2)
     {
      while(iF+ne>=mF)
       {
        mF+=100;
        nfv=(int*)realloc((void*)nfv,mF*sizeof(int));
        if(nfv==(int*)NULL){printf("Error allocating %d bytes for nfv\n",mF*sizeof(int));fflush(stdout);abort();}
        fv=(int**)realloc((void*)fv,mF*sizeof(int*));
        if(fv==(int**)NULL){printf("Error allocating %d bytes for fv\n",mF*sizeof(int*));fflush(stdout);abort();}
        for(i=mF-100;i<mF;i++){nfv[i]=0;fv[i]=(int*)NULL;}
        ft=(int*)realloc((void*)ft,mF*sizeof(int));
        if(ft==(int*)NULL){printf("Error allocating %d bytes for ft\n",mF*sizeof(int));fflush(stdout);abort();}
       }

      nfv[iF]=ne;
      fv[iF]=(int*)malloc(nfv[iF]*sizeof(int));
      if(fv[iF]==(int*)NULL){printf("Error allocating %d bytes for fv[%d]\n",nfv[iF]*sizeof(int),iF);fflush(stdout);abort();}

      ft[iF]=0;
      (fv[iF])[0]=ev[0+2*iE0];
      (fv[iF])[1]=ev[1+2*iE0];
      if(verbose){printf("  Face %d (%d vertices)\n",iF,ne);fflush(stdout);}
      if(verbose){printf("     vertex %d = %d\n",0,(fv[iF])[0]);fflush(stdout);}
      if(verbose){printf("     vertex %d = %d\n",1,(fv[iF])[1]);fflush(stdout);}

      np=2;

      ie=iE0;
      je=iE0;
      trip=0;
      while(np<ne&&trip<2)
       {
        if(ev[0+2*ie]==(fv[iF])[np-1]&&je!=ie)
         {
          (fv[iF])[np]=ev[1+2*ie];
          if(verbose){printf("     vertex %d = %d\n",np,(fv[iF])[np]);fflush(stdout);}
          np++;
          je=ie;
          trip=0;
         }else if(ev[1+2*ie]==(fv[iF])[np-1]&&je!=ie)
         {
          (fv[iF])[np]=ev[0+2*ie];
          if(verbose){printf("     vertex %d = %d\n",np,(fv[iF])[np]);fflush(stdout);}
          np++;
          je=ie;
          trip=0;
         }
        ie++;
        if(ie>iE0+ne-1){ie=iE0;trip++;}
       }
      iF++;
     }
   }

  *piV=iV;
  *piV0=iV0;
  *pmV=mV;
  *pv=v;
  *pvt=vt;

  *piE=iE;
  *pmE=mE;
  *pev=ev;
  *pet=et;

  *piF=iF;
  *pmF=mF;
  *pnfv=nfv;
  *pfv=fv;
  *pft=ft;

  return;
 }

void PFDualPolytopes(FILE *fid,int m, int k, int *piV,int *piV0,int *pmV,double **pv,int **pvt,int *piE,int *pmE,int **pev,int **pet, int *piF,int *pmF,int **pnfv,int ***pfv, int **pft, MFErrorHandler e)
 {
  int mV,iV,iV0,nV;
  double *v;
  int *vt;
  int mE,iE;
  int *ev;
  int *et;
  int mF,iF;
  int *nfv;
  int **fv;
  int *ft;
  int i,j,l,n;
  int ipoly,nv,ne,nf;
  int verbose=0;
  int iv,ie,iE0,jv,je,np,v0,v1;
  int trip;
  int sing,bnd;
  int ic,nc;
  int mVertexIndex;
  int *nVertexIndex;
  int **vertexIndex;
  int face;
  char c;
  double R;

  iV=*piV;
  iV0=*piV0;
  nV=*piV;
  mV=*pmV;
  v=*pv;
  vt=*pvt;
  iE=*piE;
  mE=*pmE;
  ev=*pev;
  et=*pet;
  mF=*pmF;
  iF=*piF;
  nfv=*pnfv;
  fv=*pfv;
  ft=*pft;

  mVertexIndex=0;
  nVertexIndex=NULL;
  vertexIndex=NULL;

  while(!feof(fid))
   {
    fscanf(fid,"Polyhedron %d",&ipoly);
    c=fgetc(fid);
    ungetc(c,fid);
    if(c=='R')
     {
      fscanf(fid,"R=%lf, %d vertices, %d edges, %d faces,  boundary %d, singular %d\n",&R,&nv,&ne,&nf,&bnd,&sing);
     }else{
      fscanf(fid,"%d vertices, %d edges, %d faces,  boundary %d, singular %d\n",&nv,&ne,&nf,&bnd,&sing);
     }

    if(verbose){printf("Polyhedron %d, %d vertices, %d edges, %d faces\n",ipoly,nv,ne,nf);fflush(stdout);}

    if(nv>mVertexIndex)
     {
      mVertexIndex+=20;
      nVertexIndex=(int*)realloc(nVertexIndex,mVertexIndex*sizeof(int));
      vertexIndex=(int**)realloc(vertexIndex,mVertexIndex*sizeof(int*));
     }

/*  Just store the last "vertex" (i.e. the center) */

    iV=ipoly;
    if(iV>=nV)nV=iV+1;
    if(nV>=mV)
     {
      mV=nV+100;
      v=(double*)realloc((void*)v,m*mV*sizeof(double));
      vt=(int*)realloc((void*)vt,mV*sizeof(int));
     }

    if(nv-1+iF>=mF)
     {
      mF+=nv+100;
      nfv=(int*)realloc((void*)nfv,mF*sizeof(int));
      fv=(int**)realloc((void*)fv,mF*sizeof(int*));
      for(iv=0;iv<nv-1;iv++)fv[iv]=(int*)malloc((nVertexIndex[iv]+1)*sizeof(int));

      ft=(int*)realloc((void*)ft,mF*sizeof(int));
     }

/* Example:

Polyhedron 6, 7 vertices, 6 edges, 6 faces, boundary 0, singular 0
Vertex 0 (0.196000,0.969965,-0.155275), 2 [105,113]
Vertex 1 (0.084337,0.989945,-0.127861), 2 [34,53]
Vertex 2 (0.127175,0.989930,-0.085236), 2 [34,45]
Vertex 3 (0.096403,0.980696,-0.180078), 2 [53,123]
Vertex 4 (0.180886,0.980284,-0.098668), 2 [45,113]
Vertex 5 (0.154125,0.970009,-0.196737), 2 [105,123]
Vertex 6 (0.140710,0.979950,-0.141065), 0 [ ]
Edge 0 (0,4), 1 [113]
Edge 1 (0,5), 1 [105]
Edge 2 (1,2), 1 [34]
Edge 3 (1,3), 1 [53]
Edge 4 (2,4), 1 [45]
Edge 5 (3,5), 1 [123]
Face 105 neighbor 15
Face 53 neighbor 8
Face 34 neighbor 1
Face 45 neighbor 7
Face 113 neighbor 16
Face 123 neighbor 17

*/

    for(iv=0;iv<nv;iv++)
     {
      fscanf(fid,"Vertex %*d (%lf",&(v[0+m*iV]));
      for(i=1;i<m;i++)fscanf(fid,",%lf",&(v[i+m*iV]));
/*    fscanf(fid,"), %*d [%*[ 0-9+-.,]]\n"); */

      fscanf(fid,"), %d [",&(nVertexIndex[iv]));

      vertexIndex[iv]=NULL;
      if(nVertexIndex[iv]>0)
       {
        vertexIndex[iv]=(int*)malloc(nVertexIndex[iv]*sizeof(int));
        ic=0;
        fscanf(fid,"%d",&((vertexIndex[iv])[ic]));
        for(ic=1;ic<nVertexIndex[iv];ic++)
         {
          fscanf(fid,",%d",&((vertexIndex[iv])[ic]));
         }
        fscanf(fid,"]\n");
       }else
        fscanf(fid,"%*[ 0-9+-.,]]\n");

      vt[iV]=0;
      if(verbose)
       {
        printf("  Vertex %d (%lf",iv,v[0+m*iV]);
        for(i=1;i<m;i++)printf(",%lf",v[i+m*iV]);
        printf(") %d [",nVertexIndex[iv]);
        for(i=0;i<nVertexIndex[iv];i++)
         {
          if(i>0)printf(",");
          printf("%d",vertexIndex[iv][i]);
         }
        printf("]\n");fflush(stdout);
       }
     }

    if(nf+iE>=mE)
     {
      mE+=nf+100;
      ev=(int*)realloc((void*)ev,2*mE*sizeof(int));
      et=(int*)realloc((void*)et,mE*sizeof(int));
     }

    for(iv=0;iv<nv-1;iv++)
     {
      nfv[iF+iv]=nVertexIndex[iv]+1;
      ft[iF+iv]=0;
      fv[iF+iv]=(int*)malloc(nfv[iF+iv]*sizeof(int));

      fv[iF+iv][0]=iV;
     }

    for(ie=0;ie<ne;ie++)
     {
      if(k>1)
        fscanf(fid,"Edge %*d (%d,%d), %*d [%*[ 0-9+-.,]]\n",&v0,&v1);
      else if(k==1)
        fscanf(fid,"Edge %*d (%d,%d), %*d []\n",&v0,&v1);
       else
        fscanf(fid,"Edge %*d (%d,%d), %*d [%*[ 0-9+-.,]]\n",&v0,&v1);
     }

    for(i=0;i<nf;i++)
     {
      fscanf(fid,"Face %d neighbor %d\n",&face,&iV0);
      if(verbose){printf("Face %d neighbor %d\n",face,iV0);fflush(stdout);}
      if(iV>iV0 && iV0>-1)
       {
        ev[0+4*iE]=iV0;
        ev[1+4*iE]=iV;    
        et[iE]=0;  
        if(verbose){printf("  Dual Edge %d (%d,%d)\n",iE,ev[0+4*iE],ev[1+4*iE]);fflush(stdout);}
        iE++;
        for(iv=0;iv<nv-1;iv++)
         {
          fv[iF+iv][0]=iV;
          for(jv=0;jv<nVertexIndex[iv];jv++)
           {
            if(vertexIndex[iv][jv]==face)fv[iF+iv][jv+1]=iV0;
           }
         }
       }
     }

    for(iv=0;iv<nv-1;iv++)fv[iF+iv][0]=iV;

    if(verbose)
     {
      for(iv=0;iv<nv-1;iv++)
       {
        printf("   Vertex %d",iv);
        printf("    index [");
        for(jv=0;jv<nVertexIndex[iv];jv++)
         {
          if(jv>0)printf(",");
          printf("%d",vertexIndex[iv][jv]);
         }
        printf("]  neighbors [");
        for(jv=0;jv<nfv[iF+iv];jv++)
         {
          if(jv>0)printf(",");
          printf("%d",fv[iv][jv]);
         }
        printf("]\n");fflush(stdout);
       }
     }

    for(iv=0;iv<nv;iv++)
      if(vertexIndex[iv]!=NULL)free(vertexIndex[iv]);

    iF=iF+nv-1;
   }

  free(nVertexIndex);
  free(vertexIndex);

  *piV=nV;
  *piV0=iV0;
  *pmV=mV;
  *pv=v;
  *pvt=vt;

  *piE=iE;
  *pmE=mE;
  *pev=ev;
  *pet=et;

  *piF=iF;
  *pmF=mF;
  *pnfv=nfv;
  *pfv=fv;
  *pft=ft;

  return;
 }

void MFPlotfileDualToDX(FILE *fid,char *name, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPlotfileDualToDX"};
  int mV,iV,iV0;
  double *v=(double*)NULL;
  int *vt=(int*)NULL;
  int mE,iE;
  int *ev=(int*)NULL;
  int *et=(int*)NULL;
  int mF,iF;
  int *nfv=(int*)NULL;
  int **fv=(int**)NULL;
  int *ft=(int*)NULL;
  int i,j,l,n,m,k;
  int verbose=0;

  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}

  iV=0;
  iV0=0;
  mV=0;
  v=NULL;
  vt=NULL;

  iE=0;
  mE=0;
  ev=NULL;
  et=NULL;

  iF=0;
  mF=0;
  fv=NULL;
  ft=NULL;
  nfv=NULL;

  fscanf(fid,"Dimension of vertices, %d\n",&m);
  fscanf(fid,"Dimension of manifold, %d\n",&k);

  if(verbose){printf("Add in Edges\n");fflush(stdout);}
  PFDualPolytopes(fid,m,k,&iV,&iV0,&mV,&v,&vt,&iE,&mE,&ev,&et,&iF,&mF,&nfv,&fv,&ft,e);
  if(verbose){printf("  There are now %d pts, %d edges %d faces\n",iV,iE,iF);fflush(stdout);}
  if(verbose){printf("Dump it to %s\n",name);fflush(stdout);}
  MFMeshToDXFile2(name,iV,m,v,vt,iE,ev,et,iF,nfv,fv,ft,e);
  if(v!=(double*)NULL)free(v);
  if(nfv!=(int*)NULL)free(nfv);
  if(fv!=(int**)NULL)
   {
    for(i=0;i<iF;i++)
     {
      if(fv[i]!=(int*)NULL)free(fv[i]);
     }
    free(fv);
    free(ft);
    free(nfv);
   }
  return;
 }

#if 0

#define NumberToAllocateToVertices 10
#define NumberToAllocateToPolys 100

void MFPlotfileDualToEasyMesh(FILE *fid,char *name, MFErrorHandler e)
 {
  static char RoutineName[]={"MFPlotfileDualToEasyMesh"};
  int iV,iV0;
  int mE,iE;
  int *ev=(int*)NULL;
  int *eR=(int*)NULL;
  int mF,iF;
  int *nfv=(int*)NULL;
  int **fv=(int**)NULL;
  int i,j,l,n,m,k;
  int ipoly,nv,ne,nf;
  int verbose=0;
  int iv,ie,iE0,je,np,v0,v1;
  int trip;

  int mpoly;
  int *polyNVertices=NULL;
  double **polyVertex=NULL;
  int ***polyVertexIndices=NULL;
  int **polyNVertexIndices=NULL;
  int *polyNEdges=NULL;
  int **polyEdge=NULL;
  int mface;
  int *faces=NULL;
  char c;

  double v;

  fscanf(fid,"Dimension of vertices, %d\n",&m);
  fscanf(fid,"Dimension of manifold, %d\n",&k);

  if(k!=2)
   {
    fprintf("Error: %s, k(%2) must be 2\n",RoutineName,k);fflush(stdout);
    return;
   }

  mpoly=0;
  mface=0;

  while(!feof(fid))
   {

/* Input Example:
 *
 *     Polyhedron 6, 7 vertices, 6 edges, 6 faces, boundary 0, singular 0
 *     Vertex 0 (0.196000,0.969965,-0.155275), 2 [105,113]
 *     Vertex 1 (0.084337,0.989945,-0.127861), 2 [34,53]
 *     Vertex 2 (0.127175,0.989930,-0.085236), 2 [34,45]
 *     Vertex 3 (0.096403,0.980696,-0.180078), 2 [53,123]
 *     Vertex 4 (0.180886,0.980284,-0.098668), 2 [45,113]
 *     Vertex 5 (0.154125,0.970009,-0.196737), 2 [105,123]
 *     Vertex 6 (0.140710,0.979950,-0.141065), 0 [ ]
 *     Edge 0 (0,4), 1 [113]
 *     Edge 1 (0,5), 1 [105]
 *     Edge 2 (1,2), 1 [34]
 *     Edge 3 (1,3), 1 [53]
 *     Edge 4 (2,4), 1 [45]
 *     Edge 5 (3,5), 1 [123]
 *     Face 105 neighbor 15
 *     Face 53 neighbor 8
 *     Face 34 neighbor 1
 *     Face 45 neighbor 7
 *     Face 113 neighbor 16
 *     Face 123 neighbor 17
 */

    fscanf(fid,"Polyhedron %d, %d vertices, %d edges, %d faces,  boundary %d, singular %d\n",&ipoly,&nv,&ne,&nf,&bnd,&sing);
    c=fgetc(fid);
    ungetc(c,fid);
    if(c=='R')
     {
      fscanf(fid,"R=%lf, %d vertices, %d edges, %d faces,  boundary %*d, singular %*d\n",&R,&nv,&ne,&nf);
     }else{
      fscanf(fid,"%d vertices, %d edges, %d faces,  boundary %*d, singular %*d\n",&nv,&ne,&nf);
     }

    if(verbose){printf("Polyhedron %d, %d vertices, %d edges, %d faces\n",ipoly,nv,ne,nf);fflush(stdout);}

    if(ipoly>=mpoly)
     {
      mpoly=ipoly+NumberToAllocateToPoly;

      polyNVertices=(int*)realloc(polyNVertices,mpoly*sizeof(int));
      polyVertex=(double**)realloc(polyVertex,mpoly*sizeof(double*));
      polyNVertexIndices=(int*)realloc(polyNVertexIndices,mpoly*sizeof(int));
      polyVertexIndices=(int***)realloc(polyVertexIndices,mpoly*sizeof(int**));
      polyNEdges=(int*)realloc(polyNEdges,mpoly*sizeof(int));
      polyEdge=(int**)realloc(polyEdge,mpoly*sizeof(int*));
      polyNFaces=(int*)realloc(polyNFaces,mpoly*sizeof(int));
      polyFace=(int**)realloc(polyFace,mpoly*sizeof(int*));
     }

    polyNVertices[ipoly]=nv;
    polyVertex[ipoly]=(double*)malloc(m*sizeof(double));
    polyVertexNIndices[ipoly]=(int*)malloc(nv*sizeof(int));
    polyVertexIndices[ipoly]=(int**)malloc(nv*sizeof(int*));

    polyNEdges[ipoly]=ne;
    polyEdge[ipoly]=(int*)malloc(4*ne*sizeof(int));

    polyNFaces[ipoly]=nf;
    polyFace[ipoly]=(int*)malloc(2*nf*sizeof(int));

/*  Just store the last "vertex" (i.e. the center), but keep the face indices for the others */

    iV=ipoly;

    for(iv=0;iv<nv-1;iv++)
     {
      fscanf(fid,"Vertex %*d (%lf",&v);
      for(i=1;i<m;i++)fscanf(fid,",%lf",&v);

      fscanf(fid,"), %d [",&(polyNVertexIndices[ipoly][iv]));
      polyVertexIndices[ipoly][iv]=(int*)malloc(polyNVertexIndices[ipoly,iv]*sizeof(int));

      if(nVertexIndex[ipoly][iv]>0)
       {
        ic=0;
        fscanf(fid,"%d",&((polyVertexIndices[ipoly][iv])[ic]));
        for(ic=1;ic<nVertexIndex[ipoly][iv];ic++)
         {
          fscanf(fid,"%d",&((PolyVertexIndices[ipoly][iv])[ic]));
         }
        fscanf(fid,"]\n");
       }else
        fscanf(fid,"%*[ 0-9+-.,]]\n");
     }

/* Edges are next, keep both the endpoints and the faces. This is a quadedge, though without the orientation and rotation */

    polyEdge[ipoly]=(int**)malloc(ne*sizeof(int*));
    polyEdge[ipoly]=(int**)malloc(ne*sizeof(int*));
    for(ie=0;ie<ne;ie++)
     {
      polyEdge[ipoly][ie]=(int*)malloc(4*sizeof(int));
      fscanf(fid,"Edge %*d (%d",&(polyEdge[ipoly][ie][0]));
      fscanf(fid,",%d)",&(polyEdgep[ipoly][ie][1]));
      fscanf(fid,"%*d [%d]\n",&(polyEdge[ipoly][ie][2]));
     }

/* Finally the faces. */

    for(i=0;i<nf;i++)
     {
      fscanf(fid,"Face %d neighbor %d\n",&face,&iV0);
      if(face>=mface)
       {
        mface=face+20;
        faces=(int*)realloc(faces,mface*sizeof(int));
       }
      faces[2*face  ]=iV;
      faces[2*face+1]=iV0;
      if(verbose){printf("Face %d neighbor %d\n",face,iV0);fflush(stdout);}
     }

    for(ie=0;ie<ne;ie++)
     {
      i=polyEdge[ipoly][ie][2];
      polyEdge[ipoly][ie][2]=faces[2*i  ];
      polyEdge[ipoly][ie][3]=faces[2*i+1];
     }

    for(iv=0;iv<nv;iv++)
     {
      i=polyVertexIndex[ipoly][iv][0];
      polyVertexIndex[ipoly][iv][0]=faces[2*i+1];
      i=polyVertexIndex[ipoly][iv][1];
      polyVertexIndex[ipoly][iv][1]=faces[2*i+1];
     }

   }

/* Build dual data structure */

/* Quadedge data structure
 *
 * Each quadedge is (v0,v1,f0,f1,next.orig=<qe,r,f>,next.dest=<qe,r,f>,next.left=<qe,r,f>,next.right=<qe,r,f>)
 *
 */

  nv=0;
  for(ipoly=0;ipoly<npoly;ipoly++)
   for(iv=0;iv<polyNVertices[ipoly];iv++)
    if(ipoly<polyVertexIndex[ipoly][iv][0] && Ipoly<polyVertexIndex[ipoly][iv][1])nv++;

  ne=0;
  for(ipoly=0;ipoly<npoly;ipoly++)
   for(ie=0;ie<polyNEdges[ipoly];ie++)
    if(polyEdge[ipoly][ie][2]<polyEdge[ipoly][ie][3])ne++;

  nf=0;
  for(ipoly=0;ipoly<npoly;ipoly++)
   for(iv=0;iv<polyNVertices[ipoly];iv++)
     if(ipoly<polyVertexIndex[ipoly][iv][0]&&ipoly<polyVertexIndex[ipoly][iv][1])nf++;

/* Now write three files
 *   nodes: name.n      
 *   sides: name.s
 *   elements: name.e
 */

/* nodes: name.n
 *
 *     # first line:       <number of nodes>
 *     # following lines:  <node number:> <x> <y> <marker>
 */

  strcpy(fullname,name);
  strcat(fullname,".n");
  fout=fopen(fullname,"w");
  fprintf("%d\n",nv);
  for(ipoly=0;ipoly<npoly;ipoly++)
   {
    fprintf(fout,"%d",ipoly);
    for(j=0;j<m;j++)
     {
      if(j>0)fprintf(fout," ");
      fprintf(fout,"%lf",polyVertex[ipoly][j]);
     }
    fprintf(fout,"\n");
   }
  fclose(fout);

/* edges (sides): name.s
 *
 *    # first line:       <number of sides>
 *    # following lines:  <side number:> <c> <d> <ea> <eb> <marker> 
 */

  strcpy(fullname,name);
  strcat(fullname,".s");
  fout=fopen(fullname,"w");

/*
 * Each edge appears twice, the test eliminates one
 */

  fprintf("%d\n",ne);

  ne=0;
  for(ipoly=0;ipoly<npoly;ipoly++)
   {
    for(ie=0;ie<polyNEdges[ipoly];ie++)
     {
      if(polyEdge[ipoly][ie][0]<polyEdge[ipoly][ie][1])
       {
        fprintf(fout,"%d %d %d %d %d\n",ne,polyEdge[ipoly][ie][0],polyEdge[ipoly][ie][1],polyEdge[ipoly][ie][2],polyEdge[ipoly][ie][3]);
        ne++;
       }
     }
   }
  fclose(fout);

/* elements: name.e
 *
 *     # first line        <number of elements>
 *     # following lines:  <element number:> <i> <j> <k> <ei> <ej> <ek> <si> <sj> <sk> <xV> <yV> <marker>
 */

  strcpy(fullname,name);
  strcat(fullname,".e");
  fout=fopen(fullname,"w");

  fprintf("%d\n",nf);

  nf=0;
  for(ipoly=0;ipoly<npoly;ipoly++)
   {
    for(iv=0;iv<polyNVertices[ipoly];iv++)
     {
      if(ipoly<polyVertexIndex[ipoly][iv][0]&&ipoly<polyVertexIndex[ipoly][iv][1])
       {
        fprintf(fout,"%d %d %d %d %d %d %d %lf %lf\n",nf,
                                        ipoly,polyVertexIndex[ipoly][iv][0],polyVertexIndex[ipoly][iv][1],
                                        ei,ej,ek,
                                        polyVertex[ipoly][0],polyVertex[ipoly][1]);
        nf++;
       }
     }
   }
  fclose(fout);

  free(nVertexIndex);
  free(vertexIndex);

  free(nVEdgeIndex);
  free(edgeIndex);

  if(v!=(double*)NULL)free(v);
  if(nfv!=(int*)NULL)free(nfv);
  if(fv!=(int**)NULL)
   {
    for(i=0;i<iF;i++)
     {
      if(fv[i]!=(int*)NULL)free(fv[i]);
     }
    free(fv);
    free(nfv);
   }
  return;
 }
#endif

#ifdef __cplusplus
}
#endif
