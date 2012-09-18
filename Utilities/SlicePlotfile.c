#include <stdio.h>
#include <MFAtlas.h>
#include <MFDraw.h>
#include <sh.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

struct MFModifiedPolytopeSt{
                            int k;                /* Dimension of the vertices */
                            int n;                /* Number of vertices */
                            double R;             /* Radius */
                            int m;                /* Space for vertices */
                            double *v;            /* vertices */
                            int *nIndices;        /* Number of indices for each vertex */
                            int *mIndices;        /* Space for indices for each vertex */
                            int **indices;        /* Indices of each vertex */
                            int *mark;            /* Mark on each vertex */

                            int nFaces;
                            int mFaces;
                            int *face;          /* Index of face */
                            int *nFaceV;        /* Number of vertices on each face */
                            MFNVector *faceN;   /* Normal of face */
                            double *faceO;      /* Origin of face */
                            int dimV;
                           };
typedef struct MFModifiedPolytopeSt* MFModifiedPolytope;

static char MFAtlasErrorMsg[256]="";

void MFWriteModifiedPolytopeToPlotFile(int i, FILE *fid, MFModifiedPolytope P, MFErrorHandler e);
MFModifiedPolytope MFReadPolytopeFromPlotfile(FILE *fid, int k, int dimV, MFErrorHandler e);
void MFFreeModifiedPolytope(MFModifiedPolytope P, MFErrorHandler e);
MFModifiedPolytope MFCreateSectionOfPolytope(MFModifiedPolytope P, int index, MFNVector nrm, double on, MFErrorHandler e);

int main(int argc, char *argv[])
 {
  static char RoutineName[]={"main"};
  FILE *fin;
  FILE *fout;
  char file[4096];
  int n,k;
  MFModifiedPolytope P,Q;
  double on;
  double  d;
  MFNVector nrm;
  MFErrorHandler e;
  int i;

  if(argc<2){printf("Usage %s infilename o n[0] n[1] n[2] ...\n",argv[0]);fflush(stdout);return 8;}

  strcpy(file,argv[1]);
  strcat(file,".plotfile");
  fin=fopen(file,"r");
  if(fin==(FILE*)NULL){printf("Error opening file %s\n",file);fflush(stdout);return 8;}

  strcpy(file,"sliced");
  strcat(file,argv[1]);
  strcat(file,".plotfile");
  fout=fopen(file,"w");
  if(fout==(FILE*)NULL){printf("Error opening file %s\n",file);fflush(stdout);return 8;}

  fscanf(fin,"Dimension of vertices, %d\n",&n);
  fscanf(fin,"Dimension of manifold, %d\n",&k);

  sscanf(argv[1],"%lf",&on);
  printf("  origin=%lf\n",on);fflush(stdout);
  nrm=MFCreateNVector(n,e);
  for(i=0;i<n;i++)
   {
    d=atoi(argv[3+i]);
    printf("  nrm[%d]=%lf\n",i,d);fflush(stdout);
    MFNVSetC(nrm,i,d,e);
   }

  e=MFCreateErrorHandler();

  fprintf(fout,"Dimension of vertices, %d\n",n);
  fprintf(fout,"Dimension of manifold, %d\n",k);

  i=0;
  while(!feof(fin))
   {
    P=MFReadPolytopeFromPlotfile(fin,k,n,e);
    Q=MFCreateSectionOfPolytope(P,-1,nrm,on,e);
    MFFreeModifiedPolytope(P,e);
    MFWriteModifiedPolytopeToPlotFile(i,fout,Q,e);
    MFFreeModifiedPolytope(Q,e);
    i++;
   }

  fclose(fin);
  fclose(fout);
  MFFreeErrorHandler(e);

  return 0;
 }

MFModifiedPolytope MFReadPolytopeFromPlotfile(FILE *fid, int k, int dimV, MFErrorHandler e)
 {
  static char RoutineName[]={"MFReadPolytopeFromPlotfile"};
  MFModifiedPolytope P;
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

  if(verbose){printf("     Read Polytope\n");fflush(stdout);}
  P=(MFModifiedPolytope)malloc(sizeof(struct MFModifiedPolytopeSt));
  P->k=k;
  P->dimV=dimV;

  fscanf(fid,"Polyhedron %d, ",&ipoly);
  if(verbose){printf("          Polytope %d\n",ipoly);fflush(stdout);}
  c=fgetc(fid);
  ungetc(c,fid);
  if(c=='R')
   {
    fscanf(fid,"R=%lf, %d vertices, %d edges, %d faces,  boundary %d, singular %d\n",&R,&nv,&ne,&nf,&bnd,&sing);
   }else{
    fscanf(fid,"%d vertices, %d edges, %d faces,  boundary %d, singular %d\n",&nv,&ne,&nf,&bnd,&sing);
    R=0;
   }
  fscanf(fid,"  %d vertices, %d edges, %d faces,  boundary %d, singular %d\n",nv,ne,nf,bnd,sing);

  P->n=nv;
  P->m=nv;
  P->v=(double*)malloc(P->n*dimV*sizeof(double));
  P->nIndices=(int*)malloc(P->n*sizeof(int));
  P->mIndices=(int*)malloc(P->n*sizeof(int));
  P->indices=(int**)malloc(P->n*sizeof(int*));
  P->mark=(int*)malloc(P->n*sizeof(int));
  for(i=0;i<P->n;i++)
   {
    P->nIndices[i]=0;
    P->mIndices[i]=0;
    P->indices[i]=(int*)NULL;
    P->mark[i]=0;
   }

  if(verbose){printf("Polyhedron %d, %d vertices, %d edges, %d faces\n",ipoly,nv,ne,nf);fflush(stdout);}

  for(iv=0;iv<nv;iv++)
   {
    fscanf(fid,"Vertex %d (%lf",&iV,&(P->v[0+P->dimV*iv]));
    if(verbose){printf("Vertex %d (%lf",iV,(P->v[0+P->dimV*iv]));fflush(stdout);}
    for(i=1;i<dimV;i++){fscanf(fid,",%lf",&(P->v[i+P->dimV*iv]));if(verbose){printf(",%lf",(P->v[i+P->dimV*iv]));fflush(stdout);}}
    fscanf(fid,"), %d ",&(P->nIndices[iv]));
    if(verbose){printf("), %d ",(P->nIndices[iv]));fflush(stdout);}
    P->mIndices[iv]=P->nIndices[iv];
    if(P->nIndices[iv]!=0)
     {
      P->indices[iv]=(int*)malloc((P->nIndices[i])*sizeof(int));

      fscanf(fid,"[%d",&(P->indices[iv][0]));
      if(verbose){printf(    "[%d", (P->indices[iv][0]));fflush(stdout);}

      for(i=1;i<P->nIndices[iv];i++)
       {
        fscanf(fid,",%d",&(P->indices[iv][i]));
        if(verbose){printf(    ",%d", (P->indices[iv][i]));fflush(stdout);}
       }

      fscanf(fid,"]\n");
      if(verbose){printf(    "]\n");fflush(stdout);}
     }else{
      fscanf(fid," [ ]\n");
      if(verbose){printf(    "\n");fflush(stdout);}
     }

    if(verbose)
     {
      printf("  Vertex %d (%lf",iv,P->v[0+P->dimV*iV]);
      for(i=1;i<P->dimV;i++)printf(",%lf",P->v[i+P->dimV*iV]);
      printf(")\n");fflush(stdout);
     }
   }


  for(ie=0;ie<ne;ie++)
   {
    if(k>1)
     {
      fscanf(fid,"Edge %d (%d,%d), %*d [%*[ 0-9+-.,]]\n",&iE,&v0,&v1);
      if(verbose){printf(    "Edge %d (%d,%d)\n",iE,v0,v1);fflush(stdout);}
     }else if(k==1)
     {
      fscanf(fid,"Edge %d (%d,%d), %*d []\n",&iE,&v0,&v1);
      if(verbose){printf(    "Edge %d (%d,%d)\n",iE,v0,v1);fflush(stdout);}
     }else{
      fscanf(fid,"Edge %d (%d,%d), %*d [%*[ 0-9+-.,]]\n",&iE,&v0,&v1);
      if(verbose){printf(    "Edge %d (%d,%d)\n",iE,v0,v1);fflush(stdout);}
     }
    if(verbose){printf("  Edge %d (%d,%d)\n",ie,v0,v1);fflush(stdout);}
   }

  P->nFaces=nf;
  P->mFaces=nf;
  P->face  =(int*)malloc(P->mFaces*sizeof(int));
  P->nFaceV=(int*)malloc(P->mFaces*sizeof(int));
  P->faceN =(MFNVector*)malloc(P->mFaces*sizeof(MFNVector));
  P->faceO =(double*)malloc(P->mFaces*sizeof(double));

  for(i=0;i<nf;i++)
   {
    fscanf(fid,"Face %d neighbor %*d\n",&(P->face[i]));
    if(verbose){printf(    "Face %d neighbor\n",(P->face[i]));fflush(stdout);}
    P->nFaceV[i]=0;
    P->faceN [i]=(MFNVector)NULL;
    P->faceO [i]=0.;
   }
  if(verbose){printf("done\n");fflush(stdout);}

  return P;
 }

int MFModifiedPolytopeIntersectIndexSets(MFModifiedPolytope P,int v1,int v2,int *inter, MFErrorHandler e)
 {
  static char RoutineName[]={"MFModifiedPolytopeIntersectIndexSets"};

/* Returns the number of entries two index sets have in common  */
/*   Assumes that the index sets are sorted in increasing order */

  int n;
  int n1,n2;
  int *inter1;
  int *inter2;
  int i1,i2;

  n1=P->nIndices[v1];
  inter1=P->indices[v1];
  n2=P->nIndices[v2];
  inter2=P->indices[v2];

  n=0;
  i1=0;i2=0;
  while(i1<n1&&i2<n2)
   {
    if(inter1[i1]==inter2[i2])
     {
      inter[n]=inter1[i1];
      n++;
      i1++;i2++;
     }else if(inter1[i1]<inter2[i2])
     {
      i1++;
     }else if(inter1[i1]>inter2[i2])
     {
      i2++;
     }
   }

  return n;
 }

MFModifiedPolytope MFCreateSectionOfPolytope(MFModifiedPolytope P, int index, MFNVector nrm, double on, MFErrorHandler e)
 {
  static char RoutineName[]={"MFCreateSectionOfPolytope"};
  MFModifiedPolytope Q;
  int i,j;
  int k;
  int n,m;
  double *d;
  double t;
  int *inter;

  Q=(MFModifiedPolytope)malloc(sizeof(struct MFModifiedPolytopeSt));
  inter=(int*)malloc(P->k*sizeof(int));

  k=P->k;

  Q->k=P->k;
  Q->R=P->R;

  d=(double*)malloc(P->n*sizeof(double));

  t=fabs(on);
  for(j=0;j<P->dimV;j++)if(fabs(MFNV_C(nrm,j,e))>t)t=fabs(MFNV_C(nrm,j,e));

/* Mark the vertices of P. */

  for(i=0;i<P->n;i++)
   {

    d[i]=-on;
    for(j=0;j<P->dimV;j++)d[i]+=P->v[j+P->dimV*i]*MFNV_C(nrm,j,e);
    d[i]=d[i]/t;

/*  printf("   (%lf",P->v[0+P->dimV*i]);fflush(stdout);
    for(j=1;j<P->dimV;j++){printf(",%lf",P->v[j+P->dimV*i]);fflush(stdout);}

    printf(").(%lf",MFNV_C(nrm,0,e)/t);fflush(stdout);
    for(j=1;j<P->dimV;j++){printf(",%lf",MFNV_C(nrm,j,e)/t);fflush(stdout);}
    printf(") - %lf=%lf\n",on,d[i]);fflush(stdout);
 */
   }

/* Count the number of edges that cross (the number of vertices in the Q). */

  Q->n=0;
  for(i=0;i<P->n;i++)
   {
    for(j=0;j<i;j++)
     {
      if(MFModifiedPolytopeIntersectIndexSets(P,i,j,inter,e)==P->k-1 && d[i]*d[j]<=0.)Q->n++;
     }
   }

/* Allocate space for the new vertices. */

  Q->m=Q->n;
  Q->dimV=P->dimV;
  Q->v=(double*)malloc(Q->n*Q->dimV*sizeof(double));
  Q->nIndices=(int*)malloc(Q->n*sizeof(int));
  Q->mIndices=(int*)malloc(Q->n*sizeof(int));
  Q->indices=(int**)malloc(Q->n*sizeof(int*));
  Q->mark=(int*)malloc(Q->n*sizeof(int));

  for(i=0;i<Q->n;i++)
   {
    Q->nIndices[i]=Q->k;
    Q->mIndices[i]=Q->k;
    Q->indices[i]=(int*)malloc(Q->k*sizeof(int));
    Q->mark[i]=0;
   }

/* Fill in the values. */
 
  n=0;
  for(i=0;i<P->n;i++)
   {
    for(j=0;j<i;j++)
     {
      if(MFModifiedPolytopeIntersectIndexSets(P,i,j,inter,e)==P->k-1 && d[i]*d[j]<=0.)
       {
        Q->indices[n][0]=index;
        for(m=1;m<Q->nIndices[n];m++)Q->indices[n][m]=inter[m-1];

        t=d[i]/(d[i]-d[j]);

        for(m=0;m<Q->dimV;m++)
          Q->v[m+Q->dimV*n]=(1-t)*P->v[m+Q->dimV*i]+t*P->v[m+Q->dimV*j];

        n++;
       }
     }
   }

  Q->nFaces=P->nFaces+1;
  Q->mFaces=Q->nFaces;
  Q->face=(int*)malloc(Q->mFaces*sizeof(int));
  Q->nFaceV=(int*)malloc(Q->mFaces*sizeof(int));
  Q->faceO=(double*)malloc(Q->mFaces*sizeof(double));
  Q->faceN=(MFNVector*)malloc(Q->mFaces*sizeof(MFNVector));
  for(i=0;i<P->nFaces;i++)
   {
    Q->face[i]=P->face[i];
    Q->faceN[i]=P->faceN[i];
    if(Q->faceN[i]!=(MFNVector)NULL)MFRefNVector(Q->faceN[i],e);
    Q->faceO[i]=P->faceO[i];
    Q->nFaceV[i]=P->nFaceV[i];
   }
  Q->face[P->nFaces]=index;
  Q->faceN[P->nFaces]=nrm; MFRefNVector(nrm,e);
  Q->faceO[P->nFaces]=on;
  Q->nFaceV[P->nFaces]=0;

  free(d);
  free(inter);

  return Q;
 }

void MFWriteModifiedPolytopeToPlotFile(int chart, FILE *fid, MFModifiedPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFWriteModifiedPolytopeToPlotFile"};
  int d;
  double *x;
  MFKVector s;
  double *ps;
  MFNVector u;
  int i,j,k,l;
  int n,m,ne,ie;
  int nv;
  int c;
  int *index;
  int *indx;
  int verbose=0;
  MFNKMatrix Phi;
  MFNVector col;
  int bnd,sing;
  double *dv;

  if(fid==NULL){printf("%s - fid is NULL %s(%d)\n",RoutineName,__LINE__,__FILE__);fflush(stdout);abort();}
  if(P->n<P->k)return;

  indx=(int*)malloc(P->k*sizeof(int));

#ifndef MFNOSAFETYNET
  if(indx==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Out of memory trying to allocate %d bytes",P->k*sizeof(int));
    MFSetError(e,12,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  ne=0;
  for(i=0;i<P->n-1;i++)
    for(j=i+1;j<P->n;j++)
     {
      k=MFModifiedPolytopeIntersectIndexSets(P,i,j,indx,e);
      if(k==P->k-1)ne++;
     }

  m=P->nFaces;
  bnd=0;
  sing=0;

  nv=P->n+1;
  fprintf(fid,"Polyhedron %d, R=%lf, %d vertices, %d edges, %d faces, boundary %d, singular %d\n",chart,1.,nv,ne,m,bnd,sing);
  fflush(fid);
  for(i=0;i<P->n;i++)
   {
    fprintf(fid,"Vertex %d (%lf",i,P->v[0+P->dimV*i]);
    for(j=1;j<P->dimV;j++)fprintf(fid,",%lf",P->v[j+P->dimV*i]);
    fprintf(fid,")");
    m=P->nIndices[i];
    index=P->indices[i];
    fprintf(fid,", %d [%d",m,index[0]);
    for(j=1;j<P->nIndices[i];j++)fprintf(fid,",%d",index[j]);
    fprintf(fid,"]\n");
   }
  fprintf(fid,"Vertex %d (%lf",P->n,0.);
  for(j=1;j<P->dimV;j++)fprintf(fid,",0.");
  fprintf(fid,"), %d [ ]\n",0);
  fflush(fid);

  ie=0;
  for(i=0;i<P->n-1;i++)
   {
    for(j=i+1;j<P->n;j++)
     {
      m=MFModifiedPolytopeIntersectIndexSets(P,i,j,indx,e);
      if(m==P->k-1)
       {
        fprintf(fid,"Edge %d (%d,%d)",ie,i,j);
        fprintf(fid,", %d [",m);
        if(m>0)fprintf(fid,"%d",indx[0]);
        for(l=1;l<m;l++)fprintf(fid,",%d",indx[l]);
        fprintf(fid,"]\n");
        ie++;
       }
     }
   }
  free(indx);

  m=P->nFaces;
  for(i=0;i<m;i++)
   {
    j=P->face[i];
    fprintf(fid,"Face %d neighbor %d\n",j,0);
   }

  fflush(fid);

#ifdef MFALLOWVERBOSE
  if(verbose){printf("done %s\n",RoutineName);fflush(stdout);}
#endif

  return;
 }

void MFFreeModifiedPolytope(MFModifiedPolytope P, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFreeModifiedPolytope"};
  int i;


#ifdef MFNOCONFIDENCE
  if(P==NULL)
   {
    sprintf(MFAtlasErrorMsg,"Polytope (argument 1) is NULL");
    MFSetError(e,4,RoutineName,MFAtlasErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif
  if(P->v!=NULL)
   {
    free(P->v);
    P->v=NULL;
   }
  if(P->nIndices!=NULL){free(P->nIndices);P->nIndices=NULL;}
  if(P->mIndices!=NULL){free(P->mIndices);P->mIndices=NULL;}

  if(P->mark!=NULL){free(P->mark);P->mark=NULL;}

  for(i=0;i<P->n;i++)
   {
    if(P->indices[i]!=NULL){free(P->indices[i]);P->indices[i]=NULL;}
   }
  if(P->indices!=NULL){free(P->indices);P->indices=NULL;}
  if(P->face!=NULL){free(P->face);P->face=NULL;}
  if(P->nFaceV!=NULL){free(P->nFaceV);P->nFaceV=NULL;}

  if(P->faceN!=NULL)
   {
    for(i=0;i<P->mFaces;i++) /* was nFaces */
     {
      if(P->faceN[i]!=NULL)MFFreeNVector(P->faceN[i],e);
      P->faceN[i]=NULL;
     }
    free(P->faceN);P->faceN=NULL;
   }
  if(P->faceO!=NULL){free(P->faceO);P->faceO=NULL;}
  free(P);

  return;
 }
