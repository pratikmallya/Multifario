#include <stdio.h>
#include <MFAtlas.h>
#include <MFDraw.h>
#include <sh.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef __cplusplus
 extern "C" {
#endif

void MFPlotfileToVBM(FILE*,char*);
void PFPolytopes(FILE*,int,int,int*,int*,int*,double**,int**,int*,int*,int**,int**,int*,int*,int**,int***,int**);
void MFMeshToVBMFile(char*,int,int,int,double*,int*,int,int*,int*,int,int*,int**,int*);

int main(int argc, char *argv[])
 {
  FILE *fid;
  char file[4096];

  if(argc<2){printf("Usage %s filename\n",argv[0]);fflush(stdout);return 8;}

  strcpy(file,argv[1]);
  strcat(file,".plotfile");
  fid=fopen(file,"r");
  if(fid==(FILE*)NULL){printf("Error opening file %s\n",file);fflush(stdout);return 8;}

  MFPlotfileToVBM(fid,argv[1]);

  fclose(fid);

  return 0;
 }

void MFPlotfileToVBM(FILE *fid,char *name)
 {
  static char RoutineName[]={"PlotfileToVBM"};
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

  printf("%s\n",RoutineName);fflush(stdout);

  iV=0;
  iV0=0;
  mV=0;

  iE=0;
  mE=0;

  iF=0;
  mF=0;

  fscanf(fid,"Dimension of vertices, %d\n",&m);
  fscanf(fid,"Dimension of manifold, %d\n",&k);

  printf("Add in Polytopes\n");fflush(stdout);
  PFPolytopes(fid,m,k,&iV,&iV0,&mV,&v,&vt,&iE,&mE,&ev,&et,&iF,&mF,&nfv,&fv,&ft);
  printf("  There are now %d pts, %d edges %d faces\n",iV,iE,iF);fflush(stdout);
  printf("Dump it to %s\n",name);fflush(stdout);
  MFMeshToVBMFile(name,k,iV,m,v,vt,iE,ev,et,iF,nfv,fv,ft);
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

void PFPolytopes(FILE *fid,int m, int k, int *piV,int *piV0,int *pmV,double **pv,int **pvt,int *piE,int *pmE,int **pev,int **pet, int *piF,int *pmF,int **pnfv,int ***pfv, int **pft)
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
  int c;

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
      fscanf(fid,"R=%*lf, %d vertices, %d edges, %d faces,  boundary %d, singular %d\n",&nv,&ne,&nf,&bnd,&sing);
     }else{
      fscanf(fid,"%d vertices, %d edges, %d faces,  boundary %d, singular %d\n",&nv,&ne,&nf,&bnd,&sing);
     }

    if(verbose){printf("Polyhedron %d, %d vertices, %d edges, %d faces\n",ipoly,nv,ne,nf);}
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
       else
        fscanf(fid,"Edge %*d (%d,%d), %*d []\n",&v0,&v1);
      ev[0+2*iE]=iV0+v0;
      ev[1+2*iE]=iV0+v1;
      et[iE]=0;
      if(verbose){printf("  Edge %d (%d,%d)\n",iE,ev[0+2*iE],ev[1+2*iE]);fflush(stdout);}
      iE++;
     }

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

    for(i=0;i<nf;i++)fscanf(fid,"Face %*d neighbor %*d\n");
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

#ifdef __cplusplus
}
#endif
