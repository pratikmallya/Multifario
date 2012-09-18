/*
    %W%
    %D% %T%

    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
#include <stdio.h>
#include <MFAtlas.h>
#include <MFDraw.h>
#include <sh.h>
#include <stdlib.h>
#include <math.h>
#include <multifarioConfig.h>
#include <MFErrorHandler.h>
#include <string.h>
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#ifdef __cplusplus
 extern "C" {
#endif

#ifdef HAVE_LIBTIFF
static int tiffavail=1;
#else
static int tiffavail=0;
#endif

int main(int argc, char *argv[])
 {
  FILE *fid;
  int ipoly;
  double *v=(double*)NULL;
  double *vmax=(double*)NULL;
  double *vmin=(double*)NULL;
  double xmin,xmax;
  double ymin,ymax;
  double zmin,zmax;
  int *p=(int*)NULL;
  float *x=(float*)NULL;
  float *y=(float*)NULL;
  float *z=(float*)NULL;
  float x0,y0,z0,x1,y1,z1;
  double tx,ty,tz;
  char name[1024]="";
  int npixels;
  char format[1024]="";
  int c;
  int n,k,nv,ne,nf,np;
  int i,iv,ie,je;
  int iF;
  int mv=0;
  int me=0;
  int v0,v1;
  int full=255;
  int zero=0;
  int w=0;
  MFErrorHandler e;
  int verbose=0;
#ifndef HAVE_GETOPT_H
  int optind=1;
#endif
  double R;

  npixels=1024;
#ifdef HAVE_LIBTIFF
  strcpy(format,"tiff");
#else
  strcpy(format,"ps");
#endif
  w=0;

#ifdef HAVE_GETOPT_H
  while((c=getopt(argc,argv,"n:f:w:h?"))!=EOF)
   {
    switch(c)
     {
      case 'n':
       sscanf(optarg,"%d",&npixels);
       break;
      case 'f':
       sscanf(optarg,"%s",format);
       if(!strcmp(format,"tiff")&&!tiffavail)
        {
         fprintf(stderr,"You requested tiff output, but tiff was not found on your system. Using PostScript output instead\n");
         fflush(stderr);
         strcpy(format,"ps");
        }
       break;
      case 'w':
       sscanf(optarg,"%d",&w);
       break;
      default:
       printf("Usage:\n");
       printf("   DrawPlotfile options Atlas\n");
       printf("      options: \n");
       printf("        -n npixels\n");
       printf("             sets the size of the output image\n");
       printf("        -f format\n");
       printf("             sets the format of the output image (tiff, ps)\n");
       printf("        -w width\n");
       printf("             sets the width (in pixels) of lines in the output image\n");
       printf("        -? \n");
       printf("        -h \n");
       printf("             Prints this message\n");
       return 0;
     }
   }

  if(argc<2)
   {
    fprintf(stderr,"%s, no plotfile specified\n",argv[0]);

    printf("Usage:\n");
    printf("   DrawPlotfile options Atlas\n");
    printf("      options: \n");
    printf("        -n npixels\n");
    printf("             sets the size of the output image\n");
    printf("        -f format\n");
    printf("             sets the format of the output image (tiff, ps)\n");
    printf("        -w width\n");
    printf("             sets the width (in pixels) of lines in the output image\n");
    printf("        -? \n");
    printf("        -h \n");
    printf("             Prints this message\n");
    printf("\n  A \"view\" file may be provided. It must have the same name as the plotfile, with extension \".view\".\n");
    printf("  This file contains a single line with the real values (separated by white space):\n");
    printf("\n a b x0 x1 y0 y1 z0 z1\n\n");
    printf("  a                 is the angle of the viewer around the equator (in degrees)\n");
    printf("  b                 is the angle of the viewer up from the equator toward the north pole (in degrees)\n");
    printf("  x0,x1 y0,y1 z0,z1 define a \"region of interest\". The view will be chosen so that this cube is visible\n");
    return 8;
   }
#else
  if(argc<2)
   {
    fprintf(stderr,"%s, no plotfile specified\n",argv[0]);

    printf("Usage:\n");
    printf("   DrawPlotfile Atlas\n");
    fflush(stdout);
    return 8;
   }
#endif

  printf("Options: npixels %d\n",npixels);
  printf("         plotfile  %s\n",argv[optind]);fflush(stdout);
  printf("         format  \"%s\"\n",format);fflush(stdout);

  strcpy(name,argv[optind]);
  strcat(name,".plotfile");

  fid=fopen(name,"r");
  if(fid==(FILE*)NULL)
   {
    fprintf(stderr,"%s, could not open plot file %s\n",argv[0],name);
    return 0;
   }

  e=MFCreateErrorHandler();

  shSetOutputResolution(npixels,npixels);
  if(format[0]!=0x0)shSetOutputFormat(format);

  shSetOutputFilename(argv[optind]);
  MFDrawInitializeFromFile(argv[optind],e);
  shSetLineWidth(w);

  fscanf(fid,"Dimension of vertices, %d\n",&n);
  if(verbose){printf("Dimension of vertices, %d\n",n);fflush(stdout);}
  if(n<1){printf("%d dimensional vertices can't be plotted\n",n);fflush(stdout);return 12;}
  if(n>3){printf("%d dimensional vertices can't be plotted, using the first three coordinates\n",n);fflush(stdout);}
  fscanf(fid,"Dimension of manifold, %d\n",&k);
  if(verbose){printf("Dimension of manifold, %d\n",k);fflush(stdout);}
  if(k>2){printf("%d dimensional manifolds can't be plotted\n",k);fflush(stdout);return 12;}

  if(k>1)
   {
    shlinc(&full,&full,&full);
    shtric(&zero,&zero,&full);
   }else{
    shlinc(&zero,&zero,&zero);
   }

  vmax=(double*)malloc(n*sizeof(double));
  vmin=(double*)malloc(n*sizeof(double));
  for(i=0;i<n;i++){vmin[i]=1.e308;vmax[i]=-1.e308;}
  xmin=1.e308;xmax=-1.e308;
  ymin=1.e308;ymax=-1.e308;
  zmin=1.e308;zmax=-1.e308;

/* Sample stanza from plotfile 

Polyhedron 0, R=0.100000, 6 vertices, 5 edges, 5 faces, boundary 0, singular 0
Vertex 0 (-0.023740,1.000000,0.054894), 2 [11,17]
Vertex 1 (0.000000,1.000000,-0.070711), 2 [5,25]
Vertex 2 (0.057057,1.000000,-0.013653), 2 [5,7]
Vertex 3 (0.038199,1.000000,0.044094), 2 [7,11]
Vertex 4 (-0.064157,1.000000,-0.006554), 2 [17,25]
Vertex 5 (0.000000,1.000000,0.000000), 0 [ ]
Edge 0 (0,3), 1 [11]
Edge 1 (0,4), 1 [17]
Edge 2 (1,2), 1 [5]
Edge 3 (1,4), 1 [25]
Edge 4 (2,3), 1 [7]
Face 7 neighbor 2
Face 25 neighbor 5
Face 17 neighbor 4
Face 11 neighbor 3
Face 5 neighbor 1

*/
  while(!feof(fid))
   {
    fscanf(fid,"Polyhedron %d, ",&ipoly);
    c=fgetc(fid);
    ungetc(c,fid);
    if(c=='R')
     {
      fscanf(fid,"R=%lf, %d vertices, %d edges, %d faces,  boundary %*d, singular %*d\n",&R,&nv,&ne,&nf);
     }else{
      fscanf(fid,"%d vertices, %d edges, %d faces,  boundary %*d, singular %*d\n",&nv,&ne,&nf);
     }

    if(verbose){printf("Polyhedron %d, %d vertices, %d edges, %d faces\n",ipoly,nv,ne,nf);}
    if(nv>=mv)
     {
      mv=nv;
      v=(double*)realloc((void*)v,n*mv*sizeof(double));
      x=(float*) realloc((void*)x,  mv*sizeof(float));
      y=(float*) realloc((void*)y,  mv*sizeof(float));
      z=(float*) realloc((void*)z,  mv*sizeof(float));
     }
    if(ne>=me)
     {
      me=ne;
      p=(int*)   realloc((void*)p,2*me*sizeof(int));
     }
    for(iv=0;iv<nv;iv++)
     {
      fscanf(fid,"Vertex %*d (%lf",&(v[0+n*iv]));
      for(i=1;i<n;i++)fscanf(fid,",%lf",&(v[i+n*iv]));
      fscanf(fid,"), %*d [%*[ 0-9+-.,]]\n");
      if(verbose)
       {
        printf("Vertex %d (%lf",iv,v[0+n*iv]);
        for(i=1;i<n;i++)printf(",%lf",v[i+n*iv]);
        printf(")\n");fflush(stdout);
       }
      for(i=0;i<n;i++){if(v[i+n*iv]<vmin[i])vmin[i]=v[i+n*iv];if(v[i+n*iv]>vmax[i])vmax[i]=v[i+n*iv];}
      if(0&&n>3)
       {
        tx=v[n-2+n*iv];
        ty=v[n-1+n*iv];
        tz=0.;for(i=0;i<n-2;i++)tz+=v[i+n*iv]*v[i+n*iv];
        v[0+n*iv]=tx;
        v[1+n*iv]=ty;
        v[2+n*iv]=sqrt(tz)/(n-2.);
       }
      if(v[0+n*iv]<xmin)xmin=v[0+n*iv];if(v[0+n*iv]>xmax)xmax=v[0+n*iv];
      if(v[1+n*iv]<ymin)ymin=v[1+n*iv];if(v[1+n*iv]>ymax)ymax=v[1+n*iv];
      if(v[2+n*iv]<zmin)zmin=v[2+n*iv];if(v[2+n*iv]>zmax)zmax=v[2+n*iv];
     }

    for(ie=0;ie<ne;ie++)
     {
      if(k>1)
        fscanf(fid,"Edge %*d (%d,%d), %*d [%*[ 0-9+-.,]]\n",&v0,&v1);
       else if(k==1)
        fscanf(fid,"Edge %*d (%d,%d), %*d []\n",&v0,&v1);
       else
        fscanf(fid,"Edge %*d (%d,%d), %*d [%*[ 0-9+-.,]]\n",&v0,&v1);
      if(verbose){printf("Edge %d (%d,%d)\n",ie,v0,v1);fflush(stdout);}
      x0=v[0+v0*n];
      if(n>1)y0=v[1+v0*n];
       else y0=0.;
      if(n>2)z0=v[2+v0*n];
       else z0=0.;
      x1=v[0+v1*n];
      if(n>1)y1=v[1+v1*n];
       else y1=0.;
      if(n>2)z1=v[2+v1*n];
       else z1=0.;
      shline(&x0,&y0,&z0,&x1,&y1,&z1);
      p[0+2*ie]=v0;
      p[1+2*ie]=v1;
     }

    if(k==2 && ne==nv-1)
     {
      np=0;
      iv=p[0];
      x[np]=v[0+iv*n];
      if(n>1)y[np]=v[1+iv*n];
       else y[np]=0.;
      if(n>2)z[np]=v[2+iv*n];
       else z[np]=0.;
      if(verbose){printf("   Point %d = (%f,%f,%f)\n",np,x[np],y[np],z[np]);fflush(stdout);}
      np++;

      iv=p[1];
      x[np]=v[0+iv*n];
      if(n>1)y[np]=v[1+iv*n];
       else y[np]=0.;
      if(n>2)z[np]=v[2+iv*n];
       else z[np]=0.;
      if(verbose){printf("   Point %d = (%f,%f,%f)\n",np,x[np],y[np],z[np]);fflush(stdout);}
      np++;

      je=0;

      ie=1;
      while(np<ne)
       {
        if(verbose){printf("     ie=%d, ne=%d, np=%d\n",ie,ne,np);fflush(stdout);}
        if(p[0+2*ie]==iv&&je!=ie)
         {
          iv=p[1+2*ie];
          x[np]=v[0+iv*n];
          if(n>1)y[np]=v[1+iv*n];
           else y[np]=0.;
          if(n>2)z[np]=v[2+iv*n];
           else z[np]=0.;
          if(verbose){printf("   Point %d = (%f,%f,%f)\n",np,x[np],y[np],z[np]);fflush(stdout);}
          np++;
          je=ie;
         }else if(p[1+2*ie]==iv&&je!=ie)
         {
          iv=p[0+2*ie];
          x[np]=v[0+iv*n];
          if(n>1)y[np]=v[1+iv*n];
           else y[np]=0.;
          if(n>2)z[np]=v[2+iv*n];
           else z[np]=0.;
          if(verbose){printf("   Point %d = (%f,%f,%f)\n",np,x[np],y[np],z[np]);fflush(stdout);}
          np++;
          je=ie;
         }
        ie++;
        if(ie>ne-1)ie=0;
       }
     }
    if(verbose&&k==2){printf("call shpg %d\n",np);fflush(stdout);}

    if(k==2)shpg(&np,x,y,z);

    for(iF=0;iF<nf;iF++)
     {
      fscanf(fid,"Face %*d neighbor %*d\n");
      if(verbose){printf("   Face %d\n",iF);fflush(stdout);}
     }
   }

  if(n>3)
   {
    for(i=0;i<n;i++){printf("Direction %d, interval [%lf,%lf]\n",i,vmin[i],vmax[i]);fflush(stdout);
   }
  printf("Plot variables in [(%lf,%lf,%lf),(%lf,%lf,%lf)]\n",xmin,ymin,zmin,xmax,ymax,zmax);fflush(stdout);}

  free(vmin);
  free(vmax);

  MFDrawDisplay(e);
  MFDrawClose(e);
  fclose(fid);

  MFFreeErrorHandler(e);

  return 0;
 }

#ifdef __cplusplus
}
#endif
