/* 
    @(#)DrawAtlasTS.c	1.7
    03/02/26 16:38:10
   
    PROGRAM NAME:  multifario

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#include <multifarioConfig.h>
#include <MFAtlas.h>
#include <MFImplicitMF.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <MFPrint.h>
#include <MFDraw.h>
#include <MFEnumPolytope.h>
#include <MFEnumDualPolytope.h>
#include <MFDX.h>
#include <MFErrorHandler.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <sh.h>
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <unistd.h>

#ifdef __cplusplus
 extern "C" {
#endif

extern int optind;
extern char *optarg;

void line(float,float,float,float,float,float);
void arrow(float,float,float,float,float,float,int);
void string(float,float,float,char*);

#define YY .3

int main(int argc, char *argv[])
 {
  MFAtlas S;
  int n;
  MFEnumDualPolytope P;
  FILE *fid;
  char name[1024]="";
  int pendula,periodic,doatlas,npixels;
  char format[10]="";
  int zero=0;
  int one=1;
  int two=2;
  int side=1;
  int gray;
  int full=255;
  float x,y,z,dx,dy,dz;
  int c;
  MFErrorHandler e;

  if(argc<2){printf("Usage %s filename\n\n",argv[0]);fflush(stdout);
             printf("Options: -n number of pixels in output image\n");
             printf("Options: -f format (\"tiff\" if libtiff is installed or \"ps\")\n");fflush(stdout);return 8;}

  periodic=0;
  doatlas=1;
  npixels=1024;
  format[0]=0x0;
  while((c=getopt(argc,argv,"mn:f:p"))!=EOF)
   {
    switch(c)
     {
      case 'p':
       periodic=1;
       break;
      case 'n':
       sscanf(optarg,"%d",&npixels);
       break;
      case 'm':
       doatlas=0;
       break;
      case 'f':
       sscanf(optarg,"%s",format);
       break;
     }
   }

  printf("%s %s\n",argv[0],argv[optind]);

  strcpy(name,argv[optind]);
  pendula= (strstr(name,"Pendula")!=(char*)NULL);
  strcpy(name,argv[optind]);
  strcat(name,".atlas");

  printf("Options: Periodic %d\n",periodic);
  printf("         npixels  %d\n",npixels);
  printf("         doatlas  %d\n",doatlas);
  printf("         atlas  %s\n",argv[optind]);fflush(stdout);
  printf("         format  %s\n",format);fflush(stdout);

  fid=fopen(name,"r");
  if(fid==(FILE*)NULL)
   {
    printf("Error, could not open file %s, %s\n",name,strerror(errno));
    fflush(stdout);
    return 12;
   }
  printf("Reading Atlas %s\n",argv[optind]);fflush(stdout);
  S=(MFAtlas)NULL;
  if(doatlas)S=MFReadAtlas(fid,e);
  fclose(fid);
  printf("Done reading Atlas, %d charts\n",MFAtlasNumberOfCharts(S,e));fflush(stdout);

  shSetOutputResolution(npixels,npixels);
  if(format[0]!=0x0)shSetOutputFormat(format);
  strcpy(name,argv[optind]);
  strcat(name,"TS");
  shSetOutputFilename(name);

  e=MFCreateErrorHandler();

  MFDrawInitializeFromFile(argv[optind],e);
  {
   float a,d,s;
   int n;
   a=.8;d=.4;s=0.4;n=5;
   shsrfp(&a,&d,&s,&n);
  }
  if(pendula)
   {
    if(periodic)
     {
      MFPendulaPeriodic(1,e);
      side=-1;
      x=-.0001; y=0.; z=0.;
      dx=-1.;dy=0.;dz=0.;
      shpln(&one,&x,&y,&z,&dx,&dy,&dz,&side);
 
      side=-1;
      x=.1*3.1415926+.0001; y=0.; z=0.;
      dx= 1.;dy=0.;dz=0.;
      shpln(&two,&x,&y,&z,&dx,&dy,&dz,&side);
     }

    shlinc(&zero,&zero,&zero);

    arrow(0., 0.,0., 3.,  0.,0.,2);
    arrow(0., 0.,0., 0., 12.,0.,2);

    arrow(0.,0.,0., 2.5,2.5,0.,2);
    arrow(0.,0.,0.,-2.5,2.5,0.,2);

    arrow(0.,0.,0.,0.,0.,26.,1);

    line(-2., 0.,0.,2., 0.,0.);

    if(periodic)shnpln(&two);

    line(-1.,10.,0., 2.,10.,0.);
    line(-1.,10.,0.,-1., 0.,0.);

    line( 0., 0.,0., 0.,10.,0.);
    line( 1., 0.,0., 1.,10.,0.);
    line( 2., 0.,0., 2.,10.,0.);

    line(-1., 0.,25.,2., 0.,25.);

    line(-1., 0.,0.,-1., 0.,25.);
    line( 1., 0.,0., 1., 0.,25.);
    line( 2., 0.,0., 2., 0.,25.);

/*  string(0.,0.,6.1,"/Helvetica findfont 15.000000 scalefont setfont");
    string(3.3,0.,0.,"Ic");
    string(0.,12.2,0.,"I");
    string( 2.9,2.9,0.,"I0");
    string(-2.6,2.6,0.,"I1");
    string(0.,0.,16.1,"T");
    string(0.,0.,6.1,"/Symbol findfont 15.000000 scalefont setfont");
    if(!periodic)string( 2.,10.6,0.,"2kp");
    string( 1.,10.6,0.,"kp");
    if(!periodic)string(-1.,10.6,0.,"-kp");*/
   }

  if(doatlas&& S!=(MFAtlas)NULL)MFDrawAtlasTS(S,e);
  MFDrawDisplay(e);

  if(S!=(MFAtlas)NULL)MFFreeAtlas(S,e);
  MFDrawClose(e);

  MFFreeErrorHandler(e);

  return(0);
 }

void arrow(float Ic0,float In0,float w0,float Ic1,float In1,float w1, int dir)
 {
  double Pi=3.1415926;
  double kappa=.1;
  float d;
  float dx,dy,dz;
  float nx,ny,nz;

  dx=2*Pi*kappa*(Ic1-Ic0);
  dy=.1*(In1-In0);
  dy=YY*(In1-In0);
/*dz=.4*(w1-w0);*/
  dz=.2*(w1-w0);
  d=.01/sqrt(dx*dx+dy*dy+dz*dz);
  dx=d*dx;
  dy=d*dy;
  dz=d*dz;

  line(Ic0,In0,w0,Ic1,In1,w1);
  switch(dir)
   {
    case 0:
     nx= 0.;
     ny= dz;
     nz=-dy;
     break;
    case 1:
     nx= dz;
     ny= 0.;
     nz=-dy;
     break;
    case 2:
     nx= dy;
     ny=-dx;
     nz= 0.;
     break;
   }
  dx=dx/Pi/kappa;
/*dy=dy/.1*/;
  dy=dy/YY;
/*dz=dz/.4;*/
  dz=dz/.2;

  nx=nx/Pi/kappa;
/*ny=ny/.1;*/
  ny=ny/YY;
/*nz=nz/.4;*/
  nz=nz/.2;

  dx=dx*6.;
  dy=dy*6.;
  dz=dz*6.;

  line(Ic1,In1,w1,Ic1-dx+nx,In1-dy+ny,w1-dz+nz);
  line(Ic1-dx+nx,In1-dy+ny,w1-dz+nz,Ic1-dx-nx,In1-dy-ny,w1-dz-nz);
  line(Ic1-dx-nx,In1-dy-ny,w1-dz-nz,Ic1,In1,w1);
  return;
 }

void line(float Ic0,float In0,float w0,float Ic1,float In1,float w1)
 {
  float x0,y0,z0;
  float x1,y1,z1;
  double Pi=3.1415926;
  double kappa=.1;

  x0=Pi*kappa*Ic0;
  x0=2*Pi*kappa*Ic0;
  y0=.1*In0;
  y0=YY*In0;
  z0=.4*Pi*w0;
  z0=.2*w0;
  x1=2*Pi*kappa*Ic1;
  y1=.1*In1;
  y1=YY*In1;
  z1=.4*Pi*w1;
  z1=.2*w1;

  shline(&x0,&y0,&z0,&x1,&y1,&z1);
  return;
 }

void string(float Ic,float In,float w,char *str)
 {
  float x,y,z;
  double Pi=3.1415926;
  double kappa=.1;

  x=Pi*kappa*Ic;
  x=2*Pi*kappa*Ic;
  y=.1*In;
  y=YY*In;
  z=.2*w;
  shstr(x,y,z,str);

  return;
 }

#ifdef __cplusplus
}
#endif
