#include <MFAtlas.h>

struct IsoFlowData;

double IMFIsoSurfaceFlowEvaluateG  (struct IsoFlowData *f,double *u);
void   IMFIsoSurfaceFlowEvaluateDG (struct IsoFlowData *f,double *u,double *gradg);
void   IMFIsoSurfaceFlowEvaluateDDG(struct IsoFlowData *f,double *u,double **dgradg);
void   IMFIsoSurfaceFlowEvaluateW  (struct IsoFlowData *f,double *u,double *w);
void   IMFIsoSurfaceFlowEvaluateDW (struct IsoFlowData *f,double *u,double **dw);

MFAtlas MFCoverFlowOnManifold(MFImplicitMF M, IMFFlow F,char *name, int nInitial, MFNVector *u0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, double dmin, int maxInterp, int maxCharts, MFErrorHandler e);

IMFFlow IMFCreateIosurfaceFlow(double (*g)(double*),void (*dg)(double*,double*),void (*ddg)(double*,double**),void (*w)(double*,double*),void (*dw)(double*,double**),MFErrorHandler e);

double gTorus(double *x);
void dgTorus(double *x, double *dg);
void ddgTorus(double *x, double **dg);

void w(double *x, double *w);
void dw(double *x, double **dw);

static void FF(int *n,double *z,int *m,double *f, void *data, MFErrorHandler e)
 {
  f[0]=gTorus(z);
  return;
 }

static void dF(int *n,double *z,int *m,double *f, void *data, MFErrorHandler e)
 {
  dgTorus(z,f);

  return;
 }

static void ddF(int *n,double *z,int *m,double *f, void *data, MFErrorHandler e)
 {
  static double *g[3]={NULL,NULL,NULL};

  g[0]=(double*)realloc(g[0],3*sizeof(double));
  g[1]=(double*)realloc(g[1],3*sizeof(double));
  g[2]=(double*)realloc(g[2],3*sizeof(double));

  ddgTorus(z,g);

  f[0]=g[0][0];
  f[1]=g[0][1];
  f[2]=g[0][2];
  f[3]=g[1][0];
  f[4]=g[1][1];
  f[5]=g[1][2];
  f[6]=g[2][0];
  f[7]=g[2][1];
  f[8]=g[2][2];

  return;
 }

static void findFixedPoint(double *x0,double *x, MFErrorHandler e);
void MFSVD(int n, double *A, double *vL, double *s, double *vR, MFErrorHandler e);
void MFEV (int n, double *A, double *vL, double *sr, double *si, double *vR, MFErrorHandler e);

/* ----------------------------------------------------------------------------------------------- */

#define Rxy .8
#define Rz .5

int main(int argc, char *argv[])
 {
  MFAtlas A;
  MFImplicitMF M;
  IMFFlow F;
  char *name;
  int nInitial;
  MFNVector u0[10000];
  MFKVector p0;
  MFNVector u;
  MFNRegion Omega;
  double eps;
  double dt;
  double tmax;
  int maxInterp;
  int maxCharts;
  MFErrorHandler e;
  int i,j,k,n;
  double x,y,z;
  double theta,phi;
  double X[3];
  int N;
  int fxpt;
  double x0[3],xf[3];
  double v[3];
  double gradg[3],gradgrad;
  double J[9],Q[9],dFlow[9],Flow[3];
  double vL[9],vR[9],sr[3],si[3];
  MFNKMatrix Phi;
  double phi0[3];
  double phi1[3];
  MFNVector tV;

  double a[4];
  double *dv[3];
  dv[0]=(double*)malloc(3*sizeof(double));
  dv[1]=(double*)malloc(3*sizeof(double));
  dv[2]=(double*)malloc(3*sizeof(double));

  e=MFCreateErrorHandler();

  F=IMFCreateIosurfaceFlow(gTorus,dgTorus,ddgTorus,w,dw,e);

  M=MFIMFCreateAlgebraicSubroutine(3,2,FF,dF,ddF,NULL,e);
  MFIMFSetR(M,.03,e);

  n=0;

#if 1
  for(fxpt=0;fxpt<4;fxpt++)
   {
    switch(fxpt)
     {
      case 0:
       theta=140;
       phi  =122;
       break;
      case 1:
       theta=300;
       phi  =230;
       break;
      case 2:
       theta=32;
       phi  =65;
       break;
      case 3:
       theta=230;
       phi  =300;
       break;
     }
    printf("(theta,phi)=(%lf,%lf)\n",theta,phi);fflush(stdout);
    printf("\n");fflush(stdout);

    u0[n]=MFIMFVectorFactory(M,e);
    x0[0]=cos(theta)*(Rxy+Rz*cos(phi));
    x0[1]=sin(theta)*(Rxy+Rz*cos(phi));
    x0[2]=                Rz*sin(phi);
    findFixedPoint(x0,xf,e);

    MFNVSetC(u0[n],0,xf[0],e);
    MFNVSetC(u0[n],1,xf[1],e);
    MFNVSetC(u0[n],2,xf[2],e);

    Phi=MFIMFTangentSpace(M,u0[n],e);

    tV=MFMColumn(Phi,0,e);
    phi0[0]=MFNV_C(tV,0,e);
    phi0[1]=MFNV_C(tV,1,e);
    phi0[2]=MFNV_C(tV,2,e);
    MFFreeNVector(tV,e);

    tV=MFMColumn(Phi,1,e);
    phi1[0]=MFNV_C(tV,0,e);
    phi1[1]=MFNV_C(tV,1,e);
    phi1[2]=MFNV_C(tV,2,e);
    MFFreeNVector(tV,e);

    dw(xf,dv);

    a[0+2*0]=phi0[0]*(dv[0][0]*phi0[0]+dv[0][1]*phi0[1]+dv[0][2]*phi0[2])
            +phi0[1]*(dv[1][0]*phi0[0]+dv[1][1]*phi0[1]+dv[1][2]*phi0[2])
            +phi0[2]*(dv[2][0]*phi0[0]+dv[2][1]*phi0[1]+dv[2][2]*phi0[2]);
    a[0+2*1]=phi0[0]*(dv[0][0]*phi1[0]+dv[0][1]*phi1[1]+dv[0][2]*phi1[2])
            +phi0[1]*(dv[1][0]*phi1[0]+dv[1][1]*phi1[1]+dv[1][2]*phi1[2])
            +phi0[2]*(dv[2][0]*phi1[0]+dv[2][1]*phi1[1]+dv[2][2]*phi1[2]);
    a[1+2*0]=phi1[0]*(dv[0][0]*phi0[0]+dv[0][1]*phi0[1]+dv[0][2]*phi0[2])
            +phi1[1]*(dv[1][0]*phi0[0]+dv[1][1]*phi0[1]+dv[1][2]*phi0[2])
            +phi1[2]*(dv[2][0]*phi0[0]+dv[2][1]*phi0[1]+dv[2][2]*phi0[2]);
    a[1+2*1]=phi1[0]*(dv[0][0]*phi1[0]+dv[0][1]*phi1[1]+dv[0][2]*phi1[2])
            +phi1[1]*(dv[1][0]*phi1[0]+dv[1][1]*phi1[1]+dv[1][2]*phi1[2])
            +phi1[2]*(dv[2][0]*phi1[0]+dv[2][1]*phi1[1]+dv[2][2]*phi1[2]);
 
    dgTorus(xf,gradg);
    gradgrad=gradg[0]*gradg[0]+gradg[1]*gradg[1]+gradg[2]*gradg[2];
    IMFEvaluateFlow(F,u0[n],(MFKVector)NULL, Flow, e);
    IMFEvaluateDerivativeOfFlow(F,u0[n],(MFKVector)NULL, J, e);

    for(i=0;i<3;i++)
     {
      for(j=0;j<3;j++)
       {
        Q[i+3*j]=0;
        if(i==j)Q[i+3*j]=1.;
        Q[i+3*j]=Q[i+3*j]-gradg[i]*gradg[j]/gradgrad;
       }
     }

    printf("Q.gradg 0, %lf\n",Q[0+3*0]*gradg[0]+Q[1+3*0]*gradg[1]+Q[2+3*0]*gradg[2]);fflush(stdout);
    printf("        1, %lf\n",Q[0+3*1]*gradg[0]+Q[1+3*1]*gradg[1]+Q[2+3*1]*gradg[2]);fflush(stdout);
    printf("        2, %lf\n",Q[0+3*2]*gradg[0]+Q[1+3*2]*gradg[1]+Q[2+3*2]*gradg[2]);fflush(stdout);
    printf("\n");fflush(stdout);

    for(i=0;i<3;i++)
     {
      for(j=0;j<3;j++)
       {
        dFlow[i+3*j]=0;
        for(k=0;k<3;k++)
         dFlow[i+3*j]+=Q[i+3*k]*J[k+3*j];
       }
      J[i+3*j]=dFlow[i+3*j];
     }

    MFEV(3, dFlow,vL,sr,si,vR,e);
    printf("ev 0, %lf+i*%lf, (%lf,%lf,%lf)\n",sr[0],si[0],vR[0+3*0],vR[1+3*0],vR[2+3*0]);fflush(stdout);
    printf("   1, %lf+i*%lf, (%lf,%lf,%lf)\n",sr[1],si[1],vR[0+3*1],vR[1+3*0],vR[2+3*1]);fflush(stdout);
    printf("   2, %lf+i*%lf, (%lf,%lf,%lf)\n",sr[2],si[2],vR[0+3*2],vR[1+3*0],vR[2+3*2]);fflush(stdout);
    printf("dg (%lf,%lf,%lf)\n",gradg[0]/sqrt(gradgrad),gradg[1]/sqrt(gradgrad),gradg[2]/sqrt(gradgrad));fflush(stdout);
    printf("\n");fflush(stdout);
  
    printf("ev.gradg 0, %lf\n",vR[0+3*0]*gradg[0]+vR[1+3*0]*gradg[1]+vR[2+3*0]*gradg[2]);fflush(stdout);
    printf("         1, %lf\n",vR[0+3*1]*gradg[0]+vR[1+3*1]*gradg[1]+vR[2+3*1]*gradg[2]);fflush(stdout);
    printf("         2, %lf\n",vR[0+3*2]*gradg[0]+vR[1+3*2]*gradg[1]+vR[2+3*2]*gradg[2]);fflush(stdout);
    printf("\n");fflush(stdout);

    for(i=0;i<9;i++)dFlow[i]=J[i];

    MFSVD(3, dFlow,vL,sr,vR,e);
    printf("svd 0, %lf, (%lf,%lf,%lf)\n",sr[0],vR[0+3*0],vR[1+3*0],vR[2+3*0]);fflush(stdout);
    printf("    1, %lf, (%lf,%lf,%lf)\n",sr[1],vR[0+3*1],vR[1+3*0],vR[2+3*1]);fflush(stdout);
    printf("    2, %lf, (%lf,%lf,%lf)\n",sr[2],vR[0+3*2],vR[1+3*0],vR[2+3*2]);fflush(stdout);
    printf("\n");fflush(stdout);
  
    printf("sv.gradg 0, %lf\n",vR[0+3*0]*gradg[0]+vR[1+3*0]*gradg[1]+vR[2+3*0]*gradg[2]);fflush(stdout);
    printf("         1, %lf\n",vR[0+3*1]*gradg[0]+vR[1+3*1]*gradg[1]+vR[2+3*1]*gradg[2]);fflush(stdout);
    printf("         2, %lf\n",vR[0+3*2]*gradg[0]+vR[1+3*2]*gradg[1]+vR[2+3*2]*gradg[2]);fflush(stdout);
    printf("\n");fflush(stdout);
  
    for(i=0;i<3;i++)
     {
      gradgrad=0.;
      for(j=0;j<3;j++)
       {
        for(k=0;k<3;k++)
         {
          gradgrad+=gradg[j]*J[j+3*k]*vR[k+3*i];
         }
       }
      printf("%d gradG^T dFlow vR %lf\n",i,gradgrad);fflush(stdout);
     }
    printf("\n");

    MFEV(2, a,vL,sr,si,vR,e);
    printf("  phi0 (%lf,%lf,%lf)\n",phi0[0],phi0[1],phi0[2]);
    printf("  phi1 (%lf,%lf,%lf)\n",phi1[0],phi1[1],phi1[2]);
    printf("  nrm  (%lf,%lf,%lf)\n",gradg[0]/sqrt(gradgrad),gradg[1]/sqrt(gradgrad),gradg[2]/sqrt(gradgrad));
    printf("\n");fflush(stdout);
    printf("in TS\n");
    printf("ev 0, %lf+i*%lf, (%lf,%lf) (%lf,%lf,%lf)\n",
                         sr[0],si[0],
                         vR[0+2*0],vR[1+2*0],
                         phi0[0]*vR[0+2*0]+phi1[0]*vR[1+2*0],
                         phi0[1]*vR[0+2*0]+phi1[1]*vR[1+2*0],
                         phi0[2]*vR[0+2*0]+phi1[2]*vR[1+2*0]);fflush(stdout);
    printf("   1, %lf+i*%lf, (%lf,%lf) (%lf,%lf,%lf)\n",
                         sr[1],si[1],
                         vR[0+2*1],vR[1+2*1],
                         phi0[0]*vR[0+2*1]+phi1[0]*vR[1+2*1],
                         phi0[1]*vR[0+2*1]+phi1[1]*vR[1+2*1],
                         phi0[2]*vR[0+2*1]+phi1[2]*vR[1+2*1]);fflush(stdout);
  
    printf("Guess %d, (%lf,%lf,%lf) g=%le, w.gradg=%le\n",fxpt,xf[0],xf[1],xf[2],
                                                          gTorus(xf),
                                                          v[0]*gradg[0]+v[1]*gradg[1]+v[2]*gradg[2]);fflush(stdout);
    printf("\n");fflush(stdout);
 
    printf("----------------------------------------------\n");fflush(stdout);
    printf("\n");fflush(stdout);
    switch(fxpt)
     {
      case 0:
      case 1:

/*  ev 0 > 0, ev 1 < 0 */
/*
       u=MFIMFVectorFactory(M,e);
       for(i=1;i<2;i++)
        {
         u0[n]=MFIMFVectorFactory(M,e);

         x=xf[0]+( .03*i*(phi0[0]*vR[0+2*0]+phi1[0]*vR[1+2*0]) + 0.00*(phi0[0]*vR[0+2*1]+phi1[0]*vR[1+2*1]));
         y=xf[1]+( .03*i*(phi0[1]*vR[0+2*0]+phi1[0]*vR[1+2*0]) + 0.00*(phi0[0]*vR[0+2*1]+phi1[0]*vR[1+2*1]));
         z=xf[2]+( .03*i*(phi0[2]*vR[0+2*0]+phi1[0]*vR[1+2*0]) + 0.00*(phi0[0]*vR[0+2*1]+phi1[0]*vR[1+2*1]));

         MFNVSetC(u,0,x,e);
         MFNVSetC(u,1,y,e);
         MFNVSetC(u,2,z,e);
         MFIMFProject(M,u,Phi,u0[n],e);

         printf("Initial Point %d, ",n);MFPrintNVector(stdout,u0[n],e);printf("\n");fflush(stdout);
         n++;

         u0[n]=MFIMFVectorFactory(M,e);

         x=xf[0]+(-.03*i*(phi0[0]*vR[0+2*0]+phi1[0]*vR[1+2*0]) + 0.00*(phi0[0]*vR[0+2*1]+phi1[0]*vR[1+2*1]));
         y=xf[1]+(-.03*i*(phi0[1]*vR[0+2*0]+phi1[0]*vR[1+2*0]) + 0.00*(phi0[0]*vR[0+2*1]+phi1[0]*vR[1+2*1]));
         z=xf[2]+(-.03*i*(phi0[2]*vR[0+2*0]+phi1[0]*vR[1+2*0]) + 0.00*(phi0[0]*vR[0+2*1]+phi1[0]*vR[1+2*1]));

         MFNVSetC(u,0,x,e);
         MFNVSetC(u,1,y,e);
         MFNVSetC(u,2,z,e);
         MFIMFProject(M,u,Phi,u0[n],e);

         printf("Initial Point %d, ",n);MFPrintNVector(stdout,u0[n],e);printf("\n");fflush(stdout);
         n++;
        }
       MFFreeNVector(u,e);
 */
       break;

      case 2:
      case 3:

/*  spiral out */

       u=MFIMFVectorFactory(M,e);
       N=150;
       for(i=0;i<N;i++)
        {
         theta=2*3.14159*(0.+i/(N-1.));
         x=xf[0]+.3*( cos(theta)*phi0[0] + sin(theta)*phi1[0]);
         y=xf[1]+.3*( cos(theta)*phi0[1] + sin(theta)*phi1[1]);
         z=xf[2]+.3*( cos(theta)*phi0[2] + sin(theta)*phi1[2]);
         
         MFNVSetC(u,0,x,e);
         MFNVSetC(u,1,y,e);
         MFNVSetC(u,2,z,e);

         u0[n]=MFIMFVectorFactory(M,e);
         MFIMFProject(M,u,Phi,u0[n],e);

         printf("Initial Point %d, ",n);MFPrintNVector(stdout,u0[n],e);printf("\n");fflush(stdout);
         n++;
        }

       N=0; /* Was 20 */
       for(i=-N;i<N;i++)
        {
         if(fxpt==2)
           theta=2*3.14159*(0.5+.025*i/(N-1.));
          else
           theta=2*3.14159*(0.0+.025*i/(N-1.));
          
         x=xf[0]+.5*( cos(theta)*phi0[0] + sin(theta)*phi1[0]);
         y=xf[1]+.5*( cos(theta)*phi0[1] + sin(theta)*phi1[1]);
         z=xf[2]+.5*( cos(theta)*phi0[2] + sin(theta)*phi1[2]);
         
         MFNVSetC(u,0,x,e);
         MFNVSetC(u,1,y,e);
         MFNVSetC(u,2,z,e);

         u0[n]=MFIMFVectorFactory(M,e);
         MFIMFProject(M,u,Phi,u0[n],e);

         printf("Initial Point %d, ",n);MFPrintNVector(stdout,u0[n],e);printf("\n");fflush(stdout);
         n++;
        }
 
       MFFreeNVector(u,e);
       break; 
     }
   }

/*
  N=20;
  for(i=-N/2;i<N;i++)
   {
    u0[n]=MFIMFVectorFactory(M,e);
    theta=2.*3.151926*(310./360.+.1*i/(N-1.));
    phi  =2.*3.151926*(270./360.+0.*i/(N-1.));
    x=cos(theta)*(Rxy+Rz*cos(phi));
    y=sin(theta)*(Rxy+Rz*cos(phi));
    z=Rz*sin(phi);
    MFNVSetC(u0[n],0,x,e);
    MFNVSetC(u0[n],1,y,e);
    MFNVSetC(u0[n],2,z,e);

    X[0]=x; X[1]=y; X[2]=z;
    printf("Line Guess %d, g=%le, theta=%lf, phi=%lf\n",n,gTorus(X),theta,phi);fflush(stdout);
    n++;
   }
 */
#endif

#if 0
  N=7;
  MFIMFSetR(M,.03,e);
  for(i=0;i<N;i++)
   {
    for(j=0;j<20*N;j++)
     {
      u0[n]=MFIMFVectorFactory(M,e);
      theta=2.*3.151926*j/(20*N-1.);
      phi  =2.*3.151926*i/(N-1.);
      x=cos(theta)*(Rxy+Rz*cos(phi));
      y=sin(theta)*(Rxy+Rz*cos(phi));
      z=Rz*sin(phi);
      MFNVSetC(u0[n],0,x,e);
      MFNVSetC(u0[n],1,y,e);
      MFNVSetC(u0[n],2,z,e);
  
      X[0]=x; X[1]=y; X[2]=z;
      printf("Guess %d, g=%le\n",n,gTorus(X));fflush(stdout);
      n++;
    }
  }

  N=80;
  for(i=0;i<N;i++)
   {
    u0[n]=MFIMFVectorFactory(M,e);
    theta=2.*3.151926*0.25;
    phi  =2.*3.151926*(1.*i/(N-1.));
    x=cos(theta)*(Rxy+Rz*cos(phi));
    y=sin(theta)*(Rxy+Rz*cos(phi));
    z=Rz*sin(phi);
    MFNVSetC(u0[n],0,x,e);
    MFNVSetC(u0[n],1,y,e);
    MFNVSetC(u0[n],2,z,e);

    X[0]=x; X[1]=y; X[2]=z;
    printf("Guess %d, g=%le\n",n,gTorus(X));fflush(stdout);
    n++;
  }

  N=90;
  for(i=0;i<N;i++)
   {
    u0[n]=MFIMFVectorFactory(M,e);
    theta=2.*3.151926*0.75;
    phi  =2.*3.151926*(9.*i/(N-1.));
    x=cos(theta)*(Rxy+Rz*cos(phi));
    y=sin(theta)*(Rxy+Rz*cos(phi));
    z=Rz*sin(phi);
    MFNVSetC(u0[n],0,x,e);
    MFNVSetC(u0[n],1,y,e);
    MFNVSetC(u0[n],2,z,e);

    X[0]=x; X[1]=y; X[2]=z;
    printf("Guess %d, g=%le\n",n,gTorus(X));fflush(stdout);
    n++;
   }

  N=80;
  for(i=0;i<N;i++)
   {
    u0[n]=MFIMFVectorFactory(M,e);
    theta= 2.*3.151926*(.2+.8*i/(N-1.));
    phi  = 2.*3.151926*(.1-.5);
    x=cos(theta)*(Rxy+Rz*cos(phi));
    y=sin(theta)*(Rxy+Rz*cos(phi));
    z=Rz*sin(phi);
    MFNVSetC(u0[n],0,x,e);
    MFNVSetC(u0[n],1,y,e);
    MFNVSetC(u0[n],2,z,e);

    X[0]=x; X[1]=y; X[2]=z;
    printf("Guess %d, g=%le\n",n,gTorus(X));fflush(stdout);
    n++;
   }

  for(i=0;i<N;i++)
   {
    u0[n]=MFIMFVectorFactory(M,e);
    theta= 2.*3.151926*(.8+.1*i/(N-1.));
    phi  = 2.*3.151926*(.1+.25);
    x=cos(theta)*(Rxy+Rz*cos(phi));
    y=sin(theta)*(Rxy+Rz*cos(phi));
    z=Rz*sin(phi);
    MFNVSetC(u0[n],0,x,e);
    MFNVSetC(u0[n],1,y,e);
    MFNVSetC(u0[n],2,z,e);

    X[0]=x; X[1]=y; X[2]=z;
    printf("Guess %d, g=%le\n",n,gTorus(X));fflush(stdout);
    n++;
   }
#endif

  p0=(MFKVector)NULL;

  Omega=MFNRegionCreateHyperCube(3,1.1,e);

  eps=.4;
  dt=.01;
  tmax= 2000.;
  maxInterp=500000;
  maxCharts=-1;

  printf("Cover\n");fflush(stdout);
  A=MFCoverFlowOnManifold(M,F,"CoverTorus", n, u0, p0, Omega, eps, dt, tmax, 0.125, maxInterp, maxCharts, e);
  printf("done Cover\n");fflush(stdout);

  MFAtlasPageOutAllCharts(A,e);
  MFAtlasClosePlotfile(A,e);

  return 0;
 }

/* ----------------------------------------------------------------------------------------------- */

void IsoFlow(double *u, double *p, double *F,  void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"IsoFlow"};
  struct IsoFlowData *flow;

  double g;
  double gradg[3];
  double w[3];
  double s;

  flow=(struct IsoFlowData*)data;

  g=IMFIsoSurfaceFlowEvaluateG((struct IsoFlowData *)data,u);
  IMFIsoSurfaceFlowEvaluateDG ((struct IsoFlowData *)data,u,gradg);
  IMFIsoSurfaceFlowEvaluateW  ((struct IsoFlowData *)data,u,w);

  s=g+(gradg[0]*w[0]+gradg[1]*w[1]+gradg[2]*w[2])/(gradg[0]*gradg[0]+gradg[1]*gradg[1]+gradg[2]*gradg[2]);

  F[0] = w[0] - s * gradg[0];
  F[1] = w[1] - s * gradg[1];
  F[2] = w[2] - s * gradg[2];

  if(0){printf("IsoFlow: u=(%lf,%lf,%lf), g=%lf, s=%lf, gradg=(%lf,%lf,%lf), w=(%lf,%lf,%lf), F=(%lf,%lf,%lf)\n",u[0],u[1],u[2],g,s,gradg[0],gradg[1],gradg[2],w[0],w[1],w[2],F[0],F[1],F[2]);fflush(stdout);}

  if(gradg[0]!=gradg[0])exit(12);
  if(gradg[1]!=gradg[1])exit(12);
  if(gradg[2]!=gradg[2])exit(12);

  return;
 }

void dIsoFlow(double *u, double *p, double *dF,  void *data, MFErrorHandler e)
 {
  static char RoutineName[]={"dIsoFlow"};

  double g;
  double gradg[3];
  double *dgradg[3];
  double w[3];
  double *dw[3];
  double s,sx,sy,sz;

  dgradg[0]=(double*)malloc(3*sizeof(double));
  dgradg[1]=(double*)malloc(3*sizeof(double));
  dgradg[2]=(double*)malloc(3*sizeof(double));

  dw[0]=(double*)malloc(3*sizeof(double));
  dw[1]=(double*)malloc(3*sizeof(double));
  dw[2]=(double*)malloc(3*sizeof(double));

  g=IMFIsoSurfaceFlowEvaluateG((struct IsoFlowData *)data,u);
  IMFIsoSurfaceFlowEvaluateDG ((struct IsoFlowData *)data,u,gradg);
  IMFIsoSurfaceFlowEvaluateDDG((struct IsoFlowData *)data,u,dgradg);
  IMFIsoSurfaceFlowEvaluateW  ((struct IsoFlowData *)data,u,w);
  IMFIsoSurfaceFlowEvaluateDW ((struct IsoFlowData *)data,u,dw);

  s=g+(gradg[0]*w[0]+gradg[1]*w[1]+gradg[2]*w[2])/(gradg[0]*gradg[0]+gradg[1]*gradg[1]+gradg[2]*gradg[2]);

  sx=gradg[0]
    +(dgradg[0][0]* w[0]   +dgradg[1][0]* w[1]   +dgradg[2][0]* w[2]
     + gradg[0]   *dw[0][0]+ gradg[1]   *dw[1][0]+ gradg[2]   *dw[2][0])
                                                /(gradg[0]*gradg[0]+gradg[1]*gradg[1]+gradg[2]*gradg[2])
    -2*(gradg[0]*w[0]+gradg[1]*w[1]+gradg[2]*w[2])*(gradg[0]*dgradg[0][0]+gradg[1]*dgradg[1][0]+gradg[2]*dgradg[2][0])
                                                /(gradg[0]*gradg[0]+gradg[1]*gradg[1]+gradg[2]*gradg[2])
                                                /(gradg[0]*gradg[0]+gradg[1]*gradg[1]+gradg[2]*gradg[2]);

  sy=gradg[1]
    +(dgradg[0][1]* w[0]   +dgradg[1][1]* w[1]   +dgradg[2][1]* w[2]
     + gradg[0]   *dw[0][1]+ gradg[1]   *dw[1][1]+ gradg[2]   *dw[2][1])
                                                /(gradg[0]*gradg[0]+gradg[1]*gradg[1]+gradg[2]*gradg[2])
    -2*(gradg[0]*w[0]+gradg[1]*w[1]+gradg[2]*w[2])*(gradg[0]*dgradg[0][1]+gradg[1]*dgradg[1][1]+gradg[2]*dgradg[2][1])
                                                /(gradg[0]*gradg[0]+gradg[1]*gradg[1]+gradg[2]*gradg[2])
                                                /(gradg[0]*gradg[0]+gradg[1]*gradg[1]+gradg[2]*gradg[2]);

  sz=gradg[2]
    +(dgradg[0][2]* w[0]   +dgradg[1][2]* w[1]   +dgradg[2][2]* w[2]
     + gradg[0]   *dw[0][2]+ gradg[1]   *dw[1][2]+ gradg[2]   *dw[2][2])
                                                /(gradg[0]*gradg[0]+gradg[1]*gradg[1]+gradg[2]*gradg[2])
    -2*(gradg[0]*w[0]+gradg[1]*w[1]+gradg[2]*w[2])*(gradg[0]*dgradg[0][2]+gradg[1]*dgradg[1][2]+gradg[2]*dgradg[2][2])
                                                /(gradg[0]*gradg[0]+gradg[1]*gradg[1]+gradg[2]*gradg[2])
                                                /(gradg[0]*gradg[0]+gradg[1]*gradg[1]+gradg[2]*gradg[2]);
  dF[0+3*0]=dgradg[0][0] - sx * gradg[0] - s * dgradg[0][0];
  dF[0+3*1]=dgradg[0][1] - sy * gradg[0] - s * dgradg[0][1];
  dF[0+3*2]=dgradg[0][2] - sz * gradg[0] - s * dgradg[0][2];

  dF[1+3*0]=dgradg[1][0] - sx * gradg[1] - s * dgradg[1][0];
  dF[1+3*1]=dgradg[1][1] - sy * gradg[1] - s * dgradg[1][1];
  dF[1+3*2]=dgradg[1][2] - sz * gradg[1] - s * dgradg[1][2];

  dF[2+3*0]=dgradg[2][0] - sx * gradg[2] - s * dgradg[2][0];
  dF[2+3*1]=dgradg[2][1] - sy * gradg[2] - s * dgradg[2][1];
  dF[2+3*2]=dgradg[2][2] - sz * gradg[2] - s * dgradg[2][2];

  free(dgradg[0]);
  free(dgradg[1]);
  free(dgradg[2]);

  free(dw[0]);
  free(dw[1]);
  free(dw[2]);

  return;
 }

struct IsoFlowData
 {
  double (*g)(double*);
  void (*dg)(double*,double*);
  void (*ddg)(double*,double**);
  void (*w)(double*,double*);
  void (*dw)(double*,double**);
 };

double IMFIsoSurfaceFlowEvaluateG(struct IsoFlowData *f, double *u)
 {
  if(0){printf("Evaluate IsoSurfaceFlow g=0x%8.8x\n",f->g);fflush(stdout);}
  return f->g(u);
 }

void   IMFIsoSurfaceFlowEvaluateDG(struct IsoFlowData *f, double*u, double *gradg)
 {
  f->dg(u,gradg);

  return;
 }

void   IMFIsoSurfaceFlowEvaluateDDG(struct IsoFlowData *f, double*u, double **dgradg)
 {
  f->ddg(u,dgradg);

  return;
 }

void   IMFIsoSurfaceFlowEvaluateW(struct IsoFlowData *f, double*u, double *w)
 {
  f->w(u,w);

  return;
 }

void   IMFIsoSurfaceFlowEvaluateDW(struct IsoFlowData *f, double*u, double **dw)
 {
  f->dw(u,dw);

  return;
 }

/* ----------------------------------------------------------------------------------------------- */

IMFFlow IMFCreateIosurfaceFlow(double (*g)(double*),void (*dg)(double*,double*),void (*ddg)(double*,double**),void (*w)(double*,double*),void (*dw)(double*,double**),MFErrorHandler e)
 {
  static char RoutineName[]={"IMFCreateIosurfaceFlow"};
  struct IsoFlowData *data;

  data=(struct IsoFlowData*)malloc(sizeof(struct IsoFlowData));

  data->g  =  g;
  data->dg = dg;
  data->ddg=ddg;
  data->w  =  w;
  data->dw = dw;

  return IMFCreateFlow(3,0,IsoFlow,dIsoFlow,NULL,NULL,NULL,(void*)data,NULL,e);
 }

/* ----------------------------------------------------------------------------------------------- */

double gTorus(double *X)
 {
  double result;
  double x,y,z;

  x=X[0]; y=X[1]; z=X[2];

  result=(sqrt(x*x+y*y)-Rxy)*(sqrt(x*x+y*y)-Rxy)+z*z-Rz*Rz;

  return result;
 }

void dgTorus(double *X, double *dg)
 {
  double x,y,z;

  x=X[0]; y=X[1]; z=X[2];

/*
 *   (sqrt(x*x+y*y)-Rxy)*(sqrt(x*x+y*y)-Rxy)+z*z-Rz*Rz;
 */

  dg[0]=2*x*(1.-Rxy/sqrt(x*x+y*y) );
  dg[1]=2*y*(1.-Rxy/sqrt(x*x+y*y) );
  dg[2]=2*z;

/*
  {
   double eps=1.e-4;
   double t;
   double gp,gm;
   int i;

   for(i=0;i<3;i++)
    {
     t=X[i];

     X[i]=t+eps;
     gp=gTorus(X);

     X[i]=t-eps;
     gm=gTorus(X);
     X[i]=t;

     if(fabs(dg[i]-(gp-gm)/2/eps) > eps)
      {
       printf("error in dg[%d]=%10.7le, approx=%10.7le, error=%14.7le\n",i,dg[i],(gp-gm)/eps/2.,fabs(dg[i]-(gp-gm)/2/eps));fflush(stdout);
      }
    }
  }

 */
 
  return;
 }

void ddgTorus(double *X, double **ddg)
 {
  double x,y,z;

  x=X[0]; y=X[1]; z=X[2];

  ddg[0][0]=2*  (1.-Rxy/sqrt(x*x+y*y) ) + 2*x*x*Rxy/pow(x*x+y*y,1.5);
  ddg[0][1]=2*x*y*Rxy/pow(x*x+y*y,1.5);
  ddg[0][2]=0.;

  ddg[1][0]=2*x*y*Rxy/pow(x*x+y*y,1.5);
  ddg[1][1]=2*  (1.-Rxy/sqrt(x*x+y*y) ) + 2*y*y*Rxy/pow(x*x+y*y,1.5);
  ddg[1][2]=0.;

  ddg[2][0]=0.;
  ddg[2][1]=0.;
  ddg[2][2]=2.;

/*
  {
   double eps=1.e-4;
   double t;
   double gp[3],gm[3];
   int i,j;

   for(i=0;i<3;i++)
    {
     t=X[i];

     X[i]=t+eps;
     dgTorus(X,gp);

     X[i]=t-eps;
     dgTorus(X,gm);
     X[i]=t;

     for(j=0;j<3;j++)
      {
       if(fabs(ddg[i][j]-(gp[j]-gm[j])/2/eps)>eps)
        {
         printf("error in ddg[%d,%d]=%14.7le, approx=%14.7le, error=%14.7le\n",i,j,ddg[i][j],(gp[j]-gm[j])/eps/2.,fabs(ddg[i][j]-(gp[j]-gm[j])/2/eps));fflush(stdout);
        }
      }
    }
  }
 */

  return;
 }

#define DIRECTION -1

void w(double *x, double *w)
 {
  w[0]=DIRECTION*(x[2]);
  w[1]=DIRECTION*(x[0]-x[1]);
  w[2]=DIRECTION*(x[2]+x[1]);

  return;
 }

void dw(double *x, double **dw)
 {
  dw[0][0]=DIRECTION*( 0.);
  dw[0][1]=DIRECTION*( 0.);
  dw[0][2]=DIRECTION*( 1.);

  dw[1][0]=DIRECTION*( 1.);
  dw[1][1]=DIRECTION*(-1.);
  dw[1][2]=DIRECTION*( 0.);

  dw[2][0]=DIRECTION*( 0.);
  dw[2][1]=DIRECTION*( 1.);
  dw[2][2]=DIRECTION*( 1.);

/*
  {
   double eps=1.e-4;
   double t;
   double wp[3],wm[3];
   double dwA[3];
   int i,j;

   for(i=0;i<3;i++)
    {
     t=x[i];

     x[i]=t+eps;
     w(x,wp);

     x[i]=t-eps;
     w(x,wm);
     x[i]=t;

     dwA[0]=(wp[0]-wm[0])/eps/2;
     dwA[1]=(wp[1]-wm[1])/eps/2;
     dwA[2]=(wp[2]-wm[2])/eps/2;

     for(j=0;j<3;j++)
      {
       if(fabs(dw[j][i]-dwA[j])>eps)
        {
         printf("error in dw[%d][%d]=%14.7le, approx=%14.7le, error=%14.7le\n",j,i,dw[j][i],dwA[j],fabs(dw[j][i]-dwA[j]));fflush(stdout);
        }
      }
    }
  }
 */

  return;
 }

void findFixedPoint(double *x0,double *y, MFErrorHandler e)
 {
  double J[16];
  double f[4];

  double v[3];
  double *dv[3];
  double g;
  double dg[3];
  double *ddg[3];

  double x[4];
  double err;
  int    itimes;
  int    verbose=0;

/* Fixed Points: 
 *
 *     g(x,y,z)=0
 *     dg(x,y,z)=a*v(x,y,z)
 *
 */

  dv[0]=(double*)malloc(3*sizeof(double));
  dv[1]=(double*)malloc(3*sizeof(double));
  dv[2]=(double*)malloc(3*sizeof(double));
  ddg[0]=(double*)malloc(3*sizeof(double));
  ddg[1]=(double*)malloc(3*sizeof(double));
  ddg[2]=(double*)malloc(3*sizeof(double));

         w(x0,  v);
        dw(x0, dv);
  g=gTorus(x0    );
   dgTorus(x0, dg);
  ddgTorus(x0,ddg);

  x[0]=x0[0];
  x[1]=x0[1];
  x[2]=x0[2];
  if(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]>1.e-4)
    x[3]=(dg[0]*v[0]+ dg[1]*v[1]+ dg[2]*v[2])/(v[0]*v[0]+ v[1]*v[1]+ v[2]*v[2]);
   else
    x[3]=1;

   f[0]=g;
   f[1]=dg[0]-x[3]*v[0];
   f[2]=dg[1]-x[3]*v[1];
   f[3]=dg[2]-x[3]*v[2];
   err=sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]+f[3]*f[3]);

/*
   {
    double t;
    double f0p,f1p,f2p,f3p;
    double f0m,f1m,f2m,f3m;
    double J0,J1,J2,J3;
    double eps=1.e-4;
    int i;

           w(x,  v);
          dw(x, dv);
    g=gTorus(x    );
     dgTorus(x, dg);
    ddgTorus(x,ddg);
 
    J[0+4*0]= dg   [0];               J[0+4*1]= dg   [1];               J[0+4*2]= dg   [2];               J[0+4*3]=   0.;
    J[1+4*0]=ddg[0][0]-x[3]*dv[0][0]; J[1+4*1]=ddg[0][1]-x[3]*dv[0][1]; J[1+4*2]=ddg[0][2]-x[3]*dv[0][2]; J[1+4*3]=-v[0];
    J[2+4*0]=ddg[1][0]-x[3]*dv[1][0]; J[2+4*1]=ddg[1][1]-x[3]*dv[1][1]; J[2+4*2]=ddg[1][2]-x[3]*dv[1][2]; J[2+4*3]=-v[1];
    J[3+4*0]=ddg[2][0]-x[3]*dv[2][0]; J[3+4*1]=ddg[2][1]-x[3]*dv[2][1]; J[3+4*2]=ddg[2][2]-x[3]*dv[2][2]; J[3+4*3]=-v[2];

    for(i=0;i<4;i++)
     {
      t=x[i];

      x[i]=t+eps;

             w(x,  v);
            dw(x, dv);
      g=gTorus(x    );
       dgTorus(x, dg);
      ddgTorus(x,ddg);

      f0p=g;
      f1p=dg[0]-x[3]*v[0];
      f2p=dg[1]-x[3]*v[1];
      f3p=dg[2]-x[3]*v[2];
  
      x[i]=t-eps;

             w(x,  v);
            dw(x, dv);
      g=gTorus(x    );
       dgTorus(x, dg);
      ddgTorus(x,ddg);

      f0m=g;
      f1m=dg[0]-x[3]*v[0];
      f2m=dg[1]-x[3]*v[1];
      f3m=dg[2]-x[3]*v[2];

      J0=.5*(f0p-f0m)/eps;
      J1=.5*(f1p-f1m)/eps;
      J2=.5*(f2p-f2m)/eps;
      J3=.5*(f3p-f3m)/eps;
 
      printf("J[0,%d]=%14.7le, approx=%14.7le, diff=%14.7le\n",i,J[0+4*i],J0,fabs(J[0+4*i]-J0));fflush(stdout);
      printf("J[1,%d]=%14.7le, approx=%14.7le, diff=%14.7le\n",i,J[1+4*i],J1,fabs(J[1+4*i]-J1));fflush(stdout);
      printf("J[2,%d]=%14.7le, approx=%14.7le, diff=%14.7le\n",i,J[2+4*i],J2,fabs(J[2+4*i]-J2));fflush(stdout);
      printf("J[3,%d]=%14.7le, approx=%14.7le, diff=%14.7le\n",i,J[3+4*i],J3,fabs(J[3+4*i]-J3));fflush(stdout);
 
      x[i]=t;
     }
   }
 */

  itimes=0;
  while(err>1.e-7 && itimes<20)
   {
    if(verbose){printf("%d %le f=(%lf,%lf,%lf,%lf)\n",itimes,err,f[0],f[1],f[2],f[3]);fflush(stdout);}

    J[0+4*0]= dg   [0];               J[0+4*1]= dg   [1];               J[0+4*2]= dg   [2];               J[0+4*3]=   0.;
    J[1+4*0]=ddg[0][0]-x[3]*dv[0][0]; J[1+4*1]=ddg[0][1]-x[3]*dv[0][1]; J[1+4*2]=ddg[0][2]-x[3]*dv[0][2]; J[1+4*3]=-v[0];
    J[2+4*0]=ddg[1][0]-x[3]*dv[1][0]; J[2+4*1]=ddg[1][1]-x[3]*dv[1][1]; J[2+4*2]=ddg[1][2]-x[3]*dv[1][2]; J[2+4*3]=-v[1];
    J[3+4*0]=ddg[2][0]-x[3]*dv[2][0]; J[3+4*1]=ddg[2][1]-x[3]*dv[2][1]; J[3+4*2]=ddg[2][2]-x[3]*dv[2][2]; J[3+4*3]=-v[2];

    MFSolveFull(4,J,f,e);

    x[0]=x[0]-f[0];
    x[1]=x[1]-f[1];
    x[2]=x[2]-f[2];
    x[3]=x[3]-f[3];

           w(x,  v);
          dw(x, dv);
    g=gTorus(x    );
     dgTorus(x, dg);
    ddgTorus(x,ddg);
    
    f[0]=g;
    f[1]=dg[0]-x[3]*v[0];
    f[2]=dg[1]-x[3]*v[1];
    f[3]=dg[2]-x[3]*v[2];
    err=sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]+f[3]*f[3]);

    itimes++;
   }

  y[0]=x[0];
  y[1]=x[1];
  y[2]=x[2];

  free(dv[0]);
  free(dv[1]);
  free(dv[2]);
  free(ddg[0]);
  free(ddg[1]);
  free(ddg[2]);

  return;
 }
