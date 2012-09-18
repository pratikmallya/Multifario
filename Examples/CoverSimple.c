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

  N=2;
  for(i=0;i<N;i++)
   {
    u0[n]=MFIMFVectorFactory(M,e);
    theta=2.*3.151926*0.;
    phi  =2.*3.151926*i/(N-1.);
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

  p0=(MFKVector)NULL;

  Omega=MFNRegionCreateHyperCube(3,1.1,e);

  eps=.4;
  dt=.01;
  tmax= 2000.;
  maxInterp=500000;
  maxCharts=-1;

  printf("Cover\n");fflush(stdout);
  A=MFCoverFlowOnManifold(M,F,"CoverSimple", n, u0, p0, Omega, eps, dt, tmax, 0.125, maxInterp, maxCharts, e);
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

void w(double *x, double *w)
 {
  double theta,phi;

  theta=atan2(x[1],x[0]);
  phi  =asin(x[2]/Rz);

/*
  x=cos(theta)*(Rxy+Rz*cos(phi));
  y=sin(theta)*(Rxy+Rz*cos(phi));
  z=Rz*sin(phi);

  w= dx/dtheta + .5*dx/dphi;
 */

  w[0]=-sin(theta)*(Rxy+Rz*cos(phi));
/* -.5*cos(theta)*Rz*sin(phi); */
  w[1]= cos(theta)*(Rxy+Rz*cos(phi));
/* -.5*sin(theta)*Rz*sin(phi); */
  w[2]= 0.;
                           /* .5*           Rz*cos(phi); */

  return;
 }

void dw(double *x, double **dw)
 {
  double theta,phi;
  double dthetadx,dphidx;
  double dthetady,dphidy;
  double dthetadz,dphidz;
  double dw0dtheta,dw0dphi;
  double dw1dtheta,dw1dphi;
  double dw2dtheta,dw2dphi;

  theta=atan2(x[1],x[0]);
  phi  =asin(x[2]/Rz);

  dthetadx=-cos(theta)*cos(theta)*x[1]/x[0]/x[0];
  dthetady= cos(theta)*cos(theta)/x[0];
  dthetadz= 0.;

  dphidx=0.;
  dphidy=0.;
  dphidz=1./Rz/cos(phi);

  dw0dtheta=-cos(theta)*(Rxy+Rz*cos(phi))+.5*sin(theta)*Rz*sin(phi);
  dw1dtheta= sin(theta)*(Rxy+Rz*cos(phi))-.5*cos(theta)*Rz*sin(phi);
  dw2dtheta= 0.;

  dw0dphi= sin(theta)*Rz*sin(phi)-.5*cos(theta)*Rz*cos(phi);
  dw1dphi=-cos(theta)*Rz*sin(phi)-.5*sin(theta)*Rz*cos(phi);
  dw2dphi=                       -.5*           Rz*sin(phi);
  dw0dphi= 0.;
  dw1dphi= 0.;
  dw2dphi= 0.;

  dw[0][0]=dw0dtheta*dthetadx+dw0dphi*dphidx;
  dw[1][0]=dw1dtheta*dthetadx+dw1dphi*dphidx;
  dw[2][0]=dw2dtheta*dthetadx+dw2dphi*dphidx;

  dw[0][1]=dw0dtheta*dthetady+dw0dphi*dphidy;
  dw[1][1]=dw1dtheta*dthetady+dw1dphi*dphidy;
  dw[2][1]=dw2dtheta*dthetady+dw2dphi*dphidy;

  dw[0][2]=dw0dtheta*dthetadz+dw0dphi*dphidz;
  dw[1][2]=dw1dtheta*dthetadz+dw1dphi*dphidz;
  dw[2][2]=dw2dtheta*dthetadz+dw2dphi*dphidz;

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

  return;
 }
