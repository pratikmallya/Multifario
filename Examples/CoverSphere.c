#include <MFAtlas.h>

struct IsoFlowData;

double IMFIsoSurfaceFlowEvaluateG  (struct IsoFlowData *f,double *u);
void   IMFIsoSurfaceFlowEvaluateDG (struct IsoFlowData *f,double *u,double *gradg);
void   IMFIsoSurfaceFlowEvaluateDDG(struct IsoFlowData *f,double *u,double **dgradg);
void   IMFIsoSurfaceFlowEvaluateW  (struct IsoFlowData *f,double *u,double *w);
void   IMFIsoSurfaceFlowEvaluateDW (struct IsoFlowData *f,double *u,double **dw);

MFAtlas MFCoverFlowOnManifold(MFImplicitMF M, IMFFlow F,char *name, int nInitial, MFNVector *u0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, double dmin, int maxInterp, int maxCharts, MFErrorHandler e);

IMFFlow IMFCreateIsosurfaceFlow(double (*g)(double*),void (*dg)(double*,double*),void (*ddg)(double*,double**),void (*w)(double*,double*),void (*dw)(double*,double**),MFErrorHandler e);

double gGenusTwo(double *x);
void dgGenusTwo(double *x, double *dg);
void ddgGenusTwo(double *x, double **dg);

double gSphere(double *x);
void dgSphere(double *x, double *dg);
void ddgSphere(double *x, double **dg);

void w(double *x, double *w);
void dw(double *x, double **dw);

void F(int *n,double *z,int *m,double *f, void *data, MFErrorHandler e)
 {
  f[0]=gSphere(z);
  return;
 }

void dF(int *n,double *z,int *m,double *f, void *data, MFErrorHandler e)
 {
  dgSphere(z,f);

  return;
 }

void ddF(int *n,double *z,int *m,double *f, void *data, MFErrorHandler e)
 {
  static double *g[3]={NULL,NULL,NULL};

  g[0]=(double*)realloc(g[0],3*sizeof(double));
  g[1]=(double*)realloc(g[1],3*sizeof(double));
  g[2]=(double*)realloc(g[2],3*sizeof(double));

  ddgSphere(z,g);

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
  MFImplicitMF MFIMFCreateAlgebraicSubroutine(int,int,
                              void (*)(int*,double*,int*,double*,void*,MFErrorHandler),
                              void (*)(int*,double*,int*,double*,void*,MFErrorHandler),
                              void (*)(int*,double*,int*,double*,void*,MFErrorHandler),
                              void *,
                              MFErrorHandler);

/* ----------------------------------------------------------------------------------------------- */

int main(int argc, char *argv[])
 {
  MFAtlas A;
  MFImplicitMF M;
  IMFFlow F;
  char *name;
  int nInitial;
  MFNVector u0[100];
  MFKVector p0;
  MFNRegion Omega;
  double eps;
  double dt;
  double tmax;
  int maxInterp;
  int maxCharts;
  MFErrorHandler e;
  int i,n;
  double x,y,z;
  double theta,phi;

  e=MFCreateErrorHandler();

  F=IMFCreateIsosurfaceFlow(gSphere,dgSphere,ddgSphere,w,dw,e);
  
  M=MFIMFCreateAlgebraicSubroutine(3,2,F,dF,ddF,NULL,e);
  MFIMFSetR(M,.05,e);

  n=0;

  for(i=0;i<10;i++)
   {
    u0[n]=MFIMFVectorFactory(M,e);
    theta=2.*3.151926*0.25;
    phi  =   3.151926*(.5+.075-.15*i/9.);
    x=cos(phi)*cos(theta);
    y=cos(phi)*sin(theta);
    z=sin(phi);
    MFNVSetC(u0[n],0,x,e);
    MFNVSetC(u0[n],1,cos(3.1415926*.25)*y+sin(3.1415926*.25)*z,e);
    MFNVSetC(u0[n],2,cos(3.1415926*.25)*z-sin(3.1415926*.25)*y,e);
    n++;

    u0[n]=MFIMFVectorFactory(M,e);
    theta=2.*3.151926*.75;
    phi  =  -3.151926*(.5+.075-.15*i/9.);
    x=cos(phi)*cos(theta);
    y=cos(phi)*sin(theta);
    z=sin(phi);
    MFNVSetC(u0[n],0,x,e);
    MFNVSetC(u0[n],1,cos(3.1415926*.25)*y+sin(3.1415926*.25)*z,e);
    MFNVSetC(u0[n],2,cos(3.1415926*.25)*z-sin(3.1415926*.25)*y,e);
    n++;
   }

  p0=(MFKVector)NULL;

  Omega=MFNRegionCreateHyperCube(3,1.1,e);

  eps=0.001;
  dt=.01;
  tmax= 100.;
  maxInterp=100;
  maxCharts=-1;

/*A=MFCoverFlowOnManifold(M,F,"CoverSphere", n, u0, p0, Omega, eps, dt, tmax, 0.7, maxInterp, maxCharts, e);*/
  A=MFCoverFlowOnManifold(M,F,"CoverSphere", n, u0, p0, Omega, eps, dt, tmax, 0.65, maxInterp, maxCharts, e);

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

IMFFlow IMFCreateIsosurfaceFlow(double (*g)(double*),void (*dg)(double*,double*),void (*ddg)(double*,double**),void (*w)(double*,double*),void (*dw)(double*,double**),MFErrorHandler e)
 {
  static char RoutineName[]={"IMFCreateIsosurfaceFlow"};
  struct IsoFlowData *data;

  data=(struct IsoFlowData*)malloc(sizeof(struct IsoFlowData));

  data->g  =  g;
  data->dg = dg;
  data->ddg=ddg;
  data->w  =  w;
  data->dw = dw;

  return IMFCreateFlow(3,0,IsoFlow,dIsoFlow,NULL,NULL,NULL,(void*)data,NULL,e);
/*IMFFlow IMFCreateFlow(int,int,MFFlowFunction,MFFlowFunction,MFFlowFunction,MFFlowFunction,MFFlowFunction,void*,MFFlowFreeData,MFErrorHandler);*/
 }

/* ----------------------------------------------------------------------------------------------- */

double gGenusTwo(double *x)
 {
  double result;
  result=((sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)*(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)+x[2]*x[2] -.3*.3)*
         ((sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)*(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)+x[2]*x[2] -.3*.3)-0.05;
  if(0){printf("in Evaluate IsoSurfaceFlow x=(%lf,%lf,%lf) result=%lf\n",x[0],x[1],x[2],result);fflush(stdout);}

  if(result!=result)exit(12);

  return result;
 }

void dgGenusTwo(double *x, double *dg)
 {
  double s;
  double numer;
  double denom;

  numer=sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8;
  denom=sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1]);

  s= 2*((sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)*(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)+x[2]*x[2] -.3*.3);

  dg[0]=s*(x[0]-.8)*(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1]);
  dg[1]=s* x[1]    *(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1]);
  dg[2]=s*2*x[2];

  if(0){printf("in Evaluate dg x=(%lf,%lf,%lf), gradg=(%lf,%lf,%lf)\n",x[0],x[1],x[2],dg[0],dg[1],dg[2]);fflush(stdout);}
  if(0){printf("    dg[0]=%lf*(%lf-.8)*%lf/%lf=%lf\n",s,x[0],numer,denom,s*(x[0]-.8)*numer/denom);fflush(stdout);}
  if(0){printf("    dg[1]=%lf* %lf    *%lf/%lf=%lf\n",s,x[1],numer,denom,s*x[1]*numer/denom);fflush(stdout);}

  return;
 }

void ddgGenusTwo(double *x, double **ddg)
 {
  double  s;
  double ds[3];

   s   = 2*((sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)*(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)+x[2]*x[2] -.3*.3);

  ds[0]= 8* (sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)*(x[0]-.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1]);
  ds[1]= 8* (sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)* x[1]    /sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1]);
  ds[2]= 4*x[2];

  ddg[0][0]=ds[0]*(x[0]-.8)*(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])
           + s             *(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])
           + s   *(x[0]-.8)*(
                 -(x[0]-.8)*(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/pow ((x[0]-.8)*(x[0]-.8)+x[1]*x[1],1.5)
                 +(x[0]-.8)/(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])
                            );

  ddg[0][1]=ds[1]*(x[0]-.8)*(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])
           + s   *(x[0]-.8)*(
                 - x[1]    *(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/pow ((x[0]-.8)*(x[0]-.8)+x[1]*x[1],1.5)
                 + x[1]    /(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])
                            );

  ddg[0][2]=ds[2]*(x[0]-.8)*(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1]);


  ddg[1][0]=ds[0]* x[1]    *(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])
           + s   * x[1]    *(x[0]-.8)*(
                  -(x[0]-.8)*(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/pow ((x[0]-.8)*(x[0]-.8)+x[1]*x[1],1.5)
                  +(x[0]-.8)/(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])
                            );

  ddg[1][1]=ds[1]* x[1]    *(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])
           + s             *(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])
           + s   * x[1]    *(
                  -x[1]    *(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/pow ((x[0]-.8)*(x[0]-.8)+x[1]*x[1],1.5)
                  +x[1]    /(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])
                            );


  ddg[1][2]=ds[2]* x[1]    *(sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1])-.8*.8)/sqrt((x[0]-.8)*(x[0]-.8)+x[1]*x[1]);

  ddg[2][0]=ds[0]*2*x[2];
  ddg[2][1]=ds[1]*2*x[2];
  ddg[2][2]=ds[2]*2*x[2]+s*2;

  return;
 }

double gSphere(double *x)
 {
  double result;
  result=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-1.;
  if(0){printf("in Evaluate IsoSurfaceFlow x=(%lf,%lf,%lf) result=%lf\n",x[0],x[1],x[2],result);fflush(stdout);}

  if(result!=result)exit(12);

  return result;
 }

void dgSphere(double *x, double *dg)
 {

  dg[0]=2*x[0];
  dg[1]=2*x[1];
  dg[2]=2*x[2];

  if(0){printf("in Evaluate dg x=(%lf,%lf,%lf), gradg=(%lf,%lf,%lf)\n",x[0],x[1],x[2],dg[0],dg[1],dg[2]);fflush(stdout);}

  return;
 }

void ddgSphere(double *x, double **ddg)
 {
  ddg[0][0]=0.;
  ddg[0][1]=0.;
  ddg[0][2]=0.;


  ddg[1][0]=0.;
  ddg[1][1]=0.;
  ddg[1][2]=0.;

  ddg[2][0]=0.;
  ddg[2][1]=0.;
  ddg[2][2]=0.;

  return;
 }

#define DIRECTION -1

void w(double *x, double *w)
 {
  w[0]=DIRECTION*(x[2]-x[1]);
  w[1]=DIRECTION*(x[0]+x[1]);
  w[2]=DIRECTION*(x[2]);

  return;
 }

void dw(double *x, double **dw)
 {
  dw[0][0]=DIRECTION*( 0.);
  dw[0][1]=DIRECTION*(-1.);
  dw[0][2]=DIRECTION*( 1.);

  dw[1][0]=DIRECTION*( 1.);
  dw[1][1]=DIRECTION*( 1.);
  dw[1][2]=DIRECTION*( 0.);

  dw[2][0]=DIRECTION*( 0.);
  dw[2][1]=DIRECTION*( 0.);
  dw[2][2]=DIRECTION*( 1.);

  return;
 }

