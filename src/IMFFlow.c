/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 */

static char *id="@(#) $Id: IMFFlow.c,v 1.8 2011/07/21 17:42:46 mhender Exp $";

static char IMFFlowErrorMsg[256];

#include <IMFFlow.h>
#include <IMFExpansion.h>
#include <IMFExpansionSpace.h>
#include <IMFExpansionPt.h>
#include <MFAtlas.h>
#include <MFErrorHandler.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct IMFFlowSt {int nu;
                  int np;
                  MFFlowFunction F;
                  MFFlowFunction dF;
                  MFFlowFunction dp;
                  MFFlowFunction ddF;
                  MFFlowFunction dddF;
                  MFFlowFreeData freeData;
                  MFNVector (*nVectorFactory)(IMFFlow,MFErrorHandler);
                  MFKVector (*kVectorFactory)(IMFFlow,MFErrorHandler);
                  int nRefs;
                  int fat;
                  IMFFlow fatF;
                  int k;
                  int dir;
                  void *data;
                 };

void IMFIntegratorDTimeDerivatives(IMFFlow,double*,MFKVector,double*,MFErrorHandler);

static MFNVector MFFlowDefaultNVectorFactory(IMFFlow this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorFactory"};

  return MFCreateNVector(IMFFlowNU(this,e),e);
 }

static MFKVector MFFlowDefaultKVectorFactory(IMFFlow this, MFErrorHandler e)
 {
  static char RoutineName[]={"MFNVectorFactory"};

  return MFCreateKVector(IMFFlowNP(this,e),e);
 }

IMFFlow IMFCreateFlow(int nu, int np, MFFlowFunction F, MFFlowFunction dF, MFFlowFunction dFdp, MFFlowFunction ddF, MFFlowFunction dddF,void *data, MFFlowFreeData freeData, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFCreateFlow"};
  IMFFlow result;

  result=(struct IMFFlowSt*)malloc(sizeof(struct IMFFlowSt));

#ifndef MFNOSAFETYNET
  if(result==NULL)
   {
    sprintf(IMFFlowErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct IMFFlowSt));
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->nu=nu;
  result->np=np;

  result->F=F;
  result->dF=dF;
  result->dp=dFdp;
  result->ddF=ddF;
  result->dddF=dddF;
  result->data=data;
  result->freeData=freeData;
  result->nVectorFactory=MFFlowDefaultNVectorFactory;
  result->kVectorFactory=MFFlowDefaultKVectorFactory;
  result->dir=1;

  result->fat=0;
  result->fatF=NULL;
  result->k=0;
  result->nRefs=1;

  return result;
 }

IMFFlow IMFCreateFatFlow(IMFFlow F,int k, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFCreateFatFlow"};
  IMFFlow result;
  int nu,m;

  m=IMFFlowNU(F,e);
  nu=m+m*k+m*k*k;

  result=(struct IMFFlowSt*)malloc(sizeof(struct IMFFlowSt));

#ifndef MFNOSAFETYNET
  if(result==NULL)
   {
    sprintf(IMFFlowErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct IMFFlowSt));
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->nu=nu;
  result->np=IMFFlowNP(F,e);

  result->F=NULL;
  result->dF=NULL;
  result->dp=NULL;
  result->ddF=NULL;
  result->dddF=NULL;
  result->data=NULL;
  result->freeData=NULL;
  result->dir=1;

  result->fat=1;
  result->fatF=F;
  IMFRefFlow(F,e);
  result->k=k;

  result->nRefs=1;

  return result;
 }

IMFFlow IMFCreateBackwardFatFlow(IMFFlow F,int k, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFCreateBackwardFlow"};
  IMFFlow result;
  int nu,m;

  m=IMFFlowNU(F,e);
  nu=m+m*k+m*k*k;

  result=(struct IMFFlowSt*)malloc(sizeof(struct IMFFlowSt));

#ifndef MFNOSAFETYNET
  if(result==NULL)
   {
    sprintf(IMFFlowErrorMsg,"Out of memory, trying to allocate %d bytes",sizeof(struct IMFFlowSt));
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return NULL;
   }
#endif

  result->nu=nu;

  result->F=NULL;
  result->dF=NULL;
  result->dp=NULL;
  result->ddF=NULL;
  result->dddF=NULL;
  result->data=NULL;
  result->freeData=NULL;
  result->dir=-1;

  result->fat=1;
  result->fatF=F;
  IMFRefFlow(F,e);
  result->k=k;

  result->nRefs=1;

  return result;
 }

void IMFFreeFlow(IMFFlow F, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFFreeFlow"};
  F->nRefs--;

  if(F->nRefs<=0)
   {
    if(F->freeData!=NULL)(F->freeData)(F->data,e);
    if(!F->fat)
     {
      F->nu=0;
      F->np=0;
      F->F=NULL;
      F->dF=NULL;
      F->dp=NULL;
      F->ddF=NULL;
      F->dddF=NULL;
      F->data=NULL;
      F->freeData=NULL;
      F->fat=0;
      F->fatF=NULL;
      F->k=0;
     }else{
      F->nu=0;
      F->np=0;
      F->F=NULL;
      F->dF=NULL;
      F->dp=NULL;
      F->ddF=NULL;
      F->dddF=NULL;
      F->data=NULL;
      F->freeData=NULL;
      F->fat=0;
      IMFFreeFlow(F->fatF,e);
      F->fatF=NULL;
      F->k=0;
     }

    free(F);
   }

  return;
 }

void IMFRefFlow(IMFFlow F, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFRefFlow"};

  F->nRefs++;
  return;
 }

void IMFEvaluateFlow(IMFFlow F, MFNVector vu, MFKVector vp, double *f, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFEvaluateFlow"};
  double *u;
  double *p;
  int i,n;

  if(vp!=NULL)p=MFKV_CStar(vp,e);
   else p=NULL;

  if(!strcmp(MFNVGetId(vu,e),"DENSE"))
   {
    u=MFNV_CStar(vu,e);
   }else if(!strcmp(MFNVGetId(vu,e),"IMFExpansionVector"))
   {
    u=IMFExpansionU(IMFExpansionNVGetE(vu,e),e);
   }else{
    sprintf(IMFFlowErrorMsg,"You cannot evaluate a flow at a point which is neither DENSE or and Expansion.");
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(!F->fat)
   {
    F->F(u,p,f,F->data,e);
   }else{
    int k,m;
    int j,l,p,q,r,v;
    double d;
    double uFuu;
    static double *uFu=NULL;
    static double *sdf=NULL;
    static double *sddf=NULL;

    n=(F->fatF)->nu;
    k=F->k;
    m=n+n*k+n*k*k;

    sdf=(double*)realloc((void*)sdf,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
    if(sdf==NULL)
     {
      sprintf(IMFFlowErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
      MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    sddf=(double*)realloc((void*)sddf,n*n*n*sizeof(double));

#ifndef MFNOSAFETYNET
    if(sddf==NULL)
     {
      sprintf(IMFFlowErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*n*sizeof(double));
      MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    uFu=(double*)realloc((void*)uFu,k*k*sizeof(double));

#ifndef MFNOSAFETYNET
    if(uFu==NULL)
     {
      sprintf(IMFFlowErrorMsg,"Out of memory, trying to allocate %d bytes",k*k*sizeof(double));
      MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    IMFEvaluateFlow(F->fatF,vu,vp,f,e);
    IMFEvaluateDerivativeOfFlow(F->fatF,vu,vp,sdf,e);
    IMFEvaluateSecondDerivativeOfFlow(F->fatF,vu,vp,sddf,e);

    d=0.;for(i=0;i<n;i++)d+=f[i]*f[i];d=sqrt(d+1.);
    for(i=0;i<n;i++)f[i]=f[i]/d;

    for(r=0;r<k;r++)
     for(j=0;j<k;j++)
      {
       uFu[r+k*j]=0.;
       for(p=0;p<n;p++)
         for(q=0;q<n;q++)
           uFu[r+k*j]+=u[n+p+n*r]*sdf[p+n*q]*u[n+q+n*j];
      }

    for(i=0;i<n*k;i++)f[n+i]=0.;

    for(j=0;j<k;j++)
     {
      for(i=0;i<n;i++)
       {
        for(p=0;p<n;p++)
          f[n+i+n*j]+=sdf[i+n*p]*u[n+p+n*j];

        for(r=0;r<k;r++)
          f[n+i+n*j]-=uFu[r+k*j]*u[n+i+n*r];

        f[n+i+n*j]=f[n+i+n*j]/d;
       }
     }

    for(i=0;i<n*k*k;i++)f[n+n*k+i]=0.;
#ifdef DOSECND
    for(l=0;l<k;l++)
     {
      for(j=0;j<k;j++)
       {
        for(i=0;i<n;i++)
         {
          for(p=0;p<n;p++)
            f[n+n*k+i+n*(j+k*l)]+=sdf[i+n*p]*u[n+n*k+p+n*(j+k*l)];
          for(p=0;p<n;p++)
            for(q=0;q<n;q++)
              f[n+n*k+i+n*(j+k*l)]+=sddf[i+n*(p+n*q)]*u[n+p+n*j]*u[n+q+n*l];
          for(r=0;r<k;r++)
            f[n+n*k+i+n*(j+k*l)]-=uFu[r+k*j]*u[n+n*k+i+n*(r+k*l)]+uFu[r+k*l]*u[n+n*k+i+n*(r+k*j)];

          for(r=0;r<k;r++)
           {
            uFuu=0.;
            for(p=0;p<n;p++)
              for(q=0;q<n;q++)
               {
                uFuu+=u[n+p+n*r]*sdf[p+n*q]*u[n+n*k+q+n*(j+k*l)];
                uFuu+=u[n+p+n*r]*sdf[q+n*p]*u[n+n*k+q+n*(j+k*l)];
                for(v=0;v<n;v++)
                  uFuu+=u[n+p+n*r]*sddf[p+n*(q+n*v)]*u[n+q+n*j]*u[n+v+n*l];
               }
            f[n+n*k+i+n*(j+k*l)]-=uFuu*u[i+n*r];
           }

          f[n+n*k+i+n*(j+k*l)]=1.*f[n+n*k+i+n*(j+k*l)]/d;
         }
       }
     }
#endif
   }

  n=F->nu;
  if(F->dir==-1)for(i=0;i<n;i++)f[i]=-f[i];

  return;
 }

void IMFEvaluateDerivativeOfFlow(IMFFlow F, MFNVector vu, MFKVector vp, double *df, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFEvaluateDerivativeOfFlow"};
  double *u;
  double *p;
  int i,n;

  if(vp!=NULL)p=MFKV_CStar(vp,e);
   else p=NULL;

  if(!strcmp(MFNVGetId(vu,e),"DENSE"))
   {
    u=MFNV_CStar(vu,e);
   }else if(!strcmp(MFNVGetId(vu,e),"IMFExpansionVector"))
   {
    u=IMFExpansionU(IMFExpansionNVGetE(vu,e),e);
   }else{
    sprintf(IMFFlowErrorMsg,"You cannot evaluate a flow at a point which is neither DENSE or and Expansion.");
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(!F->fat)
   {
    F->dF(u,p,df,F->data,e);
   }else if(1)
   {
    IMFIntegratorDTimeDerivatives(F,u,vp,df,e);
   }else{
    int i,j,l,n,k,m;
    static double *sdf=NULL;
    static double *sf=NULL;
    double d;

    n=(F->fatF)->nu;
    k=(F->fatF)->k;
    m=n+n*k+n*k*k;

    sdf=(double*)realloc((void*)sdf,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
    if(sdf==NULL)
     {
      sprintf(IMFFlowErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
      MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    sf=(double*)realloc((void*)sf,n*sizeof(double));

#ifndef MFNOSAFETYNET
    if(sf==NULL)
     {
      sprintf(IMFFlowErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
      MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
      MFErrorHandlerOutOfMemory(e);
      return;
     }
#endif

    IMFEvaluateFlow(F->fatF,vu,vp,sf,e);
    IMFEvaluateDerivativeOfFlow(F->fatF,vu,vp,sdf,e);

    d=0.;for(i=0;i<n;i++)d+=sf[i]*sf[i];d=sqrt(d);

    for(i=0;i<m*m;i++)df[i]=0.;
    for(i=0;i<n;i++)
      for(j=0;j<n;j++)df[i+m*j]=sdf[i+n*j]/d;
    for(i=0;i<n;i++)
      for(j=0;j<n;j++)
        for(l=0;l<n;l++)
          df[i+m*j]-=sf[i]*sf[l]*sdf[l+n*j]/d/d/d;
   }

  n=F->nu;
  if(F->dir==-1)for(i=0;i<n*n;i++)df[i]=-df[i];

  return;
 }

/*!\fn void IMFEvaluateParameterDerivativeOfFlow(IMFFlow F, MFNVector vu, MFKVector vp,double *dp, MFErrorHandler e)
 * \brief Evaluates the derivative of the flow at a point in phase space with respect to the parameter space variables (dFdp).
 *
 * \param F The flow.
 * \param vu A point in phase space.
 * \param vp A point in parameter space.
 * \param dp An array of length at least the dimension of the phase space times the dimension of the parameter space,
 *           to hold the derivatives. It is
 *           filled with entries in column order. That is dF^i_j = dp[i+nu*j], where nu is the dimension
 *           of the phase space.
 */
void IMFEvaluateParameterDerivativeOfFlow(IMFFlow F, MFNVector vu, MFKVector vp,double *dp, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFEvaluateParameterDerivativeOfFlow"};
  double *u;
  double *p;
  int i;
  int nu,np;

  if(vp!=NULL)p=MFKV_CStar(vp,e);
   else p=NULL;

  if(!strcmp(MFNVGetId(vu,e),"DENSE"))
   {
    u=MFNV_CStar(vu,e);
   }else if(!strcmp(MFNVGetId(vu,e),"IMFExpansionVector"))
   {
    u=IMFExpansionU(IMFExpansionNVGetE(vu,e),e);
   }else{
    sprintf(IMFFlowErrorMsg,"You cannot evaluate a flow at a point which is neither DENSE or and Expansion.");
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(!F->fat)
   {
    F->dp(u,p,dp,F->data,e);
   }else{
    sprintf(IMFFlowErrorMsg,"You cannot do this for a fat flow (not implemented).");
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    return;
   }

  nu=F->nu;
  np=F->np;
  if(F->dir==-1)for(i=0;i<nu*np;i++)dp[i]=-dp[i];

  return;
 }

/*!\fn void IMFEvaluateSecondDerivativeOfFlow(IMFFlow F, MFNVector vu, MFKVector vp,double *ddf, MFErrorHandler e)
 * \brief Evaluates the second derivative of the flow at a point in phase space (dFdudu).
 *
 * \param F The flow.
 * \param vu A point in phase space.
 * \param vp A point in parameter space.
 * \param ddf An array of length at least the cube of the dimension of the phase space to hold the derivatives. It is
 *           filled with entries in column order. That is dF^i_jk = df[i+nu*(j+nu*k)], where nu is the dimension of the phase space.
 */
void IMFEvaluateSecondDerivativeOfFlow(IMFFlow F, MFNVector vu, MFKVector vp,double *ddf, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFEvaluateSecondDerivativeOfFlow"};
  double *u;
  double *p;
  int i,n;

  if(vp!=NULL)p=MFKV_CStar(vp,e);
   else p=NULL;

  if(!strcmp(MFNVGetId(vu,e),"DENSE"))
   {
    u=MFNV_CStar(vu,e);
   }else if(!strcmp(MFNVGetId(vu,e),"IMFExpansionVector"))
   {
    u=IMFExpansionU(IMFExpansionNVGetE(vu,e),e);
   }else{
    sprintf(IMFFlowErrorMsg,"You cannot evaluate a flow at a point which is neither DENSE or and Expansion.");
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(!F->fat)
   {
    F->ddF(u,p,ddf,F->data,e);
   }else{
    sprintf(IMFFlowErrorMsg,"You cannot do this for a fat flow (not implemented).");
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    return;
   }

  n=F->nu;
  if(F->dir==-1)for(i=0;i<n*n*n;i++)ddf[i]=-ddf[i];

  return;
 }

/*!\fn void IMFEvaluateThirdDerivativeOfFlow(IMFFlow F, MFNVector vu, MFKVector vp,double *ddf,MFErrorHandler e)
 * \brief Evaluates the third derivative of the flow at a point in phase space (dFdududu).
 *
 * \param F The flow.
 * \param vu A point in phase space.
 * \param vp A point in parameter space.
 * \param dddf An array of length at least the cube of the dimension of the phase space to hold the derivatives. It is
 *           filled with entries in column order. That is dF^i_jkl = df[i+nu*(j+nu*(k+nu*l))], where nu is the dimension
 *           of the phase space.
 */
void IMFEvaluateThirdDerivativeOfFlow(IMFFlow F, MFNVector vu, MFKVector vp,double *dddf, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFEvaluateThirdDerivativeOfFlow"};
  double *u;
  double *p;
  int i,n;

  if(vp!=NULL)p=MFKV_CStar(vp,e);
   else p=NULL;

  if(!strcmp(MFNVGetId(vu,e),"DENSE"))
   {
    u=MFNV_CStar(vu,e);
   }else if(!strcmp(MFNVGetId(vu,e),"IMFExpansionVector"))
   {
    u=IMFExpansionU(IMFExpansionNVGetE(vu,e),e);
   }else{
    sprintf(IMFFlowErrorMsg,"You cannot evaluate a flow at a point which is neither DENSE or and Expansion.");
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    return;
   }

  if(!F->fat)
   {
    F->dddF(u,p,dddf,F->data,e);
   }else{
    sprintf(IMFFlowErrorMsg,"You cannot do this for a fat flow (not implemented).");
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    return;
   }

  n=F->nu;
  if(F->dir==-1)for(i=0;i<n*n*n*n;i++)dddf[i]=-dddf[i];

  return;
 }

/*!\fn int IMFFlowNU(IMFFlow F, MFErrorHandler e)
 * \brief Returns the dimension of the phase space of the flow F.
 *
 * \param F The flow.
 * \returns The dimension of the phase space.
 */
int IMFFlowNU(IMFFlow F, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFFlowNU"};
  return F->nu;
 }

/*!\fn int IMFFlowNP(IMFFlow F, MFErrorHandler e)
 * \brief Returns the dimension of the parameter space of the flow F.
 *
 * \param F The flow.
 * \returns The dimension of the parameter space.
 */
int IMFFlowNP(IMFFlow F, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFFlowNP"};
  return F->np;
 }

/*!\fn void *IMFFlowData(IMFFlow F, MFErrorHandler e)
 * \brief Returns the data block that was provided when F was constructed. BE VERY CAREFUL WITH THIS.
 *
 * \param F The flow.
 * \returns The data block which belongs to the flow.
 */
void *IMFFlowData(IMFFlow F, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFFlowData"};

  return F->data;
 }

/*!\fn MFFlowFreeData IMFFlowFreeData(IMFFlow F, MFErrorHandler e)
 * \brief Returns the routine to free the data block that was provided when F was constructed.
 *
 * \param F The flow.
 * \returns The routine to free the data block.
 */
MFFlowFreeData IMFFlowFreeData(IMFFlow F, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFFlowFreeData"};

  return F->freeData;
 }

void IMFIntegratorDTimeDerivatives(IMFFlow F, double *y, MFKVector p, double *dY, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFIntegratorDTimeDerivatives"};
  double ysave;
  int i,j;
  int n;
  double eps=1.e-3;
  static double *ydotp=NULL;
  static double *ydotm=NULL;
  MFNVector vy;

  n=IMFFlowNU(F,e);

  vy=MFCreateWrappedNVector(n,y,e);
  ydotp=(double*)realloc((void*)ydotp,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(ydotp==NULL)
   {
    sprintf(IMFFlowErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  ydotm=(double*)realloc((void*)ydotm,n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(ydotm==NULL)
   {
    sprintf(IMFFlowErrorMsg,"Out of memory, trying to allocate %d bytes",n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return;
   }
#endif

  for(i=0;i<n;i++)
   {
    ysave=y[i];
    y[i]=ysave+eps;
    IMFEvaluateFlow(F,vy,p,ydotp,e);
    y[i]=ysave-eps;
    IMFEvaluateFlow(F,vy,p,ydotm,e);
    for(j=0;j<n;j++)dY[j+n*i]=.5*(ydotp[j]-ydotm[j])/eps;
    y[i]=ysave;
   }
  MFFreeNVector(vy,e);
  return;
 }

/*!\fn IMFFlow IMFCreateBackwardFlow(IMFFlow F, MFErrorHandler e)
 * \brief Creates the reverse time flow of F.
 *
 * \param F A flow.
 * \returns The reverse time flow of F.
 */
IMFFlow IMFCreateBackwardFlow(IMFFlow F, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFCreateBackwardFlow"};
  IMFFlow result;

  result=IMFCreateFlow(F->nu,F->np,F->F,F->dF,F->dp,F->ddF,F->dddF,F->data,NULL,e);
  result->dir=-1;

  return result;
 }

/*!\fn double IMFFlowR(IMFFlow F, double eps, MFNVector u, MFKVector p, MFNKMatrix mPhi, MFerrorHandler e)
 * \brief Estimates the radius of a negihborhood of a point in phase space within which the flow F
 *        does not change more than eps.
 *
 * \param F A flow.
 * \param eps The bound on the change of F allowed within the neighborhood.
 * \param u A point in phase space.
 * \param p A point in parameter space.
 * \param mPhi An orthonormal basis for a subspace in phase space. The bound on the change in F need only hold in this
 *             subspace (a spherical ball in the subspace centered at the point in phase space).
 * \param e    An error handler.
 * \returns The size of a neighborhood in which the bound applies.
 */
double IMFFlowR(IMFFlow F, double eps, MFNVector u, MFKVector p, MFNKMatrix mPhi, MFErrorHandler e)
 {
  static char RoutineName[]={"IMFFlowR"};
  double *phi;
  double fNorm;
  double R;
  double t;
  int i,j,l,m;
  int n,k;
  static double *dF=NULL;

/* Want to choose R so that |F(u)-F(u+R*Phi s)|<eps  */
/*                                  | R*Fu Phi s |<eps  */

  n=F->nu;
  k=MFNKMatrixK(mPhi,e);

  phi=MFNKM_CStar(mPhi,e);
  dF=(double*)realloc((void*)dF,n*n*sizeof(double));

#ifndef MFNOSAFETYNET
  if(dF==NULL)
   {
    sprintf(IMFFlowErrorMsg,"Out of memory, trying to allocate %d bytes",n*n*sizeof(double));
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    MFErrorHandlerOutOfMemory(e);
    return -1.;
   }
#endif

  IMFEvaluateFlow(F,u,p,dF,e);
  fNorm=0.;for(i=0;i<n;i++)fNorm+=dF[i]*dF[i];fNorm=sqrt(fNorm);
  IMFEvaluateDerivativeOfFlow(F,u,p,dF,e);
  R=0;
  for(i=0;i<k;i++)
   {
    for(j=0;j<k;j++)
     {
      t=0.;
      for(l=0;l<n;l++)
        for(m=0;m<n;m++)t+=phi[l+n*i]*dF[l+n*m]*phi[m+n*j];
      t=t/fNorm;
      if(i==0&&j==0||t>R)R=t;
     }
   }
  R=eps/fabs(R);
/*if(R<.3)R=.3;*/
  
  return R;
 }

void MFFlowSetNVectorFactory(IMFFlow this,MFNVector (*factory)(IMFFlow,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFFlowSetNVectorFactory"};

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(IMFFlowErrorMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  this->nVectorFactory=factory;
  return;
 }

void MFFlowSetKVectorFactory(IMFFlow this,MFKVector (*factory)(IMFFlow,MFErrorHandler), MFErrorHandler e)
 {
  static char RoutineName[]={"MFFlowSetKVectorFactory"};

#ifdef MFNOCONFIDENCE
  if(this==NULL)
   {
    sprintf(IMFFlowErrorMsg,"Implicitly Defined Manifold (argument 1) is NULL");
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

  this->kVectorFactory=factory;
  return;
 }

MFNVector MFFlowNVectorFactory(IMFFlow F, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFlowNVectorFactory"};
  MFNVector Phi;
  int i;
  int rc;
  int verbose=0;

#ifdef MFNOCONFIDENCE
  if(F==NULL)
   {
    sprintf(IMFFlowErrorMsg,"Flow (argument 1) is NULL");
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  if(F->nVectorFactory==NULL)return NULL;
   else return F->nVectorFactory(F,e);
 }

MFKVector MFFlowKVectorFactory(IMFFlow F, MFErrorHandler e)
 {
  static char RoutineName[]={"MFFlowKVectorFactory"};
  MFNVector Phi;
  int i;
  int rc;
  int verbose=0;

#ifdef MFNOCONFIDENCE
  if(F==NULL)
   {
    sprintf(IMFFlowErrorMsg,"Flow (argument 1) is NULL");
    MFSetError(e,12,RoutineName,IMFFlowErrorMsg,__LINE__,__FILE__);
    return;
   }
#endif

#ifdef MFALLOWVERBOSE
  if(verbose){printf("%s\n",RoutineName);fflush(stdout);}
#endif

  if(F->kVectorFactory==NULL)return NULL;
   else return F->kVectorFactory(F,e);
 }

#ifdef __cplusplus
}
#endif
