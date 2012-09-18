/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   December 9, 2002   modified InvMF
 */

static char *id="@(#) $Id: IMFThreeDFlows.c,v 1.4 2010/05/21 12:22:04 mhender Exp $";

#include <math.h>
#include <stdio.h>
#include <IMFFlow.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

static void source(double *u, double *p, double *F,  void *data, MFErrorHandler e)
 {
  F[0]=u[0];
  F[1]=u[1];
  F[2]=u[2];

  return;
 }

static void dsource(double *u, double *p, double *dF,  void *data, MFErrorHandler e)
 {
  dF[0+3*0]=1.;
  dF[0+3*1]=0.;
  dF[0+3*2]=0.;

  dF[1+3*0]=0.;
  dF[1+3*1]=1.;
  dF[1+3*2]=0.;

  dF[2+3*0]=0.;
  dF[2+3*1]=0.;
  dF[2+3*2]=1.;

  return;
 }

static void ddsource(double *u, double *p, double *ddF,  void *data, MFErrorHandler e)
 {
  int i;

  for(i=0;i<3*3*3;i++)ddF[i]=0.;

  return;
 }

static void dddsource(double *u, double *p, double *dddF,  void *data, MFErrorHandler e)
 {
  int i;

  for(i=0;i<3*3*3*3;i++)dddF[i]=0.;
  return;
 }
/*! \fn IMFFlow IMFCreateThreeDSourceFlow(MFErrorHandler e);
 *  \brief A flow in three dimensional phase space with a source at the origin
 *
 *  \param e A place to return errors.
 *  \returns A new flow.
 */
IMFFlow IMFCreateThreeDSourceFlow(MFErrorHandler e)
 {
  return IMFCreateFlow(3,0,source,dsource,NULL,ddsource,dddsource,NULL,NULL,e);
 }

static void hyperbolic(double *u, double *p, double *F,  void *data, MFErrorHandler e)
 {
  F[0]=u[0];
  F[1]=u[1];
  F[2]=-u[2];

  return;
 }

static void dhyperbolic(double *u, double *p, double *dF,  void *data, MFErrorHandler e)
 {
  dF[0+3*0]=1.;
  dF[0+3*1]=0.;
  dF[0+3*2]=0.;

  dF[1+3*0]=0.;
  dF[1+3*1]=1.;
  dF[1+3*2]=0.;

  dF[2+3*0]=0.;
  dF[2+3*1]=0.;
  dF[2+3*2]=-1.;

  return;
 }

static void ddhyperbolic(double *u, double *p, double *ddF,  void *data, MFErrorHandler e)
 {
  int i;

  for(i=0;i<3*3*3;i++)ddF[i]=0.;

  return;
 }

static void dddhyperbolic(double *u, double *p, double *dddF,  void *data, MFErrorHandler e)
 {
  int i;

  for(i=0;i<3*3*3*3;i++)dddF[i]=0.;
  return;
 }

/*! \fn IMFFlow IMFCreateThreeDHyperbolicFlow(MFErrorHandler e);
 *  \brief A flow in three dimensional phase space with a hyperbolic fixed point at the origin with the xy-plane
 *         the unstable invariant manifold and the z-axis the stable manifold.
 *
 *  \param e A place to return errors.
 *  \returns A new flow.
 */
IMFFlow IMFCreateThreeDHyperbolicFlow(MFErrorHandler e)
 {
  return IMFCreateFlow(3,0,hyperbolic,NULL,dhyperbolic,ddhyperbolic,dddhyperbolic,NULL,NULL,e);
 }

static void saddlefocus(double *u, double *p, double *F,  void *data, MFErrorHandler e)
 {
  double eps;

  eps=p[0];

  F[0]=u[0];
  F[1]=u[1]+eps*u[2];
  F[2]=u[2]-eps*u[1];

  return;
 }

static void dpsaddlefocus(double *u, double *p, double *dF,  void *data, MFErrorHandler e)
 {
  dF[0]=0.;
  dF[1]=u[2];
  dF[2]=-u[1];

  return;
 }

static void dsaddlefocus(double *u, double *p, double *dF,  void *data, MFErrorHandler e)
 {
  double eps;

  eps=p[0];

  dF[0+3*0]=1.;
  dF[0+3*1]=0.;
  dF[0+3*2]=0.;

  dF[1+3*0]=0.;
  dF[1+3*1]=1.;
  dF[1+3*2]=eps;

  dF[2+3*0]=0.;
  dF[2+3*1]=-eps;
  dF[2+3*2]=1.;

  return;
 }

static void ddsaddlefocus(double *u, double *p, double *ddF,  void *data, MFErrorHandler e)
 {
  int i;

  for(i=0;i<3*3*3;i++)ddF[i]=0.;

  return;
 }

static void dddsaddlefocus(double *u, double *p, double *dddF,  void *data, MFErrorHandler e)
 {
  int i;

  for(i=0;i<3*3*3*3;i++)dddF[i]=0.;
  return;
 }

/*! \fn IMFFlow IMFCreateThreeDSaddleFocusFlow(MFErrorHandler e);
 *  \brief A flow in three dimensional phase space with a hyperbolic fixed point at the origin with the xy-plane
 *         the unstable invariant manifold with complex conjugate eigenvalues, and the z-axis the stable manifold.
 *         The stable eigenvalues are determined by the parameter eps, and are 1+i*eps and 1-i*eps.
 *
 *  \param e A place to return errors.
 *  \returns A new flow.
 */
IMFFlow IMFCreateThreeDSaddleFocusFlow(MFErrorHandler e)
 {
  return IMFCreateFlow(3,1,saddlefocus,dpsaddlefocus,dsaddlefocus,ddsaddlefocus,dddsaddlefocus,NULL,NULL,e);
 }

static void saddlecenter(double *u, double *p, double *F,  void *data, MFErrorHandler e)
 {
  F[0]=-u[0];
  F[1]=u[2];
  F[2]=-u[1];
  printf("saddleCenter(%lf,%lf,%lf)=(%lf,%lf,%lf)\n",u[0],u[1],u[2],F[0],F[1],F[2]);fflush(stdout);

  return;
 }

static void dsaddlecenter(double *u, double *p, double *dF,  void *data, MFErrorHandler e)
 {
  dF[0+3*0]=-1.;
  dF[0+3*1]=0.;
  dF[0+3*2]=0.;

  dF[1+3*0]=0.;
  dF[1+3*1]=0.;
  dF[1+3*2]=1.;

  dF[2+3*0]=0.;
  dF[2+3*1]=-1.;
  dF[2+3*2]=0.;

  return;
 }

static void ddsaddlecenter(double *u, double *p, double *ddF,  void *data, MFErrorHandler e)
 {
  int i;

  for(i=0;i<3*3*3;i++)ddF[i]=0.;

  return;
 }

static void dddsaddlecenter(double *u, double *p, double *dddF,  void *data, MFErrorHandler e)
 {
  int i;

  for(i=0;i<3*3*3*3;i++)dddF[i]=0.;
  return;
 }

/*! \fn IMFFlow IMFCreateThreeDSaddleCenterFlow(MFErrorHandler e);
 *  \brief A flow in three dimensional phase space with a fixed point at the origin with the xy-plane
 *         its center invariant manifold and the z-axis the stable manifold.
 *
 *  \param e A place to return errors.
 *  \returns A new flow.
 */
IMFFlow IMFCreateThreeDSaddleCenterFlow(MFErrorHandler e)
 {
  return IMFCreateFlow(3,0,saddlecenter,NULL,dsaddlecenter,ddsaddlecenter,dddsaddlecenter,NULL,NULL,e);
 }


static void lorenz(double *u, double *p, double *F,  void *data, MFErrorHandler e)
 {
  double sigma,rho,beta;

  sigma=p[0];
  rho=p[1];
  beta=p[2];
  F[0]=sigma*(u[1]-u[0]);
  F[1]=rho*u[0]-u[1]-u[0]*u[2];
  F[2]=u[0]*u[1]-beta*u[2];

  return;
 }

static void dplorenz(double *u, double *p, double *dF,  void *data, MFErrorHandler e)
 {
  double sigma,rho,beta;

  sigma=p[0];
  rho=p[1];
  beta=p[2];

  dF[0+3*0]=u[1]-u[0];
  dF[1+3*0]=0.;
  dF[2+3*0]=0.;

  dF[0+3*1]=0.;
  dF[1+3*1]=u[0];
  dF[2+3*1]=0.;

  dF[0+3*2]=0.;
  dF[1+3*2]=0.;
  dF[2+3*2]=-u[2];

  return;
 }

static void dlorenz(double *u, double *p, double *dF,  void *data, MFErrorHandler e)
 {
  double sigma,rho,beta;

  sigma=p[0];
  rho=p[1];
  beta=p[2];

  dF[0+3*0]=-sigma;
  dF[0+3*1]=sigma;
  dF[0+3*2]=0.;

  dF[1+3*0]=rho-u[2];
  dF[1+3*1]=-1.;
  dF[1+3*2]=-u[0];

  dF[2+3*0]=u[1];
  dF[2+3*1]=u[0];
  dF[2+3*2]=-beta;

  return;
 }

static void ddlorenz(double *u, double *p, double *ddF,  void *data, MFErrorHandler e)
 {
  ddF[0+3*(0+3*0)]=0.;
  ddF[0+3*(1+3*0)]=0.;
  ddF[0+3*(2+3*0)]=0.;

  ddF[1+3*(0+3*0)]=0.;
  ddF[1+3*(1+3*0)]=0.;
  ddF[1+3*(2+3*0)]=-1;

  ddF[2+3*(0+3*0)]=0.;
  ddF[2+3*(1+3*0)]=1.;
  ddF[2+3*(2+3*0)]=0.;

  ddF[0+3*(0+3*1)]=0.;
  ddF[0+3*(1+3*1)]=0.;
  ddF[0+3*(2+3*1)]=0.;

  ddF[1+3*(0+3*1)]=0.;
  ddF[1+3*(1+3*1)]=0.;
  ddF[1+3*(2+3*1)]=0.;

  ddF[2+3*(0+3*1)]=1.;
  ddF[2+3*(1+3*1)]=0.;
  ddF[2+3*(2+3*1)]=0.;

  ddF[0+3*(0+3*2)]=0.;
  ddF[0+3*(1+3*2)]=0.;
  ddF[0+3*(2+3*2)]=0.;

  ddF[1+3*(0+3*2)]=-1.;
  ddF[1+3*(1+3*2)]=0.;
  ddF[1+3*(2+3*2)]=0.;

  ddF[2+3*(0+3*2)]=0.;
  ddF[2+3*(1+3*2)]=0.;
  ddF[2+3*(2+3*2)]=0.;

  return;
 }

static void dddlorenz(double *u, double *p, double *dddF,  void *data, MFErrorHandler e)
 {
  int i;

  for(i=0;i<3*3*3*3;i++)dddF[i]=0.;
  return;
 }

/*! \fn IMFFlow IMFCreateLorenzFlow(MFErrorHandler e)
 *  \brief The Lorenz flow. The three parameers are sigmma, rho and beta.
 *
 *  \param e A place to return errors.
 *  \returns A new flow.
 */
IMFFlow IMFCreateLorenzFlow(MFErrorHandler e)
 {
  return IMFCreateFlow(3,3,lorenz,dlorenz,dplorenz,ddlorenz,dddlorenz,NULL,NULL,e);
 }


static void standardlorenz(double *u, double *p, double *F,  void *data, MFErrorHandler e)
 {
  double sigma,rho,beta;

  sigma=10.;
  rho=28.;
  beta=8./3.;
 
  F[0]=sigma*(u[1]-u[0]);
  F[1]=rho*u[0]-u[1]-u[0]*u[2];
  F[2]=u[0]*u[1]-beta*u[2];

  return;
 }

static void dstandardlorenz(double *u, double *p, double *dF,  void *data, MFErrorHandler e)
 {
  double sigma,rho,beta;

  sigma=10.;
  rho=28.;
  beta=8./3.;

  dF[0+3*0]=-sigma;
  dF[0+3*1]=sigma;
  dF[0+3*2]=0.;

  dF[1+3*0]=rho-u[2];
  dF[1+3*1]=-1.;
  dF[1+3*2]=-u[0];

  dF[2+3*0]=u[1];
  dF[2+3*1]=u[0];
  dF[2+3*2]=-beta;

  return;
 }

static void ddstandardlorenz(double *u, double *p, double *ddF,  void *data, MFErrorHandler e)
 {
  ddF[0+3*(0+3*0)]=0.;
  ddF[0+3*(1+3*0)]=0.;
  ddF[0+3*(2+3*0)]=0.;

  ddF[1+3*(0+3*0)]=0.;
  ddF[1+3*(1+3*0)]=0.;
  ddF[1+3*(2+3*0)]=-1;

  ddF[2+3*(0+3*0)]=0.;
  ddF[2+3*(1+3*0)]=1.;
  ddF[2+3*(2+3*0)]=0.;

  ddF[0+3*(0+3*1)]=0.;
  ddF[0+3*(1+3*1)]=0.;
  ddF[0+3*(2+3*1)]=0.;

  ddF[1+3*(0+3*1)]=0.;
  ddF[1+3*(1+3*1)]=0.;
  ddF[1+3*(2+3*1)]=0.;

  ddF[2+3*(0+3*1)]=1.;
  ddF[2+3*(1+3*1)]=0.;
  ddF[2+3*(2+3*1)]=0.;

  ddF[0+3*(0+3*2)]=0.;
  ddF[0+3*(1+3*2)]=0.;
  ddF[0+3*(2+3*2)]=0.;

  ddF[1+3*(0+3*2)]=-1.;
  ddF[1+3*(1+3*2)]=0.;
  ddF[1+3*(2+3*2)]=0.;

  ddF[2+3*(0+3*2)]=0.;
  ddF[2+3*(1+3*2)]=0.;
  ddF[2+3*(2+3*2)]=0.;

  return;
 }

static void dddstandardlorenz(double *u, double *p, double *dddF,  void *data, MFErrorHandler e)
 {
  int i;

  for(i=0;i<3*3*3*3;i++)dddF[i]=0.;
  return;
 }

/*! \fn IMFFlow IMFFlow IMFCreateStandardLorenzFlow(MFErrorHandler e);
 *  \brief The "Standard" Lorenz flow.  sigma=10.; rho=28.; beta=8./3.;
 *
 *  \param e A place to return errors.
 *  \returns A new flow.
 */
IMFFlow IMFCreateStandardLorenzFlow(MFErrorHandler e)
 {
  return IMFCreateFlow(3,0,standardlorenz,dstandardlorenz,NULL,ddstandardlorenz,dddstandardlorenz,NULL,NULL,e);
 }

#ifdef __cplusplus
}
#endif
