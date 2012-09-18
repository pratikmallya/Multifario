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

static char *id="@(#) $Id: ComputeSpherePacking.c,v 1.3 2011/07/21 17:43:45 mhender Exp $";

#include <MFAtlas.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <MFMultifariosMethod.h>

int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  int i,n;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0;
  FILE *fid;
  char *vars;
  char *expr;
  char *digits[]={"0","1","2","3","4","5","6","7","8","9"};
  char name[80]="";
  MFContinuationMethod H;
  MFErrorHandler e;

  e=MFCreateErrorHandler();

  n=4;
  if(argc>1)sscanf(argv[1],"%d",&n);
  if(n<1||n>9)
   {
    fprintf(stderr,"%s currently only works for 0<n<10. You supplied n=%d\n",n);
    fflush(stderr);
    return 12;
   }
  vars=(char*)malloc((2+3*n)*sizeof(char));
  strcpy(vars,"[");
  for(i=0;i<n;i++)
   {
    if(i>0)strcat(vars,",");
    strcat(vars,"v");
    strcat(vars,digits[i]);
   }
  strcat(vars,"]");

  expr=(char*)malloc((5+6*n)*sizeof(char));
  strcpy(expr,"[");
  for(i=0;i<n;i++)
   {
    if(i>0)strcat(expr,"+");
    strcat(expr,"v");
    strcat(expr,digits[i]);
    strcat(expr,"**2");
   }
  strcat(expr,"-1.]");

  printf("Variables : %s\n",vars);
  fflush(stdout);
  printf("Expression: %s\n",expr);
  fflush(stdout);

  M=MFIMFCreateAlgebraicExpression(vars,expr,e);
  MFIMFSetR(M,sqrt(3.)/2.1,e);
  Omega=MFNRegionCreateHyperCube(n,3.,e);

  u0=MFIMFVectorFactory(M,e);
  MFNVSetC(u0,0, 1.,e);
  for(i=1;i<n;i++)MFNVSetC(u0,i, 0.,e);

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",1.,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",0,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  strcpy(name,"SpherePacking");
  strcat(name,digits[n]);
  MFMultifarioSetFilename(H,name,e);

  S=MFComputeAtlas(H,M,Omega,u0,e);

  MFCloseAtlas(H,S,e);
  printf("Done computing manifold\n");fflush(stdout);

  MFFreeAtlas(S,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  MFFreeNVector(u0,e);
  MFFreeContinuationMethod(H,e);
  free(vars);
  free(expr);

  MFFreeErrorHandler(e);

  return 0;
 }
