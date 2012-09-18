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

static char *id="@(#) $Id: ComputeLine.c,v 1.3 2011/07/21 17:43:45 mhender Exp $";

#include <MFAtlas.h>
#include <MFNRegion.h>
#include <MFNVector.h>
#include <math.h>
#include <MFMultifariosMethod.h>

int main(int argc, char *argv[])
 {
  MFImplicitMF M;
  int i;
  int n,k;
  MFNRegion Omega;
  MFAtlas S;
  MFNVector u0;
  FILE *fid;
  MFContinuationMethod H;
  MFErrorHandler e;

  e=MFCreateErrorHandler();

  k=1;
  M=MFIMFCreateNSpaceWithRadius(k,.1,e);
  n=MFIMF_N(M,e);
  Omega=MFNRegionCreateHyperCube(n,.5,e);

  u0=MFIMFVectorFactory(M,e);
  for(i=0;i<k;i++)
   {
    MFNVSetC(u0,i, 0.,e);
   }

  H=MFCreateMultifariosMethod(e);
  MFMultifarioSetRealParameter(H,"epsilon",.05,e);
  MFMultifarioSetIntegerParameter(H,"maxCharts",-1,e);
  MFMultifarioSetIntegerParameter(H,"verbose",1,e);
  MFMultifarioSetIntegerParameter(H,"page",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToPlotFile",1,e);
  MFMultifarioSetIntegerParameter(H,"dumpToCenterFile",0,e);
  MFMultifarioSetFilename(H,"Line",e);

  S=MFComputeAtlas(H,M,Omega,u0,e);

/* Use this is Dump to Plotfile is set to zero */
/*
  fid=fopen("Line.atlas","w");
  MFWriteAtlas(fid,S,e);
  fclose(fid);
  printf("Done writing Atlas\n");fflush(stdout);
*/

  MFCloseAtlas(H,S,e);
  printf("Done computing Atlas\n");fflush(stdout);

  MFFreeAtlas(S,e);
  MFFreeImplicitMF(M,e);
  MFFreeNRegion(Omega,e);
  MFFreeContinuationMethod(H,e);
  MFFreeNVector(u0,e);

  MFFreeErrorHandler(e);

  return 0;
 }
