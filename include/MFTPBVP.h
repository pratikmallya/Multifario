/* 
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   February 12, 2007
 */

#include <MFErrorHandler.h>
#include <MFNVector.h>
#include <stdio.h>

#ifdef __cplusplus
 extern "C" {
#endif

int MFAtlasDetSign(int,double*,MFErrorHandler);
void MFTPBVPTestTangentRaw(int,int,double*,double*,void*,MFErrorHandler);
int MFSolveBordered(int,int,int,int,double*,double*,double*,double*,double*,double*,double*,double*,MFErrorHandler);
void MFTestSolveBordered(int,int,int,int,int,int,double*,double*,double*,double*,double*,double*,double*,double*,MFErrorHandler);
int MFAnalyzeBordered(int,int,int,int,int,double*,double*,double*,double*,MFErrorHandler);
int MFTPBVPAnalyze(MFImplicitMF,MFNVector,MFErrorHandler);
int MFSolveFull(int,double*,double*,MFErrorHandler);
void MFSetBandedMatrixElement(int,int,double*,double,int,int,int,MFErrorHandler);
void MFIncrementBandedMatrixElement(int,int,double*,double,int,int,int,MFErrorHandler);
void MFPrintBorderedBandedMatrix(FILE*,int,int,int,int,double*,double*,double*,double*,MFErrorHandler);
void MFPrintBorderedBandedMatrixByBlock(FILE*,int,int,int,int,int,int,double*,double*,double*,double*,MFErrorHandler);
void MFPrintBorderedBandedMatrixMinusFull(FILE*,int,int,int,int,double*,double*,double*,double*,double*,MFErrorHandler);
void MFTPBVPGetJac(int,int,void*,double*,double*,double*,double**,double**,double**,double**,MFErrorHandler);
double MFTPBVPGetRes(int,int,void*,double*,double*,double*,double**,double**,MFErrorHandler);
void MFTPBVPGramSchmidtNoMat(MFNSpace,int,int,int,int,int,double*,MFErrorHandler);
int MFTPBVPTestJacobian(int,int,int,int,int,int,int,double*,double*,void*,MFErrorHandler);

int MFFullNZeroSV(int,int,double*,MFErrorHandler);
int MFFullNPosEV(int,int,double*,MFErrorHandler);

static void MFFreeTPBVPData(void*,MFErrorHandler);
static int MFProjectTPBVP(int,int,MFNVector,MFNKMatrix,MFNVector,void*,int*,MFErrorHandler);
static int MFTangentTPBVP(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static int MFTangentTPBVPWithGuess(int,int,MFNVector,MFNKMatrix,MFNKMatrix,void*,MFErrorHandler);
static int MFTangentTPBVPSingle(int,int,int,double*,double*,double*,double*,void*,MFErrorHandler);
static void MFCurvatureTPBVP(int,int,double*,double*,double*,double*,double*,double*,void*,MFErrorHandler);
static double MFTPBVPCurvatureSingle(int,int,int,double*,double*,double*,void*,MFErrorHandler);
static double MFScaleTPBVP(int,int,MFNVector,MFNKMatrix,void*,MFErrorHandler);
static void MFWriteTPBVPData(FILE*,void*,MFErrorHandler);
static MFImplicitMF MFReadTPBVP(FILE*,MFErrorHandler);
static int MFSingularTPBVP(int,int,MFNVector,MFNKMatrix,MFNVector,void*,MFErrorHandler);
static int MFTPBVPProjectToSave(MFNVector,double*,void*,MFErrorHandler);
static int MFTPBVPProjectToDraw(MFNVector,double*,void*,MFErrorHandler);
static int MFTPBVPProjectForBB(MFNVector,double*,void*,MFErrorHandler);

int MFTPBVPTransformX(int,int,int,double*,double**,MFErrorHandler);

double MFEVProd(int,int,double*,MFErrorHandler);
void MFSVD(int,double*,double*,double*,double*,MFErrorHandler);

MFNVector MFNVectorFactory(MFImplicitMF,MFErrorHandler);
MFNKMatrix MFNKMatrixFactory(MFImplicitMF,MFErrorHandler);

int MFTPBVPGetNX(MFImplicitMF, MFErrorHandler);
int MFTPBVPGetNU(MFImplicitMF, MFErrorHandler);
int MFTPBVPGetNP(MFImplicitMF, MFErrorHandler);
int MFTPBVPGetNIC(MFImplicitMF, MFErrorHandler);
int MFTPBVPGetNBC(MFImplicitMF, MFErrorHandler);

void MFTPBVPEvaluateIntegralConstraints(MFImplicitMF,MFNVector,MFNVector,double*,MFErrorHandler);
void MFTPBVPEvaluateBoundaryConditions(MFImplicitMF,MFNVector,MFNVector,double*,MFErrorHandler);
void MFTPBVPNSpaceSetPeriodicParameter(MFNSpace,int,double,MFErrorHandler);

void MFTPBVPSetStability(MFImplicitMF,MFNVector,MFNKMatrix,void*,MFErrorHandler);
int MFStopTPBVP(MFImplicitMF,MFNVector,MFNKMatrix,MFNVector,MFNKMatrix,void*,MFErrorHandler);

#ifdef __cplusplus
}
#endif
