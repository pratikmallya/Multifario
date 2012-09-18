/* 
    %W%
    %D% %T%
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   November 11, 1997                     */
/*              February 2, 1999   Ported to C        */

#ifndef __MFNKMATRIX_H__
#define __MFNKMATRIX_H__
#include <MFBase.h>

#include <MFNVector.h>
#include <MFKVector.h>
#include <MFNSpace.h>
#include <MFErrorHandler.h>

struct MFNKMatrixSt;
typedef struct MFNKMatrixSt *MFNKMatrix;

#ifdef __cplusplus
 extern "C" {
#endif

MFNKMatrix MFCreateNKMatrix(int,MFNVector*,MFErrorHandler);
MFNKMatrix MFCreateNKMatrixWithData(int,int,double*,MFErrorHandler);
void MFRefNKMatrix(MFNKMatrix,MFErrorHandler);
void MFFreeNKMatrix(MFNKMatrix,MFErrorHandler);

int MFNKMatrixK(MFNKMatrix,MFErrorHandler);
int MFNKMatrixN(MFNKMatrix,MFErrorHandler);
MFNVector MFMColumn(MFNKMatrix,int,MFErrorHandler);
void MFMRow(MFNKMatrix,int,MFKVector,MFErrorHandler);
void MFNKMSetC(MFNKMatrix,int,int,double,MFErrorHandler);

void MFMVMul(MFNSpace,MFNKMatrix,MFKVector,MFNVector,MFErrorHandler);
void MFMVMulT(MFNSpace,MFNKMatrix,MFNVector,MFKVector,MFErrorHandler);

void MFWriteNKMatrix(FILE*,MFNKMatrix,MFErrorHandler);
MFNKMatrix MFReadNKMatrix(FILE*,MFErrorHandler);

void MFGramSchmidt(MFNSpace,MFNKMatrix,MFErrorHandler);
void MFGramSchmidtNoMat(int,int,double*,MFErrorHandler);

void MFNKMProjectTangentForBranchSwitch(MFNSpace,MFNKMatrix,MFNVector,MFNKMatrix,MFErrorHandler);
MFNKMatrix MFCloneNKMatrix(MFNKMatrix,MFErrorHandler);
void MFMSetColumn(MFNKMatrix,int,MFNVector,MFErrorHandler);

void MFGramSchmidtPlusOne(MFNSpace,MFNVector,MFNKMatrix,MFNKMatrix,MFErrorHandler);
void MFGramSchmidtReplace(MFNSpace,MFNKMatrix,MFNVector,MFNVector,MFErrorHandler);

void MFMMMul(MFNSpace,MFNKMatrix,MFNKMatrix,double*,MFErrorHandler);

double *MFNKM_CStar(MFNKMatrix,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
