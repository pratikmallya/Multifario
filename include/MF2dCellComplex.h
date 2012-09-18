/*
    %W%
    %D% %T%

    PROGRAM NAME:  Manifold

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   February 22, 1999                     */

#ifndef __MF2DCELLCOMPLEX_H__
#define __MF2DCELLCOMPLEX_H__
#include <stdio.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct MF2dCellComplexSt;
typedef struct MF2dCellComplexSt *MF2dCellComplex;

MF2dCellComplex MFPlotfileDual(FILE *fid, MFErrorHandler e);
void MFOutput2dCellComplexAsEasyMesh(MF2dCellComplex C, char *name,MFErrorHandler);
void MFOutput2dCellComplexAsDX(MF2dCellComplex C, char *name,MFErrorHandler);
void MFFree2dCellComplex(MF2dCellComplex C,MFErrorHandler e);

#ifdef __cplusplus
}
#endif
#endif
