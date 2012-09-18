/* 
    @(#)MFPrint.h	1.4
    02/07/24 09:33:25
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   February 24, 1999                     */

#ifndef __MFPRINT_H__
#define __MFPRINT_H__
#include <MFBase.h>

#include <stdio.h>

#include <MFKVector.h>
#include <MFNVector.h>
#include <MFNKMatrix.h>
#include <MFChart.h>
#include <MFPolytope.h>
#include <MFAtlas.h>
#include <MFBinaryTree.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

void MFPrintKVector(FILE*,MFKVector,MFErrorHandler);
void MFPrintNVector(FILE*,MFNVector,MFErrorHandler);
void MFPrintNKMatrix(FILE*,MFNKMatrix,MFErrorHandler);
void MFPrintChart(FILE*,MFChart,MFErrorHandler);
void MFPrintPolytope(FILE*,MFPolytope,MFErrorHandler);
void MFPrintPolytopeTerse(FILE*,MFPolytope,MFErrorHandler);
void MFPrintAtlas(FILE*,MFAtlas,MFErrorHandler);
void MFPrintAtlasFaceList(FILE*,MFAtlas,MFErrorHandler);
void MFPrintBinaryTree(FILE*,MFBinaryTree,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
