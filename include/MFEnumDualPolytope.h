/* 
    @(#)MFEnumDualPolytope.h	1.2
    02/04/19 14:41:03
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#ifndef __MFENUMDUALPOLYTOPE_H__
#define __MFENUMDUALPOLYTOPE_H__
#include <MFBase.h>

#include <stdio.h>
#include <stdarg.h>
#include <MFEnumPolytope.h>
#include <MFPolytope.h>
#include <MFAtlas.h>
#include <MFNVector.h>
#include <MFErrorHandler.h>

struct MFEnumDualPolytopeSt;
typedef struct MFEnumDualPolytopeSt *MFEnumDualPolytope;

#ifdef __cplusplus
 extern "C" {
#endif

struct MFEnumPolytopeSt;

MFEnumDualPolytope MFCreateDualOfEnumPolytope(struct MFEnumPolytopeSt*,MFErrorHandler);
MFEnumDualPolytope MFEnumDualOfAtlas(MFAtlas,MFErrorHandler);

void MFEnumerateDual(MFEnumDualPolytope,MFErrorHandler);
void MFPrintEnumDualPolytope(FILE*,MFEnumDualPolytope,MFErrorHandler);
void MFFreeEnumDualPolytope(MFEnumDualPolytope,MFErrorHandler);

long MFEnumDualPolytopeNumberOfVertices(MFEnumDualPolytope,MFErrorHandler);
MFNVector MFEnumDualPolytopeVertex(MFEnumDualPolytope,long,MFErrorHandler);
int MFEnumDualPolytopeNumberOfCells(MFEnumDualPolytope,int,MFErrorHandler);
int MFEnumDualPolytopeNumberOfCellFaces(MFEnumDualPolytope,int,long,MFErrorHandler);
long MFEnumDualPolytopeCellFace(MFEnumDualPolytope,int,long,int,MFErrorHandler);
int MFEnumDualPolytopeNumberOfFaceCells(MFEnumDualPolytope,int,long,MFErrorHandler);
long MFEnumDualPolytopeFaceCell(MFEnumDualPolytope,int,long,int,MFErrorHandler);

int MFEnumDualPolytopeDimension(MFEnumDualPolytope,MFErrorHandler);

int MFEnumDualPolytopeNumberOfCellIndices(MFEnumDualPolytope,int,long,MFErrorHandler);
int MFEnumDualPolytopeCellIndex(MFEnumDualPolytope,int,long,int,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
