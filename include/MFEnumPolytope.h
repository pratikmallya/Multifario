/* 
    @(#)MFEnumPolytope.h	1.2
    02/04/19 14:41:09
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#ifndef __MFENUMPOLYTOPE_H__
#define __MFENUMPOLYTOPE_H__
#include <MFBase.h>

#include <stdio.h>
#include <stdarg.h>
#include <MFPolytope.h>
#include <MFAtlas.h>
#include <MFKVector.h>
#include <MFErrorHandler.h>

struct MFEnumPolytopeSt;
typedef struct MFEnumPolytopeSt *MFEnumPolytope;

#ifdef __cplusplus
 extern "C" {
#endif

MFEnumPolytope MFEnumeratePolytope(MFPolytope,MFErrorHandler);
MFEnumPolytope MFEnumerateAtlas(MFAtlas,MFErrorHandler);
void MFEnumerate(MFEnumPolytope,MFErrorHandler);
void MFFreeEnumPolytope(MFEnumPolytope,MFErrorHandler);

long MFEnumPolytopeNumberOfVertices(MFEnumPolytope,MFErrorHandler);
MFKVector MFEnumPolytopeVertex(MFEnumPolytope,long,MFErrorHandler);
int MFEnumPolytopeNumberOfCells(MFEnumPolytope,int,MFErrorHandler);
int MFEnumPolytopeNumberOfCellFaces(MFEnumPolytope,int,long,MFErrorHandler);
long MFEnumPolytopeCellFace(MFEnumPolytope,int,long,int,MFErrorHandler);
int MFEnumPolytopeNumberOfFaceCells(MFEnumPolytope,int,long,MFErrorHandler);
long MFEnumPolytopeFaceCell(MFEnumPolytope,int,long,int,MFErrorHandler);

int MFEnumPolytopeDimension(MFEnumPolytope,MFErrorHandler);

int MFEnumPolytopeNumberOfCellIndices(MFEnumPolytope,int,long,MFErrorHandler);
int MFEnumPolytopeCellIndex(MFEnumPolytope,int,long,int,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
