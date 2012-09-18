/* 
    @(#)MFPolyhedralComplex.h	1.2
    02/04/19 14:41:09
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#ifndef __MFPOLYHEDRALCOMPLEX_H__
#define __MFPOLYHEDRALCOMPLEX_H__
#include <MFBase.h>

#include <stdio.h>
#include <stdarg.h>
#include <MFPolytope.h>
#include <MFAtlas.h>
#include <MFKVector.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

struct MFCellSt;
typedef struct MFCellSt *MFCell;
struct MFPolyhedralCellComplexSt;
typedef struct MFPolyhedralCellComplexSt *MFPolyhedralCellComplex;

#ifdef __cplusplus
}
#endif

#endif
