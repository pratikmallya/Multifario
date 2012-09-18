/* 
    @(#)MFDraw.h	1.5
    02/04/19 14:40:43
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#ifndef __MFDRAW_H__
#define __MFDRAW_H__
#include <MFBase.h>
#include <MFAtlas.h>
#include <MFErrorHandler.h>

struct MFChartStateSt;
struct MFChartSt;
struct MFEnumPolytopeSt;
struct MFEnumDualPolytopeSt;

#ifdef __cplusplus
 extern "C" {
#endif

void MFDrawInitialize(float,float,MFErrorHandler);
void MFDrawInitializeCube(float,float,float,float,float,float,float,float,MFErrorHandler);
void MFDrawInitializeNoCube(float,float,float,float,float,float,float,float,MFErrorHandler);
void MFDrawInitializeFromFile(char*,MFErrorHandler);
void MFDrawClose(MFErrorHandler);
void MFDrawClear(MFErrorHandler);
void MFDrawDisplay(MFErrorHandler);
void MFDrawAtlas(MFAtlas,MFErrorHandler);
void MFDrawAtlasOnce(MFAtlas,MFErrorHandler);
void MFDrawAtlasTS(MFAtlas,MFErrorHandler);
void MFDrawEnumPolytope(struct MFEnumPolytopeSt*,struct MFChartSt*,struct MFChartStateSt*,MFErrorHandler);
void MFDrawEnumDualPolytope(struct MFEnumDualPolytopeSt*,MFErrorHandler);
void MFDrawEnumDualPolytopeEdges(struct MFEnumDualPolytopeSt*,MFErrorHandler);
void MFPendulaPeriodic(int,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
