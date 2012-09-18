/* 
    @(#)MFDrawChartState.h	1.4
    02/04/19 14:40:50
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#ifndef __MFDRAWCHARTSTATE_H__
#define __MFDRAWCHARTSTATE_H__
#include <MFBase.h>

#include <MFChart.h>
#include <MFErrorHandler.h>

struct MFChartStateSt;
typedef struct MFChartStateSt *MFChartState;

#ifdef __cplusplus
 extern "C" {
#endif

MFChartState MFCreateChartState(MFChart,MFErrorHandler);
void MFDrawChartFromState(MFChartState,MFErrorHandler);
MFChartState MFAddChartState(MFChart,MFErrorHandler);
int MFNChartStates(MFErrorHandler);
void MFClearChartState(MFChartState,MFErrorHandler);
void MFSetChartStateColor(MFChartState,int,int,int,int,MFErrorHandler);
void MFGetChartStateColor(MFChartState,int,int*,int*,int*,MFErrorHandler);
void MFChartStateAddPolygon(MFChartState,int,float*,float*,float*,int,MFErrorHandler);
void MFChartStateAddPolygonWithNormal(MFChartState,int,float*,float*,float*,float*,int,MFErrorHandler);
void MFChartStateAddLine(MFChartState,float,float,float,float,float,float,MFErrorHandler);
MFChartState MFGetChartState(int,MFErrorHandler);
int MFChartStateChanged(MFChartState,MFErrorHandler);
void MFFreeDrawCharts(MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
