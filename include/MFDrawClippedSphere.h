/* 
    @(#)MFDrawClippedSphere.h	1.2
    02/04/19 14:40:57
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#ifndef __MFDRAWCLIPPEDSPHERE_H__
#define __MFDRAWCLIPPEDSPHERE_H__

#include <MFBase.h>
#include <MFChart.h>
#include <MFDrawChartState.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

void MFDrawClippedSphere(MFChart,double,double,double,MFChartState, MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
