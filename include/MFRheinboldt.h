/* 
    @(#)MFKVector.h	1.3
    02/04/19 14:41:38
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   Februaury 20, 1997                     */

#ifndef __MFRHEINBOLDT_H__
#define __MFRHEINBOLDT_H__
#include <MFBase.h>

#include <stdio.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

char *MFRheinboldtGetFilename(MFContinuationMethod,MFErrorHandler);
void MFRheinboldtSetFilename(MFContinuationMethod,char*,MFErrorHandler);
int MFRheinboldtSetIntegerParameter(MFContinuationMethod,char*,int,MFErrorHandler);
int MFRheinboldtSetRealParameter(MFContinuationMethod,char*,double,MFErrorHandler);
int MFRheinboldtGetIntegerParameter(MFContinuationMethod,char*,MFErrorHandler);
double MFRheinboldtGetRealParameter(MFContinuationMethod,char*,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
