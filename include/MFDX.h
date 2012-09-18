/* 
    @(#)MFDX.h	1.3
    02/04/19 14:40:33
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#ifndef __MFDX_H__
#define __MFDX_H__
#include <MFErrorHandler.h>

struct MFEnumDualPolytopeSt;

#ifdef __cplusplus
 extern "C" {
#endif

void MFDualPolytopeToDXFile(char *name,struct MFEnumDualPolytopeSt*,MFErrorHandler);
void MFAtlasToDX(MFAtlas,char *name,MFErrorHandler);
void MFAtlasToDX2(MFAtlas,char *name,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
