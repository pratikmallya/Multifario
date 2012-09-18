/* 
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   February 20 2007                     */


#ifndef __MFFORCED_H__
#define __MFFORCED_H__
#include <MFBase.h>

#include <MFNVector.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

void ForcedWrite(MFNVector, MFErrorHandler);
void ForcedCloser(MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
