#ifndef __MFLOCA_H__
#define __MFLOCA_H__
/*
    @(#)MFLOCA.h	1.2
    02/07/18 15:08:33

    PROGRAM NAME:  Manifold

    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.

    Please refer to the LICENSE file in the top directory

    author: Mike Henderson
    date:   April 21, 2002 for LOCA

*/

#include <MFNSpace.h>
#include <MFNRegion.h>
#include <MFImplicitMF.h>
#include <MFNVector.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

MFNSpace MFCreateLOCANSpace(void*,MFErrorHandler);
MFNRegion MFNRegionCreateLOCA(void*,MFErrorHandler);
MFImplicitMF MFIMFCreateLOCA(int,void*,MFErrorHandler);
MFNVector MFCreateLOCANVectorWithData(int,double*,int,double*,MFErrorHandler);
MFNVector MFCreateLOCANVector(int,int,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
