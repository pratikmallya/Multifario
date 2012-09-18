/* 
    @(#)MFAtlasFriends.h	1.5
    02/04/19 14:45:12
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   August 30, 1999                       */
#ifndef __MFATLASFRIENDS_H__
#define __MFATLASFRIENDS_H__
#include <MFAtlas.h>
#include <MFChart.h>
#include <MFNSpace.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

/* Friend access routines */

MFChart MFAtlasChart(MFAtlas,int,MFErrorHandler);
int MFAtlasNumberOfHalfSpaces(MFAtlas,MFErrorHandler);
int MFAtlasHalfSpaceIndex(MFAtlas,int,MFErrorHandler);
int MFAtlasHalfSpaceLeftChart(MFAtlas,int,MFErrorHandler);
int MFAtlasHalfSpaceRightChart(MFAtlas,int,MFErrorHandler);
MFKVector MFAtlasHalfSpaceNormal(MFAtlas,int,MFErrorHandler);
double MFAtlasHalfSpaceOrigin(MFAtlas,int,MFErrorHandler);
int MFAtlasIsHalfSpaceHyperCube(MFAtlas,int,MFErrorHandler);
int MFAtlasIsHalfSpaceInBothPolytopes(MFAtlas,int,MFErrorHandler);
int MFAtlasLeftPolytope(MFAtlas,int,MFErrorHandler);
int MFAtlasRightPolytope(MFAtlas,int,MFErrorHandler);
MFNSpace MFAtlasNSpace(MFAtlas,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
