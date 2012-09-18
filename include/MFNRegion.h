/* 
    @(#)MFNRegion.h	1.11
    02/09/27 11:30:13
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   November 11, 1997                     */
/*              February 2, 1999   Ported to C        */

#ifndef __MFNREGION_H__
#define __MFNREGION_H__
#include <MFBase.h>

#include <stdio.h>
#include <MFNVector.h>
#include <MFNSpace.h>
#include <MFErrorHandler.h>

struct MFNRegionSt;
typedef struct MFNRegionSt *MFNRegion;

#ifdef __cplusplus
 extern "C" {
#endif

int MFNRegionInterior(MFNRegion,MFNVector,MFErrorHandler);
MFNRegion MFNRegionCreateCube(double,double,double,double,double,double,MFErrorHandler);
MFNRegion MFNRegionCreateRectangle(double,double,double,double,MFErrorHandler);
MFNRegion MFNRegionCreateHyperCube(int,double,MFErrorHandler);
MFNRegion MFNRegionCreateHyperCubeByCorners(int,MFNVector,MFNVector,MFErrorHandler);
MFNRegion MFNRegionCreatePendula(int,double,double,int,double,double,double,double,double,MFErrorHandler);
MFNRegion MFNRegionCreatePendulaFour(int,int,double,double,double,double,double,double,double,double,double,MFErrorHandler);
MFNRegion MFNRegionCreateDodecahedronMinusIcosahedron(MFErrorHandler);
MFNRegion MFNRegionCreatePolygonal3dRegion(int,double*,MFErrorHandler);
MFNRegion MFNRegionCreateEdge3dRegion(double*,double*,MFErrorHandler);
void MFFreeNRegion(MFNRegion,MFErrorHandler);
void MFRefNRegion(MFNRegion,MFErrorHandler);
void MFWriteNRegion(FILE*,MFNRegion,MFErrorHandler);
MFNRegion MFReadNRegion(FILE*,MFErrorHandler);

int MFCreateTriangulatedPolygon(int,double*,int,int*,int**,MFErrorHandler);

MFNRegion MFNRegionCreatePolyhedral3dRegion(int,double*,int,int*,int**,MFErrorHandler);
MFNRegion MFNRegionCreateTPBVP(int,int,int,double*,double*,double,double,MFErrorHandler);
MFNRegion MFNRegionCreateForcedOscillator(int,int,double,double,double,double,MFErrorHandler);

MFNRegion MFNRegionCreateBaseClass(char*,MFErrorHandler);
void MFNRegionSetTest(MFNRegion,int (*)(MFNVector,void*,MFErrorHandler),MFErrorHandler);
void MFNRegionSetData(MFNRegion,void*,MFErrorHandler);
void *MFNRegionGetData(MFNRegion,MFErrorHandler);
void MFNRegionSetFreeData(MFNRegion,void (*)(void*,MFErrorHandler),MFErrorHandler);
void MFNRegionSetWriteData(MFNRegion,void (*)(FILE*,void*,MFErrorHandler),MFErrorHandler);


#ifdef __cplusplus
}
#endif

#endif
