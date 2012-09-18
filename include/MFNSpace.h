/* 
    @(#)MFNSpace.h	1.8
    02/05/24 10:30:49
   
    PROGRAM NAME:  Manifold
   
    (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
    CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
   
    Please refer to the LICENSE file in the top directory

*/
/*      author: Mike Henderson mhender@watson.ibm.com */
/*      date:   November 11, 1997                     */
/*              February 2, 1999   Ported to C        */

#ifndef __MFNSPACE_H__
#define __MFNSPACE_H__
#include <MFBase.h>

#include <stdio.h>
#include <MFNVector.h>
#include <MFErrorHandler.h>

struct MFNSpaceSt;
typedef struct MFNSpaceSt *MFNSpace;

#ifdef __cplusplus
 extern "C" {
#endif

MFNSpace MFCreateNSpace(int,MFErrorHandler);
MFNSpace MFCreateWeightedNSpace(int,int,double,MFErrorHandler);
MFNSpace MFCreateTPBVPNSpace(int,int,int,MFErrorHandler);
MFNSpace MFCreateNSpaceWithExponents(int,int*,MFErrorHandler);
MFNSpace MFCreateWeightedNSpaceWithExponents(int,int,double,int*,MFErrorHandler);
void MFRefNSpace(MFNSpace,MFErrorHandler);
void MFFreeNSpace(MFNSpace,MFErrorHandler);
MFNSpace MFCreateForcedOscillatorNSpace(MFErrorHandler);

double MFNSpaceInner(MFNSpace,MFNVector,MFNVector,MFErrorHandler);
double MFNSpaceDistance(MFNSpace,MFNVector,MFNVector,MFErrorHandler);
double MFNSpaceDistance(MFNSpace,MFNVector,MFNVector,MFErrorHandler);
struct MFNKMatrixSt;
double MFNSpaceTangentDistance(MFNSpace,struct MFNKMatrixSt*,struct MFNKMatrixSt*,MFErrorHandler);

void MFNSpaceDirection(MFNSpace,MFNVector,MFNVector,MFNVector,MFErrorHandler);
void MFNSpaceAdd(MFNSpace,MFNVector,MFNVector,MFNVector,MFErrorHandler);
void MFNSpaceScale(MFNSpace,double,MFNVector,MFNVector,MFErrorHandler);

void MFWriteNSpace(FILE*,MFNSpace,MFErrorHandler);
MFNSpace MFReadNSpace(FILE*,MFErrorHandler);

MFNSpace MFCreateNSpaceBaseClass(char*,MFErrorHandler);

void MFNSpaceSetDistance(MFNSpace,double (*distance)(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler),MFErrorHandler);
void MFNSpaceSetInnerProduct(MFNSpace,double (*inner)(MFNSpace,MFNVector,MFNVector,void*,MFErrorHandler),MFErrorHandler);
void MFNSpaceSetDirection(MFNSpace,void (*distance)(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler),MFErrorHandler);
void MFNSpaceSetAdd(MFNSpace,void (*add)(MFNSpace,MFNVector,MFNVector,MFNVector,void*,MFErrorHandler),MFErrorHandler);
void MFNSpaceSetScale(MFNSpace,void (*scale)(MFNSpace,double,MFNVector,MFNVector,void*,MFErrorHandler),MFErrorHandler);
void MFNSpaceSetWriteData(MFNSpace,void (*writedata)(FILE*,MFNSpace,void*,MFErrorHandler),MFErrorHandler);
void MFNSpaceSetFreeData(MFNSpace,void (*freedata)(void *,MFErrorHandler),MFErrorHandler);
void MFNSpaceSetData(MFNSpace,void*,MFErrorHandler);
void *MFNSpaceGetData(MFNSpace,MFErrorHandler);
char *MFNSpaceGetId(MFNSpace,MFErrorHandler);


#ifdef __cplusplus
}
#endif

#endif
