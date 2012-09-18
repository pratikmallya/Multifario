#ifndef __MFKUHNTESSELLATION_H__
#define __MFKUHNTESSELLATION_H__

/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   November 22, 2002
 */
static char *id="@(#) $Id: MFKuhnTessellation.h,v 1.1 2011/07/21 18:33:25 mhender Exp $";

#ifdef __cplusplus
extern "C" {
#endif

struct MFKuhnSimplexSt;
typedef struct MFKuhnSimplexSt *MFKuhnSimplex;

MFKuhnSimplex MFCreateKuhnSimplex(int,int*,int*,MFErrorHandler);
MFKuhnSimplex MFCreateFirstKuhnSimplex(int,MFErrorHandler);
void MFRefKuhnSimplex(MFKuhnSimplex,MFErrorHandler);
void MFFreeKuhnSimplex(MFKuhnSimplex,MFErrorHandler);
void MFKuhnSimplexGetVertex(MFKuhnSimplex,MFNVector,double,int,MFNVector,MFErroHandler);
int MFKuhnSimplicesEqual(MFKuhnSimplex,MFKuhnSimplex);
MFKuhnSimplex MFPivotKuhnSimplexAcrossFace(MFKuhnSimplex,int,MFErrorHandler);

#ifdef __cplusplus
 }
#endif
#endif
