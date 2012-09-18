/*
 *  PROGRAM NAME:  multifario
 *
 *  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
 *  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
 *
 *  Please refer to the LICENSE file in the top directory
 *
 *      author: Mike Henderson mhender@watson.ibm.com
 *      date:   May 20, 2010       Adapted from IMFComputeStableMF
 */

static char *id="@(#) $Id: MFFlowOnManifold.h,v 1.1 2010/05/21 12:22:04 mhender Exp $";

MFAtlas MFCoverFlowOnManifold(MFImplicitMF M, IMFFlow F,char *name, int nInitial, MFNVector *u0, MFKVector p0, MFNRegion Omega, double eps, double dt, double tmax, int maxInterp, int maxCharts, MFErrorHandler e);
