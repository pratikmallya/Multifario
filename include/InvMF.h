#ifndef __InvMF__
#define __InvMF__

#include <IMF.h>
#include <IMFExpansion.h>
#include <IMFFlow.h>
#include <IMFFixedPt.h>
#include <math.h>
#include <IMFSphereOnExpansion.h>
#include <IMFExpansionSpace.h>
#include <IMFExpansionPt.h>
#include <IMFIntegrateFat.h>
#include <IMFInterpolation.h>
#include <IMFFlat.h>

int (*MFIMFGetProjectForSave(MFImplicitMF this,MFErrorHandler e))(MFNVector,double*,void*,MFErrorHandler);
int (*MFIMFGetProjectForDraw(MFImplicitMF this,MFErrorHandler e))(MFNVector,double*,void*,MFErrorHandler);
int (*MFIMFGetProjectForBB  (MFImplicitMF this,MFErrorHandler e))(MFNVector,double*,void*,MFErrorHandler);
#endif
