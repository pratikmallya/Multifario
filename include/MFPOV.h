#ifndef __MFPOV_H__
#define __MFPOV_H__

#include <MFErrorHandler.h>

struct MFEnumDualPolytopeSt;

#ifdef __cplusplus
 extern "C" {
#endif

void MFDualPolytopeToPOVFile(char *name,struct MFEnumDualPolytopeSt*,MFErrorHandler);
void MFAtlasToPOV(MFAtlas,char *name,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
