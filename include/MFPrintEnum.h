/*    @(#)MFPrintEnum.h	1.1
      02/07/24 09:34:50
 
      Author: Mike Henderson
      Date:   July 18, 2002
*/

#ifndef PRINTENUM__H
#define PRINTENUM__H
#include <MFEnumPolytope.h>
#include <MFEnumDualPolytope.h>
#include <MFErrorHandler.h>
#include <stdio.h>

#ifdef __cplusplus
 extern "C" {
#endif

void MFPrintEnumPolytope(FILE*,MFEnumPolytope,MFErrorHandler);
void MFPrintEnumDualPolytope(FILE*,MFEnumDualPolytope,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
