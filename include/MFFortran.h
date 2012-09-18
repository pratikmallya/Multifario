#ifndef __MFFORTRAN_H__
#define __MFFORTRAN_H__

#include "multifarioConfig.h"

void CALLDBOSVD(double*,int*,int*,int*,int*,double*,int*,int*,double*,int*,double*,int*,double*,double*,double*,double*,int*,int*,int*);
void CALLDGEEV(char*,char*,int*,double*,int*,double*,double*,double*,int*,double*,int*,double*,int*,int*);
void CALLDGESVD(char*,char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,double*,int*,int*);
void CALLDGETRF(int*,int*,double*,int*,int*,int*);
void CALLDGETRS(char*,int*,int*,double*,int*,int*,double*,int*,int*);
void CALLDGEFA(double*,int*,int*,int*,int*);
void CALLDGESL(double*,int*,int*,int*,double*,int*);

#ifdef __cplusplus
 extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#endif
