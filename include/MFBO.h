/* 
    %W%
    %D% %T%
   
  Manifold: Please refer to the LICENSE file in the top directory 

*/
/*      author: Mike Henderson mhender@watson.ibm.com */

#ifndef __MFBO_H__
#define __MFBO_H__

#include <fortran.h>

#ifdef __cplusplus
 extern "C" {
#endif

void F77NAME(dbosl)(double*,int*,int*,int*,int*,double*,int*,int*,double*,int*,double*,int*,double*,double*,double*,double*,int*,int*);
void F77NAME(dbonv)(int*,int*,double*,int*,int*,int*,int*,double*,int*,int*,double*,int*,double*,int*,double*,double*,double*,double*,int*,int*);
void F77NAME(dbofa)(double*,int*,int*,int*,int*,double*,int*,int*,double*,int*,double*,int*,double*,int*,int*);


#ifdef __cplusplus
}
#endif

#endif
