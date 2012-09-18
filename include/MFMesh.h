#ifndef __MFMESH_F__
#define __MFMESH_F__
#include <MFEnumDualPolytope.h>
#include <MFErrorHandler.h>

#ifdef __cplusplus
 extern "C" {
#endif

MFEnumDualPolytope MFGenerateMeshOnPolygon(int,int,double*,int,int*,int,int*,double,MFErrorHandler);
MFEnumDualPolytope MFGenerateMeshInPolyhedron(int nV,double **v,int nE,int **ev, int **evs, int nF, int *nFe, int **fe,int **fes, double R,MFErrorHandler e);
MFEnumDualPolytope MFGenerateMeshOnSphere(int,double*,double,MFErrorHandler);
MFEnumDualPolytope MFGenerateMeshOnSphere3d(int,double*,double,MFErrorHandler);

#ifdef __cplusplus
}
#endif

#endif
