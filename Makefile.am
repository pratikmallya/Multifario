ACLOCAL_AMFLAGS = -I m4 --install

F77 = $(FC)
FFLAGS = $(FCFLAGS)

SUBDIRS = objs

include_HEADERS= \
                       include/Dodecahedron.h include/ECMsg.h include/ExpCmp.h include/ExpCmpI.h include/IMF.h \
                       include/IMFComputeStableMF.h include/IMFExpansion.h include/IMFExpansionPt.h include/IMFExpansionSpace.h \
                       include/IMFFixedPt.h include/IMFFlat.h include/IMFFlow.h include/IMFIntegrateFat.h include/IMFInterpolation.h \
                       include/IMFNSpace.h include/IMFSphereOnExpansion.h include/IMFThreeDFlows.h include/Icosahedron.h \
                       include/InvMF.h include/MF2dCellComplex.h include/MFAUTO.h include/MFAUTOIMF.h include/MFAtlas.h \
                       include/MFAtlasFriends.h include/MFBO.h include/MFBase.h include/MFBinaryTree.h include/MFChart.h \
                       include/MFContinuationMethod.h include/MFCounters.h include/MFDX.h include/MFDraw.h include/MFDrawChartState.h \
                       include/MFDrawClippedSphere.h include/MFEnumDualPolytope.h include/MFEnumPolytope.h include/MFErrorHandler.h \
                       include/MFFlowOnManifold.h include/MFForced.h include/MFFortran.h include/MFImplicitMF.h include/MFKVector.h \
                       include/MFKuhnTessellation.h include/MFLOCA.h include/MFListOfCharts.h include/MFManifoldComplex.h include/MFMesh.h \
                       include/MFMultifariosMethod.h include/MFNKMatrix.h include/MFNRegion.h include/MFNSpace.h include/MFNVector.h \
                       include/MFPOV.h include/MFPolyhedralComplex.h include/MFPolytope.h include/MFPrint.h include/MFPrintEnum.h \
                       include/MFRheinboldt.h include/MFTPBVP.h include/auto.h include/auto_c.h include/auto_f2c.h \
                       include/auto_mpi.h include/auto_types.h include/multifarioConfig.h include/roma.h include/sh.h include/shInternal.h

#include_HEADERS=
#                      include/Dodecahedron.h include/ECMsg.h include/ExpCmp.h include/ExpCmpI.h \
#                      include/IMF.h include/IMFExpansion.h include/IMFExpansionSpace.h include/IMFFixedPt.h \
#                      include/IMFFlow.h include/IMFIntegrateFat.h include/IMFInterpolation.h include/IMFNSpace.h \
#                      include/IMFSphereOnExpansion.h include/IMFThreeDFlows.h include/Icosahedron.h \
#                      include/MF2dCellComplex.h include/MFAUTO.h include/MFAtlas.h include/MFAtlasFriends.h \
#                      include/MFBO.h include/MFBase.h include/MFBinaryTree.h include/MFChart.h \
#                      include/MFContinuationMethod.h include/MFCounters.h include/MFDX.h include/MFDraw.h \
#                      include/MFDrawChartState.h include/MFDrawClippedSphere.h include/MFEnumDualPolytope.h \
#                      include/MFEnumPolytope.h include/MFErrorHandler.h include/MFForced.h include/MFFortran.h \
#                      include/MFImplicitMF.h include/MFKVector.h include/MFLOCA.h include/MFListOfCharts.h \
#                      include/MFMesh.h include/MFMultifariosMethod.h include/MFNKMatrix.h include/MFNRegion.h \
#                      include/MFNSpace.h include/MFNVector.h include/MFPOV.h include/MFPolytope.h \
#                      include/MFPrint.h include/MFPrintEnum.h include/MFRheinboldt.h include/MFTPBVP.h \
#                      include/auto.h include/auto_c.h include/auto_f2c.h include/auto_mpi.h include/auto_types.h \
#                      include/multifarioConfig.h include/multifarioConfig.h.in include/roma.h include/sh.h \
#                      include/shInternal.h

EXTRA_DIST = bootstrap LICENSE Visualization/genericDX.net Visualization/genericDX.cfg Doc Visualization/Torus.view multifario.dox

if DOXYGEN
include aminclude.am

docs: doxygen-doc

MOSTLYCLEANFILES = $(DX_CLEANFILES)
endif

