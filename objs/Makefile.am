ACLOCAL_AMFLAGS = -I m4 --install

F77 = $(FC)
AM_FFLAGS = $(FCFLAGS)


lib_LIBRARIES   = ../lib/libMF.a ../lib/libExpCmp.a ../lib/libSH.a ../lib/libauto2000.a ../lib/libautof2c.a

AM_CFLAGS = -g -O -I${top_srcdir}/include -DNOMFSIMPLEX -DMFALLOWVERBOSE -DMANIFOLD -DHAVE_AUTO $(LAPACK_DEFS) $(LIBTIFF_DEFS)

include:
___lib_libMF_a_SOURCES= ../src/MFChart.c ../src/MFNVector.c ../src/MFDenseNVector.c ../src/MFKVector.c \
                      ../src/MFNKMatrix.c ../src/MFNSpace.c ../src/MFTPBVPNSpace.c ../src/MFWeightedNSpace.c \
                      ../src/MFPolytope.c ../src/MFAtlas.c ../src/MFListOfCharts.c ../src/MFBinaryTree.c \
                      ../src/MFErrorHandler.c ../src/MFNRegion.c ../src/MFCSGBallsNRegion.c ../src/MFImplicitMF.c \
                      ../src/MFPlane.c ../src/MFSphere.c ../src/MFAlgebraic.c ../src/MFCircle.c \
                      ../src/MFFlat.c ../src/MFTorus.c ../src/MFPendula.c ../src/MFTPBVP.c \
                      ../src/MFSwallow.c ../src/MFNSpaceMF.c ../src/MFPrint.c ../src/MFEnumPolytope.c \
                      ../src/MFEnumDualPolytope.c ../src/MF3dPolygonMF.c ../src/MF3dEdgeMF.c ../src/MFVolumeOfChart.c \
                      ../src/MFComputePointsOn3dEdge.c ../src/MFComputePointsOn3dPolygon.c ../src/MFContinuationMethod.c \
                      ../src/MFMultifario.c ../src/MFRheinboldt.c ../src/MFPendulaNRegion.c ../src/MFTPBVPNRegion.c \
                      ../src/MFPolyNRegion.c ../src/MFPrintEnumPolytope.c ../src/IMFComputeStableMF.c ../src/IMFNSpace.c \
                      ../src/IMFExpansion.c ../src/IMFExpansionPt.c ../src/IMFExpansionSpace.c ../src/IMFFixedPt.c \
                      ../src/IMFFlat.c ../src/IMFFlow.c ../src/IMFIntegrateFat.c ../src/IMFInterpolation.c \
                      ../src/IMFSphereOnExpansion.c ../src/IMFThreeDFlows.c ../src/MFForced.c ../src/MFForcedNRegion.c \
                      ../src/MF2dCellComplex.c ../src/MFAUTONVector.c ../src/MFAUTONSpace.c ../src/MFAUTONRegion.c \
                      ../src/MFMesh.c \
                      ../src/MFAUTOIMF.c ../src/MFAUTO.c ../src/IMFIntegrateFat.c \
                      ../src/MFFortran.c \
                      ../src/malloc.c ../src/realloc.c \
                      ../src/MFFlowOnManifold.c \
                      ../Utilities/MFPlotfile.c ../Visualization/MFDraw.c \
                      ../Visualization/MFDrawChartState.c ../Visualization/MFDrawClippedSphere.c \
                      ../Visualization/MFDX.c ../Visualization/MFPOV.c
if HAVELAPACK
  ___lib_libMF_a_SOURCES += ../src/dbofa.f ../src/dbosl.f ../src/dbonv.f ../src/dbosvd.f
endif
if HAVEFITPACK
  ___lib_libMF_a_SOURCES += ../src/MFSpline.c
endif

___lib_libExpCmp_a_SOURCES= ../ExpCmpSrc/CmpExpr.c ../ExpCmpSrc/CpExeCd.c ../ExpCmpSrc/CpObjCd.c ../ExpCmpSrc/CrExpD.c \
                      ../ExpCmpSrc/CrObjCd.c ../ExpCmpSrc/CrTokCd.c ../ExpCmpSrc/ECFunction.c ../ExpCmpSrc/ECcos.c \
                      ../ExpCmpSrc/ECfabs.c ../ExpCmpSrc/ECsin.c ../ExpCmpSrc/ECsqrt.c ../ExpCmpSrc/EDacos.c \
                      ../ExpCmpSrc/EDasin.c ../ExpCmpSrc/EDatan.c ../ExpCmpSrc/EDcos.c ../ExpCmpSrc/EDfabs.c \
                      ../ExpCmpSrc/EDlog.c ../ExpCmpSrc/EDlog10.c ../ExpCmpSrc/EDsqrt.c ../ExpCmpSrc/EDtan.c \
                      ../ExpCmpSrc/EDtanh.c ../ExpCmpSrc/EvalExp.c ../ExpCmpSrc/FrObjCd.c ../ExpCmpSrc/FrTokCd.c \
                      ../ExpCmpSrc/GtErMsg.c ../ExpCmpSrc/GtFnIdNm.c ../ExpCmpSrc/GtFnIdVa.c ../ExpCmpSrc/GtIdNm.c \
                      ../ExpCmpSrc/GtIdTy.c ../ExpCmpSrc/GtInIdNm.c ../ExpCmpSrc/GtInIdVl.c ../ExpCmpSrc/GtMsgPrt.c \
                      ../ExpCmpSrc/GtRlIdNm.c ../ExpCmpSrc/GtRlIdVl.c ../ExpCmpSrc/GtUnIdNm.c ../ExpCmpSrc/IdSet.c \
                      ../ExpCmpSrc/NFnId.c ../ExpCmpSrc/NId.c ../ExpCmpSrc/NInId.c ../ExpCmpSrc/NRlId.c \
                      ../ExpCmpSrc/NUnId.c ../ExpCmpSrc/PrsTokCd.c ../ExpCmpSrc/PrtObjCd.c ../ExpCmpSrc/PrtSt.c \
                      ../ExpCmpSrc/PrtSymTb.c ../ExpCmpSrc/PrtTokCd.c ../ExpCmpSrc/PrtTokTy.c ../ExpCmpSrc/StId.c \
                      ../ExpCmpSrc/StIdFn.c ../ExpCmpSrc/StIdIn.c ../ExpCmpSrc/StIdRl.c ../ExpCmpSrc/StIdUn.c \
                      ../ExpCmpSrc/StMsgPrt.c ../ExpCmpSrc/StMthCn.c ../ExpCmpSrc/StMthFn.c ../ExpCmpSrc/TokExpr.c

___lib_libSH_a_SOURCES= ../shSrc/sh3dmask.c ../shSrc/sharrow.c ../shSrc/shchar.c ../shSrc/shclippg.c ../shSrc/shclr.c ../shSrc/shcolor.c ../shSrc/shcs.c \
                     ../shSrc/shcube.c ../shSrc/shdopg.c ../shSrc/shdraw.c ../shSrc/shend.c ../shSrc/shgetz.c ../shSrc/shinit.c \
                     ../shSrc/shlinc.c ../shSrc/shline.c ../shSrc/shline2s.c ../shSrc/shlineoff.c ../shSrc/shlit.c ../shSrc/shlnonrm.c \
                     ../shSrc/shlss.c ../shSrc/shmask.c ../shSrc/shname.c ../shSrc/shnlit.c ../shSrc/shnpln.c ../shSrc/shpause.c \
                     ../shSrc/shpers.c ../shSrc/shpg.c ../shSrc/shpg2.c ../shSrc/shpgc.c ../shSrc/shpgnrm.c ../shSrc/shpl.c \
                     ../shSrc/shpln.c ../shSrc/shplnrm.c ../shSrc/shpnt.c ../shSrc/shpntc.c ../shSrc/shputps.c ../shSrc/shputtiff.c \
                     ../shSrc/shqlit.c ../shSrc/shqnpln.c ../shSrc/shrefsh.c ../shSrc/shres.c ../shSrc/shsave.c ../shSrc/shscal.c \
                     ../shSrc/shsetp.c ../shSrc/shsphere.c ../shSrc/shsrfp.c ../shSrc/shstr.c ../shSrc/shstrs.c ../shSrc/shsync.c \
                     ../shSrc/shtri.c ../shSrc/shtric.c ../shSrc/shtype.c ../shSrc/shview.c

___lib_libauto2000_a_SOURCES= ../auto2000/autlib1.c ../auto2000/autlib2.c ../auto2000/autlib3.c ../auto2000/autlib4.c \
                        ../auto2000/autlib5.c ../auto2000/conpar.c ../auto2000/dmatrix.c ../auto2000/eispack.c \
                        ../auto2000/reduce.c ../auto2000/setubv.c ../auto2000/worker.c ../auto2000/autoGlobals.c

___lib_libautof2c_a_SOURCES= ../auto2000/cabs.c ../auto2000/d_imag.c ../auto2000/d_lg10.c ../auto2000/d_sign.c \
                        ../auto2000/i_dnnt.c ../auto2000/i_nint.c ../auto2000/pow_dd.c ../auto2000/pow_di.c \
                        ../auto2000/pow_ii.c ../auto2000/r_lg10.c ../auto2000/z_abs.c ../auto2000/z_exp.c \
                        ../auto2000/z_log.c

bin_PROGRAMS =   ../bin/AUTBVP ../bin/AUTEXP ../bin/AUTINT ../bin/AUTINTMultifario ../bin/AUTINTMultifario2d ../bin/Compute3Space \
                 ../bin/Compute4Space ../bin/ComputeCircle  ../bin/ComputeForced ../bin/ComputeLine \
                 ../bin/ComputePlane ../bin/ComputePlaneClip ../bin/ComputeTorus \
                 ../bin/DrawAtlas ../bin/DrawAtlasTS ../bin/DrawDual ../bin/DrawPlotfile ../bin/DualToDX \
                 ../bin/PlotfileDualToDX ../bin/PlotfileToDX ../bin/PlotfileToPOV ../bin/PlotfileToVBM ../bin/RheinboldtSphere \
                 ../bin/AtlasToCenters ../bin/AtlasToDX ../bin/AtlasToDual ../bin/AtlasToPOV ../bin/AutoToDX \
                 ../bin/SlicePlotfile \
                 ../bin/Wunderlich \
                 ../bin/ComputeTempl \
                 ../bin/CoverTorus \
                 ../bin/CoverSimple \
                 ../bin/CoverSphere \
                 ../bin/CoverGenusTwo \
                 ../bin/FranksTestCase

TESTS        =   ../bin/AUTBVP ../bin/AUTEXP ../bin/AUTINT ../bin/AUTINTMultifario ../bin/AUTINTMultifario2d ../bin/Compute3Space \
                 ../bin/Compute4Space ../bin/ComputeCircle  ../bin/ComputeForced ../bin/ComputeLine \
                 ../bin/ComputePlane ../bin/ComputePlaneClip ../bin/ComputeTorus \
                 ../bin/DrawAtlas ../bin/DrawAtlasTS ../bin/DrawDual ../bin/DrawPlotfile ../bin/DualToDX \
                 ../bin/PlotfileDualToDX ../bin/PlotfileToDX ../bin/PlotfileToPOV ../bin/PlotfileToVBM ../bin/RheinboldtSphere \
                 ../bin/AtlasToCenters ../bin/AtlasToDX ../bin/AtlasToDual ../bin/AtlasToPOV ../bin/AutoToDX ../bin/SlicePlotfile \xx


if HAVELAPACK
  bin_PROGRAMS+=../bin/ComputeRod ../bin/ComputeCusp ../bin/ComputeExpression ../bin/ComputeGenusTwo ../bin/ComputeSphere \
                ../bin/ComputeSpherePacking ../bin/ComputeSphereSub ../bin/ComputeTaylor24 ../bin/ComputeTaylorA \
                ../bin/ComputeTranscritical ../bin/ComputeDomokos ../bin/LorenzU0 ../bin/PeriodFourResonance \
                ../bin/ComputeOriol ../bin/MeshPolyhedron

  TESTS       +=../bin/ComputeRod ../bin/ComputeCusp ../bin/ComputeExpression ../bin/ComputeGenusTwo ../bin/ComputeSphere \
                ../bin/ComputeSpherePacking ../bin/ComputeSphereSub ../bin/ComputeTaylor24 ../bin/ComputeTaylorA \
                ../bin/ComputeTranscritical ../bin/ComputeDomokos ../bin/LorenzU0 ../bin/PeriodFourResonance \
                ../bin/ComputeOriol
endif
if HAVEFITPACK
  bin_PROGRAMS+=../bin/ComputeADOTorus
  TESTS       +=../bin/ComputeADOTorus
endif


AM_LDFLAGS  = -L${top_srcdir}/lib 

LDADD = -lMF -lExpCmp -lauto2000 -lautof2c -lSH $(TIFF_LIBRARY) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) -lm

___bin_FranksTestCase_SOURCES =          ../Examples/FranksTestCase.c
___bin_AUTBVP_SOURCES =                 ../Examples/AUTBVP.c
___bin_AUTEXP_SOURCES =                 ../Examples/AUTEXP.c
___bin_AUTINT_SOURCES =                 ../Examples/AUTINT.c
___bin_AUTINTMultifario_SOURCES =       ../Examples/AUTINTMultifario.c
___bin_AUTINTMultifario2d_SOURCES =     ../Examples/AUTINTMultifario2d.c
___bin_Compute3Space_SOURCES =          ../Examples/Compute3Space.c
___bin_Compute4Space_SOURCES =          ../Examples/Compute4Space.c
___bin_ComputeCusp_SOURCES =            ../Examples/ComputeCusp.c
___bin_ComputeExpression_SOURCES =      ../Examples/ComputeExpression.c
___bin_ComputeCircle_SOURCES =          ../Examples/ComputeCircle.c
___bin_ComputeDomokos_SOURCES =         ../Examples/ComputeDomokos.c
___bin_ComputeForced_SOURCES =          ../Examples/ComputeForced.c
___bin_ComputeGenusTwo_SOURCES =        ../Examples/ComputeGenusTwo.c
___bin_ComputeLine_SOURCES =            ../Examples/ComputeLine.c
___bin_ComputePlane_SOURCES =           ../Examples/ComputePlane.c
___bin_ComputePlaneClip_SOURCES =       ../Examples/ComputePlaneClip.c
___bin_ComputeRod_SOURCES =             ../Examples/ComputeRod.c
___bin_ComputeSphere_SOURCES =          ../Examples/ComputeSphere.c
___bin_ComputeSpherePacking_SOURCES =   ../Examples/ComputeSpherePacking.c
___bin_ComputeSphereSub_SOURCES =       ../Examples/ComputeSphereSub.c
___bin_ComputeOriol_SOURCES =           ../Examples/ComputeOriol.c
___bin_ComputeTaylor24_SOURCES =        ../Examples/ComputeTaylor24.c
___bin_ComputeTaylorA_SOURCES =         ../Examples/ComputeTaylorA.c
___bin_ComputeTorus_SOURCES =           ../Examples/ComputeTorus.c
___bin_ComputeTranscritical_SOURCES =   ../Examples/ComputeTranscritical.c
___bin_LorenzU0_SOURCES =               ../Examples/LorenzU0.c
___bin_PeriodFourResonance_SOURCES =    ../Examples/PeriodFourResonance.c
___bin_RheinboldtSphere_SOURCES =       ../Examples/RheinboldtSphere.c
___bin_Wunderlich_SOURCES =             ../Examples/Wunderlich.c
___bin_DrawAtlas_SOURCES =              ../Visualization/DrawAtlas.c
___bin_DrawAtlasTS_SOURCES =            ../Visualization/DrawAtlasTS.c
___bin_DrawDual_SOURCES =               ../Visualization/DrawDual.c
___bin_DrawPlotfile_SOURCES =           ../Visualization/DrawPlotfile.c
___bin_DualToDX_SOURCES =               ../Visualization/DualToDX.c
___bin_PlotfileDualToDX_SOURCES =       ../Visualization/PlotfileDualToDX.c
___bin_PlotfileToDX_SOURCES =           ../Visualization/PlotfileToDX.c
___bin_SlicePlotfile_SOURCES =          ../Utilities/SlicePlotfile.c
___bin_PlotfileToPOV_SOURCES =          ../Visualization/PlotfileToPOV.c
___bin_PlotfileToVBM_SOURCES =          ../Visualization/PlotfileToVBM.c
___bin_AtlasToCenters_SOURCES =         ../Utilities/AtlasToCenters.c
___bin_AtlasToDX_SOURCES =              ../Visualization/AtlasToDX.c
___bin_AtlasToDual_SOURCES =            ../Utilities/AtlasToDual.c
___bin_AtlasToPOV_SOURCES =             ../Visualization/AtlasToPOV.c
___bin_AutoToDX_SOURCES =               ../Visualization/AutoToDX.c
___bin_ComputeTempl_SOURCES =           ../Examples/ComputeTempl.c ../Examples/MFTempl.c
___bin_CoverSphere_SOURCES =            ../Examples/CoverSphere.c
___bin_CoverGenusTwo_SOURCES =          ../Examples/CoverGenusTwo.c
___bin_CoverTorus_SOURCES =             ../Examples/CoverTorus.c
___bin_CoverSimple_SOURCES =             ../Examples/CoverSimple.c
___bin_ComputeADOTorus_SOURCES =        ../Examples/ComputeADOTorus.c ../Utilities/splineInterface.f
___bin_ComputeADOTorus_LDADD =           -lMF -lExpCmp -lauto2000 -lautof2c $(LAPACK_LIBS) $(BLAS_LIBS) -lSH $(TIFF_LIBRARY) $(LIBS) $(FLIBS) -lm
___bin_MeshPolyhedron_SOURCES =         ../Examples/MeshPolyhedron.c
