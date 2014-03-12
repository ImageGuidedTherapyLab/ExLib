
PROJECT (CRKIT_TOOLS)

 FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
            "Cannot build CRKIT_TOOLS without ITK. " )
  ENDIF(ITK_FOUND)

INCLUDE_DIRECTORIES(
${ITK_SOURCE_DIR}
${ITK_BINARY_DIR}
${CRKIT_TOOLS_BINARY_DIR}
${CRKIT_SOURCE_DIR}/common/external
${CRKIT_SOURCE_DIR}/common/include
${CRKIT_BINARY_DIR}
)

LINK_LIBRARIES(
ITKNumerics
ITKBasicFilters
ITKStatistics
ITKIO
)

ADD_EXECUTABLE( crlIndexOfMaxComponent crlIndexOfMaxComponent.cxx )
ADD_EXECUTABLE( crlExtractSmallerImageFromImage crlExtractSmallerImageFromImage.cxx )
ADD_EXECUTABLE( crlChangeImageInformation crlChangeImageInformation.cxx )

INSTALL(TARGETS
  crlIndexOfMaxComponent
  crlExtractSmallerImageFromImage
  crlChangeImageInformation
  RUNTIME DESTINATION bin
)

ADD_EXECUTABLE( crlAnisotropicDiffusion crlAnisotropicDiffusion.cxx )
ADD_EXECUTABLE( crlImageStatsLabelled crlImageStatsLabelled.cxx )
ADD_EXECUTABLE( crlImageStats crlImageStats.cxx )
ADD_EXECUTABLE( crlImageStats4d crlImageStats4d.cxx )
ADD_EXECUTABLE( crlImageStats2d crlImageStats2d.cxx )
INSTALL(TARGETS
crlImageStats crlImageStats4d crlImageStats2d crlImageStatsLabelled 
RUNTIME DESTINATION bin
)

ADD_EXECUTABLE( crlExtract3DImageFrom4DImage crlExtract3DImageFrom4DImage.cxx )
ADD_EXECUTABLE( crlExtract2DImageFrom3DImage crlExtract2DImageFrom3DImage.cxx )
ADD_EXECUTABLE( crlExtractROI crlExtractROI.cxx )
ADD_EXECUTABLE( crlExtractROI2d crlExtractROI2d.cxx )
ADD_EXECUTABLE( crlExtractROItensor crlExtractROItensor.cxx )
ADD_EXECUTABLE( crlExtractROIRGB2d crlExtractROIRGB2d.cxx )
ADD_EXECUTABLE( crlInsertROI crlInsertROI.cxx )
INSTALL(TARGETS
crlInsertROI crlExtractROI crlExtract2DImageFrom3DImage 
crlExtract3DImageFrom4DImage crlExtractROI2d crlExtractROItensor 
crlExtractROIRGB2d
RUNTIME DESTINATION bin )

ADD_EXECUTABLE( crlVolumeToSlices crlVolumeToSlices.cxx )
ADD_EXECUTABLE( crlGESlicesToVolume crlGESlicesToVolume.cxx )

ADD_EXECUTABLE( crlCastSymMatFloatToDouble crlCastSymMatFloatToDouble.cxx )
ADD_EXECUTABLE( crlCastSymMatDoubleToFloat crlCastSymMatDoubleToFloat.cxx )
ADD_EXECUTABLE( crlCastVectorSixToFloatSymMat crlCastVectorSixToFloatSymMat.cxx)
ADD_EXECUTABLE( crlNormalizeImageFilter crlNormalizeImageFilter.cxx )
INSTALL(TARGETS
crlCastSymMatFloatToDouble
crlCastVectorSixToFloatSymMat
crlCastSymMatDoubleToFloat 
crlNormalizeImageFilter
RUNTIME DESTINATION bin
)

ADD_EXECUTABLE( crlChangeFormat crlChangeFormat.cxx )
ADD_EXECUTABLE( crlChangeFormat4d crlChangeFormat4d.cxx )
ADD_EXECUTABLE( crlChangeFormat2d crlChangeFormat2d.cxx )
ADD_EXECUTABLE( crlZeroImage crlZeroImage.cxx )
ADD_EXECUTABLE( crlBinaryMorphology crlBinaryMorphology.cxx )
ADD_EXECUTABLE( crlBinaryThreshold crlBinaryThreshold.cxx )
ADD_EXECUTABLE( crlThresholdClipClamp crlThresholdClipClamp.cxx )
ADD_EXECUTABLE( crlConnectedComponentFilter crlConnectedComponentFilter.cxx )
ADD_EXECUTABLE( crlConnectedComponentFilter2d crlConnectedComponentFilter2d.cxx )
ADD_EXECUTABLE( crlRootMeanSquare crlRootMeanSquare.cxx )
ADD_EXECUTABLE( crlIndexOfMaxComponent4d crlIndexOfMaxComponent4d.cxx )
ADD_EXECUTABLE( crlIndexOfMaxComponent3d crlIndexOfMaxComponent3d.cxx )
ADD_EXECUTABLE( crlSetBoundary2d crlSetBoundary2d.cxx )
ADD_EXECUTABLE( crlSignedDanielssonDistanceMap2d crlSignedDanielssonDistanceMap2d.cxx )
ADD_EXECUTABLE( crlSignedDanielssonDistanceMap3d crlSignedDanielssonDistanceMap3d.cxx )
ADD_EXECUTABLE( crlImageAddMultiplyAdd crlImageAddMultiplyAdd.cxx )
ADD_EXECUTABLE( crlImageAddMpyAdd crlImageAddMpyAdd3d.cxx )
ADD_EXECUTABLE( crlImageAddMpyAdd2d crlImageAddMpyAdd2d.cxx )
ADD_EXECUTABLE( crlImageAddMpyAdd4d crlImageAddMpyAdd4d.cxx )
ADD_EXECUTABLE( crlFlipImageFilter crlFlipImageFilter.cxx )
ADD_EXECUTABLE( crlShrinkVectorImage crlShrinkVectorImage.cxx )
ADD_EXECUTABLE( crlConvertVectorImageToImage crlConvertVectorImageToImage.cxx )
ADD_EXECUTABLE( crlConvertImageToVectorImage crlConvertImageToVectorImage.cxx )
ADD_EXECUTABLE( crlConvertVectorImageToImage2d3d crlConvertVectorImageToImage2d3d.cxx )
ADD_EXECUTABLE( crlDemonsReg2d crlDemonsReg2d.cxx )
ADD_EXECUTABLE( crlDemonsReg3d crlDemonsReg3d.cxx )
ADD_EXECUTABLE( crlDeformScalarImage crlDeformScalarImage.cxx )
ADD_EXECUTABLE( crlFixSliceMeans crlFixSliceMeans.cxx )
ADD_EXECUTABLE( crlBoundingBox crlBoundingBox.cxx )
ADD_EXECUTABLE( crlConstructVectorImage crlConstructVectorImage.cxx )
ADD_EXECUTABLE( crlConstructVectorImage2d crlConstructVectorImage2d.cxx )
ADD_EXECUTABLE( crlAtlasSample3 crlAtlasSample3.cxx )

INSTALL(TARGETS
crlChangeFormat crlChangeFormat2d crlChangeFormat4d
crlBinaryThreshold crlThresholdClipClamp crlConnectedComponentFilter
crlRootMeanSquare crlSetBoundary2d crlSignedDanielssonDistanceMap2d 
crlSignedDanielssonDistanceMap3d 
crlImageAddMultiplyAdd crlImageAddMpyAdd crlImageAddMpyAdd2d
crlImageAddMpyAdd4d crlFlipImageFilter crlShrinkVectorImage crlVolumeToSlices 
crlFixSliceMeans crlGESlicesToVolume crlBoundingBox crlIndexOfMaxComponent4d
crlIndexOfMaxComponent3d crlConstructVectorImage crlAtlasSample3 
RUNTIME DESTINATION bin
)

INSTALL(TARGETS
crlBinaryMorphology crlConvertVectorImageToImage crlConvertImageToVectorImage 
crlConvertVectorImageToImage2d3d crlConnectedComponentFilter2d 
crlConstructVectorImage2d
crlDemonsReg2d crlDemonsReg3d crlDeformScalarImage crlAnisotropicDiffusion 
crlZeroImage
RUNTIME DESTINATION bin
)

ADD_EXECUTABLE( crlImageAlgebra crlImageAlgebra.cxx )
IF(WIN32)
  SET_TARGET_PROPERTIES(crlImageAlgebra PROPERTIES COMPILE_FLAGS "/bigobj")
ENDIF(WIN32)
INSTALL(TARGETS crlImageAlgebra
RUNTIME DESTINATION bin )

ADD_EXECUTABLE( crlOverlapStats2d crlOverlapStats2d.cxx )
ADD_EXECUTABLE( crlOverlapstats3d crlOverlapStats3d.cxx )
ADD_EXECUTABLE( crlOverlapstats4d crlOverlapStats4d.cxx )
ADD_EXECUTABLE( crlRelabelImages crlRelabelImages.cxx )
ADD_EXECUTABLE( crlOrientImage crlOrientImage.cxx)

ADD_EXECUTABLE( crlExtractFromVectorImage crlExtractFromVectorImage.cxx)
ADD_EXECUTABLE( crlConstructDisplacementFieldVector crlConstructDisplacementFieldVector.cxx )
ADD_EXECUTABLE( crlChangeFormatDicomSeries crlChangeFormatDicomSeries.cxx )
ADD_EXECUTABLE( crlConvertDICOMMOSAIC crlConvertDICOMMOSAIC.cxx )
ADD_EXECUTABLE( crlDICOMPrintPatientInfo crlDICOMPrintPatientInfo.cxx )
ADD_EXECUTABLE( crlDICOMUpdateTagValue crlDICOMUpdateTagValue.cxx )
ADD_EXECUTABLE( crlDICOMSeriesTags crlDICOMSeriesTags.cxx )
ADD_EXECUTABLE( crlCheckerBoard crlCheckerBoard.cxx )
ADD_EXECUTABLE( crlImageInfo crlImageInfo.cxx )
ADD_EXECUTABLE( crlTissueVolumes crlTissueVolumes.cxx )
ADD_EXECUTABLE( crlSquarifyImage crlSquarifyImage.cxx )
ADD_EXECUTABLE( crlUnSquarifyImage crlUnSquarifyImage.cxx )
ADD_EXECUTABLE( crlProbeVoxel crlProbeVoxel.cxx )
ADD_EXECUTABLE( crlProbeVoxelsLabelled crlProbeVoxelsLabelled.cxx )

INSTALL(TARGETS 
crlOverlapStats2d crlOverlapstats3d crlOverlapstats4d
crlProbeVoxelsLabelled
crlProbeVoxel
crlSquarifyImage crlUnSquarifyImage
crlTissueVolumes
crlRelabelImages
crlOrientImage
crlExtractFromVectorImage
crlConstructDisplacementFieldVector
crlChangeFormatDicomSeries
crlImageInfo
crlConvertDICOMMOSAIC
crlDICOMPrintPatientInfo
crlDICOMUpdateTagValue
crlDICOMSeriesTags
crlCheckerBoard
RUNTIME DESTINATION bin )

ADD_EXECUTABLE( crlMultiEchoDicomToITK crlMultiEchoDicomToITK.cxx )
ADD_EXECUTABLE( crlGaussianSmooth crlGaussianSmooth.cxx )
ADD_EXECUTABLE( crlTransformDWIData crlTransformDWIData.cxx )
ADD_EXECUTABLE( crlMeanFieldMRF crlMeanFieldMRF.cxx )
ADD_EXECUTABLE( crlVectorMeanFieldMRF crlVectorMeanFieldMRF.cxx )
ADD_EXECUTABLE( crlCopyImageWithGeometry crlCopyImageWithGeometry.cxx )
ADD_EXECUTABLE( crlIdentityMeasurementFrame crlIdentityMeasurementFrame.cxx )
ADD_EXECUTABLE( crlIdentityMeasurementFrameFullProcessing crlIdentityMeasurementFrameFullProcessing.cxx )
ADD_EXECUTABLE( crlAdaptiveThreshold1 crlAdaptiveThreshold1.cxx )

ADD_EXECUTABLE( crlMultipleOtsuMaskMaker crlMultipleOtsuMaskMaker.cxx )
ADD_EXECUTABLE( crlMultipleMaskImage crlMultipleMaskImage.cxx )
ADD_EXECUTABLE( crlIndexOfMaxComponentVector crlIndexOfMaxComponentVector.cxx )
ADD_EXECUTABLE( crlCorrelator crlCorrelator.cxx )
ADD_EXECUTABLE( crlSubsampleImage2d crlSubsampleImage2d.cxx )
ADD_EXECUTABLE( crlSubsampleImage3d crlSubsampleImage3d.cxx )
ADD_EXECUTABLE( crlSubsampleImage4d crlSubsampleImage4d.cxx )
ADD_EXECUTABLE( crlMakeAverage crlMakeAverage.cxx )
ADD_EXECUTABLE( crlExportDICOM crlExportDICOM.cxx )

ADD_EXECUTABLE( crlExportDICOMRGB crlExportDICOMRGB.cxx )
ADD_EXECUTABLE( crlImageLaplace crlImageLaplace.cxx crlImageLaplace.h crlImageLaplace.txx )
ADD_EXECUTABLE( crlImageCircle crlImageCircle.cxx )
ADD_EXECUTABLE( crlGradientAnisotropicDiffusion crlGradientAnisotropicDiffusion.cxx )
ADD_EXECUTABLE( crlMaskImage crlMaskImage.cxx )
ADD_EXECUTABLE( crlDensityWeightedStats crlDensityWeightedStats.cxx )
ADD_EXECUTABLE( crlSulciDetection crlSulciDetection.cxx)
ADD_EXECUTABLE( crlRelaxationScheme crlRelaxationScheme.cxx )
ADD_EXECUTABLE( crlFuseRGBImages crlFuseRGBImages.cxx )
ADD_EXECUTABLE(crlCannyEdgeDetection crlCannyEdgeDetection.cxx)
ADD_EXECUTABLE( crlCreateLogOddsVectorAtlasFromTemplate crlCreateLogOddsVectorAtlasFromTemplate.cxx )
ADD_EXECUTABLE( crlVectorImageNorm crlVectorImageNorm.cxx )

INSTALL(TARGETS
crlMultiEchoDicomToITK 
crlGaussianSmooth
crlTransformDWIData
crlMeanFieldMRF
crlVectorMeanFieldMRF
crlCopyImageWithGeometry
crlIdentityMeasurementFrame
crlIdentityMeasurementFrameFullProcessing
crlAdaptiveThreshold1
crlMultipleOtsuMaskMaker
crlMultipleMaskImage
crlIndexOfMaxComponentVector
crlCorrelator
crlSubsampleImage2d crlSubsampleImage3d crlSubsampleImage4d
crlMakeAverage
crlExportDICOM
crlExportDICOMRGB
crlImageLaplace
crlImageCircle
crlGradientAnisotropicDiffusion
crlMaskImage
crlDensityWeightedStats
crlSulciDetection
crlRelaxationScheme
crlFuseRGBImages
crlCannyEdgeDetection
crlCreateLogOddsVectorAtlasFromTemplate
crlVectorImageNorm
RUNTIME DESTINATION bin
)

