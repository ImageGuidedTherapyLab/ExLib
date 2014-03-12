
PROJECT (ITKTensorResampler)

FIND_PACKAGE(ITK)
IF (NOT ITK_FOUND)
  MESSAGE(FATAL_ERROR "Please set ITK_DIR.")
ENDIF(NOT ITK_FOUND)
INCLUDE(${USE_ITK_FILE})

# Include from the current directory, 
# and assuming a whole tree build which defines the PROJECT CRKIT 
# and assuming a local build from this tree.
# This is to allow the TCLAP API to be included.
# It would be a good idea to set up a variables.cmake file that could be
# sourced to provide the global configuration environment.
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/crlDWICommon
  ${CRKIT_SOURCE_DIR}/common/external/
  ${CRKIT_SOURCE_DIR}/common/include/
  ${CRKIT_TOOLS_BINARY_DIR}
  ${CRKIT_BINARY_DIR}
)

SUBDIRS(crlDWICommon)
SUBDIRS(crlMFMTensorEstimate)


LINK_LIBRARIES(
ITKNumerics
ITKBasicFilters
ITKStatistics
ITKIO
crlCommon
)

ADD_EXECUTABLE (crlTensorResampler crlTensorResamplerMain.cxx ITKTensorResampler.cxx)
TARGET_LINK_LIBRARIES(crlTensorResampler
  ITKIO 
  ITKBasicFilters 
  ITKNumerics 
  ITKCommon
)
INSTALL(TARGETS crlTensorResampler RUNTIME DESTINATION bin)

ADD_EXECUTABLE (crlDeformTensorImage crlDeformTensorImage.cxx TensorDeformableResampler.cxx )
TARGET_LINK_LIBRARIES( crlDeformTensorImage
  ITKIO
  ITKBasicFilters
  ITKNumerics
  ITKCommon
)

ADD_EXECUTABLE (crlTensorLog crlTensorLog.cxx )
TARGET_LINK_LIBRARIES( crlTensorLog
  ITKIO
  ITKBasicFilters
  ITKNumerics
  ITKCommon
)

ADD_EXECUTABLE (crlTensorExp crlTensorExp.cxx )
TARGET_LINK_LIBRARIES( crlTensorExp
  ITKIO
  ITKBasicFilters
  ITKNumerics
  ITKCommon
)

ADD_EXECUTABLE (crlTensorSum crlTensorSum.cxx )
TARGET_LINK_LIBRARIES( crlTensorSum
  ITKIO
  ITKBasicFilters
  ITKNumerics
  ITKCommon
)

ADD_EXECUTABLE (crlTensorScalarMultiply crlTensorScalarMultiply.cxx )
TARGET_LINK_LIBRARIES( crlTensorScalarMultiply
  ITKIO
  ITKBasicFilters
  ITKNumerics
  ITKCommon
)

SET( tensorscalarparameter_SOURCE crlTensorScalarParameter.cxx)
ADD_EXECUTABLE (crlTensorScalarParameter ${tensorscalarparameter_SOURCE} )
TARGET_LINK_LIBRARIES( crlTensorScalarParameter
  ITKIO ITKBasicFilters ITKNumerics ITKCommon
)


ADD_EXECUTABLE ( crlDWIConvertFSLToNHDR
    crlDWIConvertFSLToNHDR.h
    crlDWIConvertFSLToNHDR.cxx
  )
TARGET_LINK_LIBRARIES( crlDWIConvertFSLToNHDR
    crlCommon
    crlDWICommon
    ITKNumerics
    ITKBasicFilters
    ITKIO
  )  


ADD_EXECUTABLE ( crlDWIConvertNHDRForFSL
    crlDWIConvertNHDRForFSL.cxx
    crlDWIConvertNHDRForFSL.h
  )
TARGET_LINK_LIBRARIES( crlDWIConvertNHDRForFSL
    crlCommon
    crlDWICommon
    ITKNumerics
    ITKBasicFilters
    ITKIO
  )  


ADD_EXECUTABLE ( crlDWIConvertNHDRForMedINRIA
    crlDWIConvertNHDRForMedINRIA.cxx
  )
TARGET_LINK_LIBRARIES( crlDWIConvertNHDRForMedINRIA
    crlCommon
    crlDWICommon
    ITKNumerics
    ITKBasicFilters
    ITKIO
  )  

  ADD_EXECUTABLE ( crlDWIRemoveGradientImage
    crlDWIRemoveGradientImage.cxx
    crlDWIRemoveGradientImage.h
  )
  TARGET_LINK_LIBRARIES( crlDWIRemoveGradientImage
    crlCommon
    crlDWICommon
    ITKNumerics
    ITKBasicFilters
    ITKIO
  )  


  ADD_EXECUTABLE ( crlDWIRobustB0
    crlDWIRobustB0.cxx
    crlDWIRobustB0.h
  )
  TARGET_LINK_LIBRARIES( crlDWIRobustB0
    crlCommon
    crlDWICommon
    ITKNumerics
    ITKBasicFilters
    ITKIO
  )

ADD_EXECUTABLE ( crlMFMTensorResample
    crlMFMTensorResample.cxx
    crlMFMTensorResampleClass.h
    crlMFMTensorResampleImageFilter.txx
    crlMFMTensorResampleImageFilter.h
  )
TARGET_LINK_LIBRARIES( crlMFMTensorResample
    crlCommon
    ITKNumerics
    ITKBasicFilters
    ITKIO
  )  


ADD_EXECUTABLE ( crlMFMSimpleModelSelection
    crlMFMSimpleModelSelection.cxx  
    crlDiffusionTensorMeasuresFilter.cxx
    crlDiffusionTensorMeasuresFilter.h    
    crlCCSegmentFilter.h
    crlCCSegmentFilter.cxx
 )

TARGET_LINK_LIBRARIES( crlMFMSimpleModelSelection
    crlCommon
    crlDWICommon
    ITKNLOPTOptimizers
    ITKNumerics
    ITKBasicFilters
    ITKIO
    ITKStatistics
  ) 



ADD_EXECUTABLE ( crlDWIAutoSegmentCC
    crlDWIAutoSegmentCC.cxx
    crlCCSegmentFilter.h
    crlCCSegmentFilter.cxx
 )

TARGET_LINK_LIBRARIES( crlDWIAutoSegmentCC
    crlCommon
    crlDWICommon
    ITKNLOPTOptimizers
    ITKNumerics
    ITKBasicFilters
    ITKIO
    ITKStatistics
  ) 


ADD_EXECUTABLE (crlTensorShapeMeasures crlTensorShapeMeasures.cxx)

ADD_EXECUTABLE (crlTensorToRGB crlTensorToRGB.cxx)

ADD_EXECUTABLE (crlTensorToEigenvector crlTensorToEigenvector.cxx)

ADD_EXECUTABLE (crlVector3ToRGB crlVector3ToRGB.cxx)

ADD_EXECUTABLE( crlDWIToDTIEstimator crlDWIToDTIEstimator.cxx )

ADD_EXECUTABLE (crlTensorClean crlTensorClean.cxx )

ADD_EXECUTABLE (crlTractGenerator crlTractGenerator.cxx )



FIND_PACKAGE(VTK)
IF(VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})

  INCLUDE_DIRECTORIES( 
    ${CRKIT_SOURCE_DIR}/common/crlVTKCommon 
    ${CRKIT_SOURCE_DIR}/common/external/gifticlib-1.0.9  )
  

  #--------------------------------------------------------------------
  # For zlib (because gifti_io.h includes itkzlib/zconf.h)
  #--------------------------------------------------------------------
  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
  # ITK provides expat and zlib.
     INCLUDE(${ITK_USE_FILE})
     SET(EXPAT_LIBRARIES ITKEXPAT)
     SET(ZLIB_LIBRARY itkzlib)
     SET(NIFTI_LIBRARIES ITKniftiio)
  ELSE(ITK_FOUND)
  # Otherwise, get it from elsewhere.
    FIND_PACKAGE(EXPAT)
    FIND_PACKAGE(ZLIB)
    FIND_LIBRARY(NIFTI_LIBRARY NAMES niftiio PATH /usr/local/nifti/lib )
    SET(NIFTI_LIBRARIES ${NIFTI_LIBRARY} znz m)
    FIND_PATH(NIFTI_INCLUDE_DIR nifti1.h /usr/local/nifti/include/nifti )
    GET_FILENAME_COMPONENT(GIFTI_LINK_DIRECTORIES ${NIFTI_LIBRARY} PATH)
    LINK_DIRECTORIES(${GIFTI_LINK_DIRECTORIES})
    INCLUDE_DIRECTORIES(${NIFTI_INCLUDE_DIR})
  ENDIF(ITK_FOUND)

  ADD_EXECUTABLE (crlMFMTractGenerator crlMFMTractGenerator.cxx crlMFMTractGenerator.h crlMFMTractGenerator.txx )
  TARGET_LINK_LIBRARIES( crlMFMTractGenerator ITKIO ITKBasicFilters ITKNumerics ITKCommon crlCommon crlVTKCommon giftiio vtkCommon vtkIO vtkGraphics  )

  ADD_EXECUTABLE (crlMFMWBConnectivity crlMFMWBConnectivity.cxx crlMFMWBConnectivity.h crlMFMWBConnectivity.txx )
  TARGET_LINK_LIBRARIES( crlMFMWBConnectivity ITKIO ITKBasicFilters ITKNumerics ITKCommon crlCommon crlVTKCommon giftiio vtkCommon vtkIO vtkGraphics )

  INSTALL(TARGETS crlMFMTractGenerator crlMFMWBConnectivity RUNTIME DESTINATION bin) 

ELSE(VTK_FOUND)
  MESSAGE(WARNING "VTK not found. Will not compile crlMFMTractGenerator and crlMFMWBConnectivity" )
ENDIF(VTK_FOUND)



INSTALL(TARGETS 
crlTensorLog
crlTensorExp
crlTensorSum
crlTensorScalarMultiply
crlTensorScalarParameter
crlDeformTensorImage
crlTensorShapeMeasures
crlTensorToRGB
crlTensorToEigenvector
crlVector3ToRGB
crlDWIToDTIEstimator
crlTensorClean
crlTractGenerator
crlDWIConvertNHDRForMedINRIA
crlDWIConvertNHDRForFSL
crlDWIConvertFSLToNHDR
crlDWIRobustB0
crlDWIRemoveGradientImage
crlMFMTensorResample
RUNTIME DESTINATION bin)
