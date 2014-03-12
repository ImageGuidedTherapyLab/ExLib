PROJECT(Utilities)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


ADD_EXECUTABLE( crlScalarImageKmeansClassifier crlScalarImageKmeansClassifier.cxx )
TARGET_LINK_LIBRARIES(crlScalarImageKmeansClassifier ITKNumerics ITKStatistics ITKIO ITKBasicFilters ITKCommon )
INSTALL_PROGRAMS( /bin FILES crlScalarImageKmeansClassifier )

ADD_EXECUTABLE( crlExpectationMaximizationMixtureModelEstimator crlExpectationMaximizationMixtureModelEstimator.cxx )
INSTALL_PROGRAMS( /bin FILES crlExpectationMaximizationMixtureModelEstimator )

TARGET_LINK_LIBRARIES(crlExpectationMaximizationMixtureModelEstimator
ITKIO ITKBasicFilters ITKCommon ITKStatistics ITKNumerics )

