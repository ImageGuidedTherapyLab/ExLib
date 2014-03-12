
PROJECT (crlSTAPLE)

 FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
            "Cannot build without ITK. " )
  ENDIF(ITK_FOUND)

INCLUDE_DIRECTORIES(
${ITK_SOURCE_DIR}
${ITK_BINARY_DIR}
${crlSTAPLE_BINARY_DIR}
${CRKIT_SOURCE_DIR}/common/include
${CRKIT_SOURCE_DIR}/common/external
${CRKIT_BINARY_DIR}
)

SET(crlSTAPLE_SRCS
crlSTAPLE.cxx
)

LINK_LIBRARIES(
ITKNumerics
ITKBasicFilters
ITKIO
crlCommon
)

ADD_EXECUTABLE( crlSTAPLE ${crlSTAPLE_SRCS} )
ADD_EXECUTABLE( crlRandomSegmentationImageGenerator crlRandomSegmentationImageGenerator.cxx )

INSTALL(TARGETS crlSTAPLE crlRandomSegmentationImageGenerator
RUNTIME DESTINATION bin )

