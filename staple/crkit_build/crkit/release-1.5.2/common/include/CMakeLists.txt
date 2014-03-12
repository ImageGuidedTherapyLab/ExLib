# DEFINE PROJECT VARIABLES
# - provides CRLCOMMON_BINARY_DIR CRLCOMMON_SOURCE_DIR
PROJECT(CRLCOMMON)

FIND_PACKAGE(ITK)
IF (NOT ITK_FOUND)
  MESSAGE(FATAL_ERROR "Please set ITK_DIR.")
ENDIF(NOT ITK_FOUND)
INCLUDE(${USE_ITK_FILE})


INCLUDE_DIRECTORIES( ../itkNLOPTOptimizers )
#---------------------------------------
# Define the source/header files
#---------------------------------------
SET(CRLCOMMON_HDRS
    crlFileName.h
    crlGeneralizedInstantiator.h
    crlGeneralizedInstantiatorMacros.h
    crlNHDRFileUtils.h
    crlTaskProgress.h
    crlVectorImageUtils.h
    crlTransformUtils.h
    crlDocumentedToolOutput.h
    crlRobustMeanCostFunction.h
    crlRobustMeanCostFunction.cxx
    crlRobustMeanImageFilter.h
    crlRobustMeanImageFilter.txx
    crlFastOps.h
    crlMFMUtils.h
    crlOptimizerUtils.h   
    crlDWIStudy.h
    crlDWIStudy.txx 
    crlDilateLabelImageFilter.h
    crlDilateLabelImageFilter.txx
   )

SET(CRLCOMMON_SRCS
    crlFileName.cxx
    crlNHDRFileUtils.cxx
    crlTaskProgress.cxx
    crlFastOps.cxx
    crlOptimizerUtils.cxx    
   )

#---------------------------------------
# Create the Doxygen  documentation
# (the path of the files must be absolute)
#---------------------------------------
OPTION(CRKIT_COMPILE_DOC "Compile the CRKit classes documentation with Doxygen" OFF )

IF(CRKIT_COMPILE_DOC)
  SET(TMP ${CRLCOMMON_HDRS} ${CRLCOMMON_SRCS} )
  SET(DOXYGEN_INPUT "" )
  FOREACH(I ${TMP})
     SET(DOXYGEN_INPUT "${DOXYGEN_INPUT} ${CMAKE_CURRENT_SOURCE_DIR}/${I}" )
  ENDFOREACH(I)

  SET( DOXYGEN_EXCLUDE  )

  SET( DOXYGEN_TARGET_NAME "doc_crlCommon" )
  SET( DOXYGEN_OUTPUT_DIR "${PROJECT_BINARY_DIR}/docs/crlCommon")
  SET( DOXYGEN_INSTALL_DIR "${PROJECT_BINARY_DIR}/docs/crlCommon" )
  INCLUDE(../../DoxygenDoc.cmake)
ENDIF(CRKIT_COMPILE_DOC)

#---------------------------------------
# Add the library
#---------------------------------------
ADD_LIBRARY (crlCommon STATIC
  ${CRLCOMMON_SRCS}
  ${CRLCOMMON_HDRS}
)
LINK_DIRECTORIES(${CRLCOMMON_BINARY_DIR})
INSTALL(TARGETS crlCommon DESTINATION bin)

