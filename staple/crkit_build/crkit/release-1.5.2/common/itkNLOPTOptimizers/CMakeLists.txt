
PROJECT(ITKNLOPT)

#
# ITK Wrapper for the NLOPT library
#
# Copyright (c) 2010-2011 Children's Hospital Boston.
# Benoit Scherrer, CRL (Computational Radiology Laboratory), Harvard Medical School
#
# This software is licensed by the copyright holder under the terms of the
# Open Software License version 3.0.
# http://www.opensource.org/licenses/osl-3.0.php
#
# Attribution Notice.
#
# This research was carried out in the Computational Radiology Laboratory of
# Children's Hospital, Boston and Harvard Medical School.
# http://www.crl.med.harvard.edu
# For more information contact: simon.warfield@childrens.harvard.edu
#
# This research work was made possible by Grant Number R01 RR021885 (Principal
# Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
# from the National Center for Research Resources (NCRR), a component of the
# National Institutes of Health (NIH).
#

# These utilities rely upon ITK.
FIND_PACKAGE(ITK)
IF (NOT ITK_FOUND)
  MESSAGE(FATAL_ERROR "Please set ITK_DIR.")
ENDIF(NOT ITK_FOUND)
INCLUDE(${USE_ITK_FILE})


#--------------------------------------------------------------------
# SET THE INCLUDE DIRECTORIES
#--------------------------------------------------------------------
INCLUDE_DIRECTORIES( ./ 
./nlopt-2.2/api 
./nlopt-2.2/util
${ITKNLOPT_BINARY_DIR} )

#--------------------------------------------------------------------
# SET THE SUBDIRECTORIES
#--------------------------------------------------------------------
SUBDIRS(nlopt-2.2)

#---------------------------------------
# ADD ITKNLOPTOptimizers AS A STATIC LIBRARY
#---------------------------------------
ADD_LIBRARY ( ITKNLOPTOptimizers
  itkNLOPTOptimizers.h
  itkNLOPTOptimizers.cxx
)

TARGET_LINK_LIBRARIES(ITKNLOPTOptimizers nlopt)
  

#---------------------------------------
# TEST PROGRAM
#---------------------------------------
ADD_EXECUTABLE ( ITKNLOPTOptimizers_test
    test.cxx
  )

TARGET_LINK_LIBRARIES( ITKNLOPTOptimizers_test
    ITKNLOPTOptimizers
    ITKNumerics
    ITKBasicFilters
    ITKIO
  )


#---------------------------------------
# DOCUMENTATION
#---------------------------------------
OPTION(ITKNLOPT_COMPILE_DOC "Compile the itkNLOPTOptimizers documentation with Doxygen" OFF )

IF(ITKNLOPT_COMPILE_DOC)
  SET(DOXYGEN_INPUT "${CMAKE_CURRENT_SOURCE_DIR}/itkNLOPTOptimizers.cxx ${CMAKE_CURRENT_SOURCE_DIR}/itkNLOPTOptimizers.h" )
  SET( DOXYGEN_EXCLUDE  )

  SET( DOXYGEN_TARGET_NAME "doc_ITKNLOPTOptimizers" )
  SET( DOXYGEN_OUTPUT_DIR "${ITKNLOPT_BINARY_DIR}/docs/itkNLOPTOptimizers")
  SET( DOXYGEN_INSTALL_DIR "${ITKNLOPT_BINARY_DIR}/docs/itkNLOPTOptimizers" )
  SET( PROJECT_NAME "itkNLOPTOptimizers" )
 

  INCLUDE(../../DoxygenDoc.cmake)
ENDIF(ITKNLOPT_COMPILE_DOC)
 
