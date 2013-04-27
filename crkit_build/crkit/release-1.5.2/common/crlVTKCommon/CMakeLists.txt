
#
# crlVTKCommon 
#
# Copyright (c) 2010-2011 Children's Hospital Boston.
#   Arne Hans, CRL (Computational Radiology Laboratory), Harvard Medical School
#   Benoit Scherrer, CRL (Computational Radiology Laboratory), Harvard Medical School
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

# These utilities rely upon VTK.

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})


#--------------------------------------------------------------------
# SET THE INCLUDE DIRECTORIES
#--------------------------------------------------------------------
INCLUDE_DIRECTORIES( 
    ./ 
    ${CRKIT_SOURCE_DIR}/common/external/gifticlib-1.0.9 
    ${CRKIT_SOURCE_DIR}/common/include
  )

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


#---------------------------------------
# ADD ITKNLOPTOptimizers AS A STATIC LIBRARY
#---------------------------------------
ADD_LIBRARY ( crlVTKCommon
  vtkGIFTIWriter.h
  vtkGIFTIWriter.cxx
  vtkGIFTIReader.h
  vtkGIFTIReader.cxx
  crlVtkMeshIO.h
  crlVtkMeshIO.cxx
)

TARGET_LINK_LIBRARIES( crlVTKCommon vtkIO vtkCommon vtkGraphics giftiio crlCommon ${ZLIB_LIBRARY} ${EXPAT_LIBRARIES} ${NIFTI_LIBRARIES} )

INSTALL(TARGETS crlVTKCommon DESTINATION bin)
