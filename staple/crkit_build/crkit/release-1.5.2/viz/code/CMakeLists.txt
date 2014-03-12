
PROJECT (Viz)

SET(VTK_NO_LOAD_CMAKE_EXTENSIONS 1)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${USE_ITK_FILE})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

# import the Qt settings from when VTK was built
# CRL Viz uses the QVTK for 2D and 3D drawing
SET(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "")
SET(DESIRED_QT_VERSION ${VTK_DESIRED_QT_VERSION} CACHE FILEPATH "")

MESSAGE(STATUS "VTK_QT_QMAKE_EXECUTABLE is ${VTK_QT_QMAKE_EXECUTABLE}.")

SET(QT_MIN_VERSION "4.6.3")
SET(QT_REQUIRED TRUE)
FIND_PACKAGE(Qt4)
IF(NOT QT4_FOUND)
    MESSAGE(SEND_ERROR "Qt ${QT_MIN_VERSION} or greater not found.  Please check the QT_QMAKE_EXECUTABLE variable.")
    MESSAGE(SEND_ERROR "qmake QT_QMAKE_EXECUTABLE ${QT_QMAKE_EXECUTABLE} .")
ENDIF(NOT QT4_FOUND)

#FIND_PACKAGE(Qt4 COMPONENTS QtCore QtGui QtOpenGL QtXml REQUIRED)
# FIND_PACKAGE(Qt4 4.5.3 COMPONENTS QtCore QtGui QtOpenGL QtXml REQUIRED)
IF(QT_USE_FILE)
  INCLUDE(${QT_USE_FILE})
ELSE(QT_USE_FILE)
  SET(QT_LIBRARIES   ${QT_QT_LIBRARY})
ENDIF(QT_USE_FILE)

# Some files get generated in the current working directory.
# However, VTK should already have set these variables for us, since it
# needs them for the QVTK widget.
INCLUDE_DIRECTORIES(
  ${QT_INCLUDES}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CRKIT_SOURCE_DIR}/common/crlVTKCommon
  ${CRKIT_SOURCE_DIR}/common/include
  ${CRKIT_SOURCE_DIR}/common/external
  ${CRKIT_SOURCE_DIR}/common/external/gifticlib-1.0.9
  ${CRKIT_SOURCE_DIR}/apps/staple/code
  ${CRKIT_SOURCE_DIR}/tools/code
)

# Some configuration settings and version info are set in the toplevel CRkit 
# directory.
# This provides the file configuration.h in the top level of the build directory
# of the entire toolkit.  This subdirectory cannot be built without this.
INCLUDE_DIRECTORIES(${CRKIT_BINARY_DIR})

OPTION(VIZ_VVLINK "Enable VVLink support in crlViz" OFF)
IF(VIZ_VVLINK)
  FIND_PATH(VVLCLIENTAPI_DIR VVLConnection.h)
  FIND_LIBRARY(VVLCLIENTAPI_LIB NAMES vvlcppapi)
  INCLUDE_DIRECTORIES(${VVLCLIENTAPI_DIR})
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_VVLINK")
ENDIF(VIZ_VVLINK)

IF(ITK_FOUND)
   SET(EXPAT_LIBRARIES ITKEXPAT)
ELSE(ITK_FOUND)
# Otherwise, get it from elsewhere.
  FIND_PACKAGE(EXPAT)
ENDIF(ITK_FOUND)


SET(VIZ_SRCS 
main.cxx 
VizMainWindow.h
VizMainWindow.cxx 
VizViewerWindow.h
VizViewerWindow.cxx
VizFourViewers.h
VizFourViewers.cxx
MultiVizViewerWindow.h
MultiVizViewerWindow.cxx
version.h
version.cxx
VizWizardSTAPLE.h
VizWizardSTAPLE.cxx
VizViewerWindowRotationDialog.h
VizViewerWindowRotationDialog.cxx
)


SET(VizHeaders 
itkImageToVTKImageFilter.h 
itkVTKImageToImageFilter.h
SceneImageView.h
)
SET (UI_SRCS
viz.ui
vizviewer.ui
vizfourviewers.ui
multivizviewer.ui
vizviewerwindowrotationdialog.ui
)

QT4_WRAP_UI(UI_HEADERS ${UI_SRCS})
QT4_WRAP_CPP(MOC_SRCS VizMainWindow.h)
QT4_WRAP_CPP(MOC2_SRCS VizViewerBase.h)
QT4_WRAP_CPP(MOC_SRCS VizViewerWindow.h)
QT4_WRAP_CPP(MOC_SRCS VizFourViewers.h)
QT4_WRAP_CPP(MOC_SRCS MultiVizViewerWindow.h)
QT4_WRAP_CPP(MOC_SRCS VizWizardSTAPLE.h)
QT4_WRAP_CPP(MOC_SRCS VizViewerWindowRotationDialog.h)

ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
SET_SOURCE_FILES_PROPERTIES(${VIZ_SRCS} PROPERTIES OBJECT_DEPENDS "${UI_HEADERS}")
SET_SOURCE_FILES_PROPERTIES(Render2D3D.cxx PROPERTIES OBJECT_DEPENDS "${UI_HEADERS}")

# It's nice to have the ui in the windows project file...just double click on
# it and designer comes up on that ui file :)
IF (${CMAKE_BUILD_TOOL} MATCHES "msdev")
  SET (VIZ_SRCS ${VIZ_SRCS} viz.ui)
ENDIF (${CMAKE_BUILD_TOOL} MATCHES "msdev")
#
#IF(APPLE)
#  # next line due to bug 3648 http://public.kitware.com/ParaView/Bug/bug.php?op=show&bugid=3648&pos=
#  SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}) 
#  ADD_EXECUTABLE (crlViz MACOSX_BUNDLE ${VIZ_SRCS} ${UI_SRCS} ${MOC_SRCS} )
#  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/crlViz.app )
#  INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/crlViz.app DESTINATION bin 
#    DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
#    FILE_PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
#  )
#ELSE(APPLE)
#  ADD_EXECUTABLE (crlViz MACOSX_BUNDLE ${VIZ_SRCS} ${UI_SRCS} ${MOC_SRCS} )
#  INSTALL(TARGETS crlViz DESTINATION bin )
#ENDIF(APPLE)
#
# 

ADD_LIBRARY(VizRender Render2D3D.cxx CrosshairActor.cxx SceneImageView.cxx SceneImage.cxx SceneTensors.cxx SceneModel.cxx SceneRGBImage.cxx ImageOperations.cxx crlTensorGlyph.cxx itkImageToVTKImageFilter.txx itkVTKImageToImageFilter.txx ${MOC2_SRCS} )
TARGET_LINK_LIBRARIES( VizRender crlVTKCommon )
INSTALL(TARGETS VizRender ARCHIVE DESTINATION bin)


ADD_EXECUTABLE (crlViz ${VIZ_SRCS} ${UI_SRCS} ${MOC_SRCS} )

TARGET_LINK_LIBRARIES( crlViz
  QVTK
  ${QT_LIBRARIES}
  VizRender
  ITKIO ITKBasicFilters ITKNumerics ITKCommon
  vtkHybrid
  vtkRendering
  vtkIO
  vtkImaging
  vtkCommon
  vtkWidgets
  vtkGraphics
  crlVTKCommon
  giftiio
  crlCommon
)
IF(VIZ_VVLINK)
  TARGET_LINK_LIBRARIES(crlViz ${VVLCLIENTAPI_LIB})
ENDIF(VIZ_VVLINK)

INSTALL(TARGETS crlViz RUNTIME DESTINATION bin )

ADD_EXECUTABLE( crlCreateTriangleModel crlCreateTriangleModel.cxx )
TARGET_LINK_LIBRARIES( crlCreateTriangleModel 
  ITKIO ITKBasicFilters ITKNumerics ITKCommon
  vtkIO vtkCommon vtkGraphics vtkImaging giftiio crlVTKCommon
)

ADD_EXECUTABLE( crlRotateTriangleModel crlRotateTriangleModel.cxx )
TARGET_LINK_LIBRARIES( crlRotateTriangleModel 
  ITKIO ITKBasicFilters ITKNumerics ITKCommon
  vtkIO vtkCommon vtkGraphics vtkImaging giftiio crlVTKCommon
)

ADD_EXECUTABLE( crlDecimateAndSmoothTriangleModel crlDecimateAndSmoothTriangleModel.cxx )
TARGET_LINK_LIBRARIES( crlDecimateAndSmoothTriangleModel
  vtkIO vtkCommon vtkGraphics vtkImaging crlVTKCommon giftiio
)

ADD_EXECUTABLE( crlTractSelector crlTractSelector.cxx )
TARGET_LINK_LIBRARIES( crlTractSelector
  ITKIO ITKBasicFilters ITKNumerics ITKCommon
  vtkIO vtkCommon vtkGraphics vtkImaging crlVTKCommon giftiio
)

ADD_EXECUTABLE( crlTractDensity crlTractDensity.cxx )
TARGET_LINK_LIBRARIES( crlTractDensity
  ITKIO ITKBasicFilters ITKNumerics ITKCommon
  vtkIO vtkCommon vtkGraphics vtkImaging crlVTKCommon giftiio
)

ADD_EXECUTABLE( crlFascicleStreamlineRatio crlFascicleStreamlineRatio.cxx )
TARGET_LINK_LIBRARIES( crlFascicleStreamlineRatio
  ITKIO ITKBasicFilters ITKNumerics ITKCommon
  vtkIO vtkCommon vtkGraphics vtkImaging crlVTKCommon giftiio
)

ADD_EXECUTABLE( crlColorSurfaceModel crlColorSurfaceModel.cxx )
TARGET_LINK_LIBRARIES( crlColorSurfaceModel
  ITKIO ITKBasicFilters ITKNumerics ITKCommon
  vtkIO vtkCommon vtkGraphics vtkImaging giftiio crlVTKCommon 
)

ADD_EXECUTABLE( crlColorSurfaceModelWithRGBImage crlColorSurfaceModelWithRGBImage.cxx )
TARGET_LINK_LIBRARIES( crlColorSurfaceModelWithRGBImage
  ITKIO ITKBasicFilters ITKNumerics ITKCommon
  vtkIO vtkCommon vtkGraphics vtkImaging giftiio crlVTKCommon
)

ADD_EXECUTABLE(crlScalarImageToRGBImage crlScalarImageToRGBImage.cxx)
TARGET_LINK_LIBRARIES(crlScalarImageToRGBImage ITKIO ITKBasicFilters ITKCommon vtkCommon vtkGraphics vtkImaging)

ADD_EXECUTABLE( crlClampModelScalars crlClampModelScalars.cxx)
TARGET_LINK_LIBRARIES( crlClampModelScalars
  vtkIO vtkCommon vtkGraphics vtkImaging giftiio crlVTKCommon
)

ADD_EXECUTABLE( crlHistogram crlHistogram.cxx )
TARGET_LINK_LIBRARIES( crlHistogram
  ITKIO ITKCommon
  vtkIO vtkCommon vtkGraphics vtkImaging crlVTKCommon giftiio
)

ADD_EXECUTABLE(crlConvertGIFTI crlConvertGIFTI.cxx )
TARGET_LINK_LIBRARIES(crlConvertGIFTI
  vtkIO vtkCommon vtkGraphics giftiio ITKIO crlVTKCommon
)

INSTALL(TARGETS 
crlCreateTriangleModel
crlRotateTriangleModel
crlDecimateAndSmoothTriangleModel
crlColorSurfaceModel
crlTractSelector
crlTractDensity
crlFascicleStreamlineRatio
crlColorSurfaceModelWithRGBImage
crlScalarImageToRGBImage
crlClampModelScalars 
crlHistogram
crlConvertGIFTI
RUNTIME DESTINATION bin BUNDLE DESTINATION bin )

# SET_TARGET_PROPERTIES(crlViz crlCreateTriangleModel crlRotateTriangleModel crlColorSurfaceModel crlColorSurfaceModelWithRGBImage crlDecimateAndSmoothTriangleModel crlClampModelScalars crlHistogram crlConvertGIFTI crlScalarImageToRGBImage PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE )

GET_DIRECTORY_PROPERTY(CRKIT_INCLUDE_DIRS INCLUDE_DIRECTORIES)
GET_DIRECTORY_PROPERTY(CRKIT_LIBRARY_DIRS LINK_DIRECTORIES)
#CONFIGURE_FILE(${CRKIT_SOURCE_DIR}/UseCRKIT.cmake.in ${CRKIT_BINARY_DIR}/UseCRKIT.cmake)
