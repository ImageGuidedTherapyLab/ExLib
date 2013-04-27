
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(project_policies
  CMP0001 # NEW: CMAKE_BACKWARDS_COMPATIBILITY should no longer be used.
  CMP0002 # NEW: Logical target names must be globally unique.
  CMP0003 # NEW: Libraries linked via full path no longer produce linker search paths.
  CMP0004 # NEW: Libraries linked may NOT have leading or trailing whitespace.
  CMP0005 # NEW: Preprocessor definition values are now escaped automatically.
  CMP0006 # NEW: Installing MACOSX_BUNDLE targets requires a BUNDLE DESTINATION.
  CMP0007 # NEW: List command no longer ignores empty elements.
  CMP0008 # NEW: Libraries linked by full-path must have a valid library file name.
  CMP0009 # NEW: FILE GLOB_RECURSE calls should not follow symlinks by default.
  CMP0010 # NEW: Bad variable reference syntax is an error.
  CMP0011 # NEW: Included scripts do automatic cmake_policy PUSH and POP.
  CMP0012 # NEW: if() recognizes numbers and boolean constants.
  CMP0013 # NEW: Duplicate binary directories are not allowed.
  CMP0014 # NEW: Input directories must have CMakeLists.txt
  )
FOREACH(policy ${project_policies})
  IF(POLICY ${policy})
    CMAKE_POLICY(SET ${policy} NEW)
  ENDIF()
ENDFOREACH()

PROJECT(CRKIT)

# This allows building of old pre-ITK code.
# Current and future code requires ITK.
OPTION( USE_ITK "Use ITK" ON )

# Standard utilities
OPTION(BUILD_COMMAND_LINE_TOOLS "Build command line tools." ON)

# CRL Viz depends on VTK 5.6.0 with QVTK and Qt 4.7.1 or later.
OPTION(BUILD_VISUALIZATION "Build visualization tool." ON)

# CRkit version number.  An even minor number corresponds to releases.
SET(CRKIT_VERSION_MAJOR "1")
SET(CRKIT_VERSION_MINOR "5")
SET(CRKIT_VERSION_PATCH "2")
# Version string should not include patch level.  The major.minor is
# enough to distinguish available features of the toolkit.
# 
# NIW: attempts to execute svnversion to generate a subversion version
# number for the build.  If it fails, it does so silently and the result
# "unknown" is added in the version no.  
EXECUTE_PROCESS(COMMAND "svnversion" "." WORKING_DIRECTORY
    "${CRKIT_SOURCE_DIR}" OUTPUT_VARIABLE SVN_VERSION RESULT_VARIABLE
    SVN_RESULT OUTPUT_STRIP_TRAILING_WHITESPACE)
# non-zero result indicates failure
IF(SVN_RESULT)
    SET(SVN_VERSION "unknown")
ENDIF(SVN_RESULT)
# A parenthesis here causes the Mac OS bundle Info.plist file to fail.
# SET(CRKIT_VERSION_STRING "${CRKIT_VERSION_MAJOR}.${CRKIT_VERSION_MINOR}.${CRKIT_VERSION_PATCH} (${SVN_VERSION})")
SET(CRKIT_VERSION_STRING "${CRKIT_VERSION_MAJOR}.${CRKIT_VERSION_MINOR}.${CRKIT_VERSION_PATCH}")
CONFIGURE_FILE(${CRKIT_SOURCE_DIR}/common/scripts/configuration.h.in ${CRKIT_BINARY_DIR}/configuration.h)

# Provide default options and paths to simplify making a build.
IF(NOT CMAKE_BUILD_TYPE)
SET( CMAKE_BUILD_TYPE Release CACHE STRING 
     "Type of build: Release or Debug or RelWithDebInfo" 
     FORCE
)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET( EXECUTABLE_OUTPUT_PATH . CACHE PATH "Path for binaries created by building to be placed prior to installation." )
SET( LIBRARY_OUTPUT_PATH . CACHE PATH "Path for libraries created by building to be placed prior to installation." )

SET(CMAKE_SKIP_RPATH TRUE)

# The directory External provides tools from third parties and are used
# under the terms of the license provided by the third party.
# These are required for building CRkit.
ADD_SUBDIRECTORY(common/external)

# The directory common/include provides a common library for tools.
# These are required for building CRkit.
#INCLUDE_DIRECTORIES ( common/include common/itkNLOPTOptimizers )
ADD_SUBDIRECTORY(common/include)
ADD_SUBDIRECTORY(common/itkNLOPTOptimizers)
ADD_SUBDIRECTORY(common/crlVTKCommon)

IF(BUILD_VISUALIZATION)
# Software developed at the Computational Radiology Laboratory
# Children's Hospital Boston
  ADD_SUBDIRECTORY(viz)
ENDIF(BUILD_VISUALIZATION)

IF(BUILD_COMMAND_LINE_TOOLS)
# Software developed at the Computational Radiology Laboratory
# Children's Hospital Boston
#
  ADD_SUBDIRECTORY(apps/staple)
  ADD_SUBDIRECTORY(apps/continuous_staple)
  ADD_SUBDIRECTORY(registration)
  ADD_SUBDIRECTORY(tools)
  ADD_SUBDIRECTORY(dti)
ENDIF(BUILD_COMMAND_LINE_TOOLS)


SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Computational Radiology Kit")
SET(CPACK_PACKAGE_VENDOR "Computational Radiology Laboratory, Children's Hospital, Harvard Medical School.")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
# SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../CopyrightLicense.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR ${CRKIT_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${CRKIT_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${CRKIT_VERSION_PATCH})
SET(CPACK_PACKAGE_NAME CRKIT)

IF (BUILD_VISUALIZATION)

# Note that we rely on the version of Qt that is used for QVTK.
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

# import the Qt settings from when VTK was built
# CRL Viz uses QVTK for 2D and 3D drawing
SET(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "")
SET(DESIRED_QT_VERSION ${VTK_DESIRED_QT_VERSION} CACHE FILEPATH "")

FIND_PACKAGE(Qt4 4.7.0)
IF(QT_USE_FILE)
  INCLUDE(${QT_USE_FILE})
ELSE(QT_USE_FILE)
  SET(QT_LIBRARIES   ${QT_QT_LIBRARY})
ENDIF(QT_USE_FILE)

EXECUTE_PROCESS(COMMAND "qmake" "-query" "QT_INSTALL_PREFIX"
    WORKING_DIRECTORY "${CRKIT_BINARY_DIR}" 
    OUTPUT_VARIABLE QT_INSTALL_PREFIX 
    RESULT_VARIABLE QT_INSTALL_PREFIX_RV OUTPUT_STRIP_TRAILING_WHITESPACE
)
# non-zero result indicates failure
IF(QT_INSTALL_PREFIX_RV)
    SET(QT_INSTALL_PREFIX "/tmp/qt-not-found")
ENDIF(QT_INSTALL_PREFIX_RV)
GET_FILENAME_COMPONENT(QT_COMPONENT_NAME ${QT_INSTALL_PREFIX} NAME)
STRING(REGEX MATCH ^qt.*$ QTMATCH ${QT_COMPONENT_NAME})
IF(NOT QTMATCH)
# If the destination doesnt start with qt, add this to the destination.
  SET(QT_COMPONENT_NAME qt/${QT_COMPONENT_NAME})
  SET(QT_DESTINATION_NAME qt)
ELSE(NOT QTMATCH)
# If the destination starts with qt, we don't alter it.
  SET(QT_DESTINATION_NAME "")
ENDIF(NOT QTMATCH)

GET_FILENAME_COMPONENT(VTK_COMPONENT_NAME ${VTK_DIR} NAME)

ENDIF(BUILD_VISUALIZATION)

IF(APPLE)
SET(CPACK_GENERATOR Bundle)
SET(CPACK_BUNDLE_NAME CRKit)

SET(ICON "CRKit")
SET(EXECUTABLE "CRKit")
SET(EXECUTABLE "crlViz")
SET(TYPEINFO "crkt")
CONFIGURE_FILE(${CRKIT_SOURCE_DIR}/common/scripts/crkit-env.in 
  ${CRKIT_BINARY_DIR}/crkit-env.sh )
INSTALL( PROGRAMS ${CRKIT_BINARY_DIR}/crkit-env.sh DESTINATION bin )

IF(BUILD_VISUALIZATION)

CONFIGURE_FILE(${VIZ_SOURCE_DIR}/code/mac/Info.plist.app ${VIZ_BINARY_DIR}/code/crlViz.app/Contents/Info.plist)
SET(CPACK_BUNDLE_PLIST ${VIZ_BINARY_DIR}/code/crlViz.app/Contents/Info.plist)
CONFIGURE_FILE(${VIZ_SOURCE_DIR}/code/mac/StartupCommand.in ${VIZ_BINARY_DIR}/code/StartupCommand)
SET(CPACK_BUNDLE_STARTUP_COMMAND ${VIZ_BINARY_DIR}/code/StartupCommand)
SET(CPACK_BUNDLE_ICON ${VIZ_SOURCE_DIR}/code/mac/staple-128x128.icns )
SET(MACOSX_BUNDLE_ICON_FILE ${VIZ_SOURCE_DIR}/code/mac/staple-128x128.icns )
SET_SOURCE_FILES_PROPERTIES( ${VIZ_SOURCE_DIR}/code/mac/staple-128x128.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources )
SET(plugin_dest_dir CRKit.app/Contents/MacOS)
SET(qtconf_dest_dir CRKit.app/Contents/Resources)

# Maybe my app should be the viz binary version of this...not the installed one.
SET(APPS "${CMAKE_INSTALL_PREFIX}/bin/crlViz")
MESSAGE("BUNDLE APPS is ${APPS}")
MESSAGE("CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX}")

# CMake has limited support for incorporating required Frameworks in the bundle
# so we install the required files by hand.
# Turn off for now, since we will used the fixup_bundle output for this.
#INSTALL(DIRECTORY ${ITK_DIR} DESTINATION Frameworks USE_SOURCE_PERMISSIONS)
#INSTALL(DIRECTORY ${VTK_DIR} DESTINATION Frameworks USE_SOURCE_PERMISSIONS)
#INSTALL(DIRECTORY ${QT_INSTALL_PREFIX} DESTINATION Frameworks/${QT_DESTINATION_NAME} USE_SOURCE_PERMISSIONS)

# If we are using Cocoa and we do a framework build, it works OK.
# If we are using Cocoa and it is not a framework build, we seem to need to
# install the qt_menu.nib by hand. The build script copies it into qt bin dir.
# It is not found there by the launcher.
#    IF QT_MENU.NIB IS A DIR THEN
IF(IS_DIRECTORY ${QT_INSTALL_PREFIX}/bin/qt_menu.nib)
  INSTALL(DIRECTORY ${QT_INSTALL_PREFIX}/bin/qt_menu.nib DESTINATION bin USE_SOURCE_PERMISSIONS)
  INSTALL(CODE "
    file(WRITE \"\${CMAKE_INSTALL_PREFIX}/${qtconf_dest_dir}/qt.conf\" \"\")
  " COMPONENT Runtime )
ENDIF(IS_DIRECTORY ${QT_INSTALL_PREFIX}/bin/qt_menu.nib)

ENDIF(BUILD_VISUALIZATION)

### SIMON FIX:
# Defined by the USE_FILE of the installed VTK, Qt and ITK.
SET(LIBRARY_DIRS 
    ${ITK_DIR}
    ${VTK_LIBRARY_DIRS} 
    ${QT_INSTALL_PREFIX}
    ${CMAKE_INSTALL_PREFIX}/Frameworks/vtk-5.6.1/lib
    ${CMAKE_INSTALL_PREFIX}/Frameworks/qt-4.7.1/lib
    ${CMAKE_INSTALL_PREFIX}/Frameworks/InsightToolkit
)

# Only do this if we are being asked to 'make package'. There is an order
# dependency - this can't run until the crlViz binary exists.
#   Currently working on getting the library dirs and crlViz to configure 
# correctly.
INSTALL(CODE "
  IF(EXISTS ${CMAKE_INSTALL_PREFIX}/bin/crlViz)
    MESSAGE(\"Bundle binary ready - doing bundle fixup.\")
    include(BundleUtilities)
    fixup_bundle(\"${APPS}\" \"\${QTPLUGINS}\" \"${LIBRARY_DIRS}\")
  ENDIF(EXISTS ${CMAKE_INSTALL_PREFIX}/bin/crlViz)
  " COMPONENT Runtime)

ELSEIF(WIN32 AND NOT UNIX)

# NSIS Cpack configuration
INCLUDE(InstallRequiredSystemLibraries)

SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/crlViz.exe")
# Defined by the USE_FILE of the installed VTK, Qt and ITK.
SET(LIBRARY_DIRS ${VTK_LIBRARY_DIRS} ${QT_LIBRARY_DIRS} ${ITK_LIBRARY_DIRS})

#   Currently working on getting the library dirs and crlViz to configure 
# correctly.
INSTALL(CODE "
  include(BundleUtilities)
  fixup_bundle(\"${APPS}\" \"\${QTPLUGINS}\" \"${LIBRARY_DIRS}\")
  " COMPONENT Runtime)

ELSE()
SET(CPACK_GENERATOR TGZ)
CONFIGURE_FILE(${CRKIT_SOURCE_DIR}/common/scripts/crkit-launcher.in 
  ${CRKIT_BINARY_DIR}/crkit-launcher.sh )
CONFIGURE_FILE(${CRKIT_SOURCE_DIR}/common/scripts/crkit-env.in 
  ${CRKIT_BINARY_DIR}/crkit-env.sh )
INSTALL( PROGRAMS ${CRKIT_BINARY_DIR}/crkit-launcher.sh 
  ${CRKIT_BINARY_DIR}/crkit-env.sh DESTINATION bin )

# CMake has limited support for incorporating required Frameworks in the bundle
# so we install the required files by hand.
INSTALL(DIRECTORY ${ITK_DIR} DESTINATION Frameworks USE_SOURCE_PERMISSIONS)

IF(BUILD_VISUALIZATION)
INSTALL(DIRECTORY ${VTK_DIR} DESTINATION Frameworks USE_SOURCE_PERMISSIONS)
INSTALL(DIRECTORY ${QT_INSTALL_PREFIX} 
    DESTINATION Frameworks/${QT_DESTINATION_NAME} USE_SOURCE_PERMISSIONS)
ENDIF(BUILD_VISUALIZATION)

ENDIF()

INCLUDE(CPack)

