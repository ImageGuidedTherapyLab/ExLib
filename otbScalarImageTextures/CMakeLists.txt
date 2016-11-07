# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.8)
#find_package(CUDA REQUIRED)
 
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(Textures)


# Find OTB.
FIND_PACKAGE(OTB REQUIRED)
IF(OTB_FOUND)
    INCLUDE(${OTB_USE_FILE})
ELSE(OTB_FOUND)
    MESSAGE(FATAL_ERROR "OTB not found. Please set OTB_DIR.")
ENDIF(OTB_FOUND)

PROJECT(TutorialsExamples)

ADD_EXECUTABLE(otbScalarImageToTexturesFilter  otbScalarImageToTexturesFilter.cxx)
TARGET_LINK_LIBRARIES(otbScalarImageToTexturesFilter  OTBCommon OTBIO)


