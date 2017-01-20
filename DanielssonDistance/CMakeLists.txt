# This is the root ITK CMakeLists file.
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

PROJECT(PBNRR)


# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)


add_executable(DanielssonDistanceMapImageFilter DanielssonDistanceMapImageFilter.cxx)
target_link_libraries(DanielssonDistanceMapImageFilter ${ITK_LIBRARIES})

add_executable(SignedMaurerDistanceMapImageFilter SignedMaurerDistanceMapImageFilter.cxx)
target_link_libraries(SignedMaurerDistanceMapImageFilter ${ITK_LIBRARIES})

add_executable(BresenhamLine BresenhamLine.cxx)
target_link_libraries(BresenhamLine ${ITK_LIBRARIES})
 
add_executable(LineProfileFilter LineProfileFilter.cxx)
target_link_libraries(LineProfileFilter ${ITK_LIBRARIES})

