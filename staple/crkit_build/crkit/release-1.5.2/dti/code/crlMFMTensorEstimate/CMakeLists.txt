#---------------------------------------
# crlMFMTensorEstimate
#---------------------------------------

INCLUDE_DIRECTORIES( ../../../common/itkNLOPTOptimizers )

If(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mssse3")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

IF(MSVC)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2")
                SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}     /Ox /Ob2 /Oi /Ot /Oy /GT /GL /GF /FD /GS- /fp:fast" CACHE STRING "Flags used by the compiler during release builds." FORCE)
                SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Ox /Ob2 /Oi /Ot /Oy /GT /GL /GF /FD /GS- /fp:fast" CACHE STRING "Flags used by the compiler during release builds." FORCE)
                SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /OPT:REF /OPT:ICF /LTCG"  CACHE STRING "Flags used by the compiler during release builds." FORCE)
                SET(CMAKE_MODULE_LINKER_FLAGS_RELEASE "${CMAKE_MODULE_LINKER_FLAGS_RELEASE} /OPT:REF /OPT:ICF /LTCG"  CACHE STRING "Flags used by the compiler during release builds." FORCE)

ENDIF(MSVC)

ADD_EXECUTABLE ( crlMFMTensorEstimate
    crlMFMTensorEstimate.cxx
    crlMFMTensorEstimateFilter.h
    crlMFMTensorEstimateFilter.txx
    crlMFMTensorEstimateCostFunction.h
    crlMFMTensorEstimateCostFunction.txx
    itkOneTensorEstimateFilter.h
    itkOneTensorEstimateFilter.txx
    GlobalDefinitions.h
  )

TARGET_LINK_LIBRARIES( crlMFMTensorEstimate
    crlCommon
    crlDWICommon
    ITKNumerics
    ITKBasicFilters
    ITKIO
    ITKNLOPTOptimizers
  )  

