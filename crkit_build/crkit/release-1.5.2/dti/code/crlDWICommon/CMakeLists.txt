
#---------------------------------------
# Define the source/header files
#---------------------------------------
SET(CRL_DWICOMMON_SRCS
    crlDWICommon.cxx
    crlDWIGradients.cxx
    crlDWIRaw.txx
    )

SET(CRL_DWICOMMON_HDRS
    crlDWICommon.h
    crlDWIGradients.h
    crlDWIGradientDirectionsDef.h
    crlDWIRaw.h
    crlDWIRaw.txx
    crlDWIUtils.h
    )

#---------------------------------------
# Add the library
#---------------------------------------
ADD_LIBRARY (crlDWICommon
  ${CRL_DWICOMMON_SRCS}
  ${CRL_DWICOMMON_HDRS}
)

