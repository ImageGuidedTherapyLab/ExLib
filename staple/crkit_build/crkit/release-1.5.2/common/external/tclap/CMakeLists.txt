PROJECT(TCLAP)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

CONFIGURE_FILE(
  ${TCLAP_SOURCE_DIR}/TCLAPConfig.cmake.in
  ${TCLAP_BINARY_DIR}/tclap/TCLAPConfig.cmake
)

# We are only using it out of the source tree, for the purpose of building
# local files.
# INSTALL(DIRECTORY ${TCLAP_SOURCE_DIR} DESTINATION ${TCLAP_BUILD_DIR}/
#       PATTERN ".svn" EXCLUDE
# )

# This seems to be not necessary, in that FINDPACKAGE triggers off a file
# called TCLAPConfig.cmake.
# CONFIGURE_FILE(
#  ${TCLAP_SOURCE_DIR}/TCLAPConfig.cmake.in
#  ${TCLAP_INSTALL_DIR}/tclap/UseTCLAP.cmake
#)

