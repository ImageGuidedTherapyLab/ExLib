# This project is intended to be built outside the Insight source tree
PROJECT(ConvertBetweenFileFormats)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

ADD_EXECUTABLE(crlConvertBetweenFileFormats
 itkCommandLineArgumentParser.cxx
 castconvert.cxx
 castconvertScalar.cxx
 castconvertDicomScalarA.cxx
 castconvertDicomScalarB.cxx
 castconvertScalar2DA.cxx
 castconvertScalar2DB.cxx
 castconvertScalar3DA.cxx
 castconvertScalar3DB.cxx
 castconvertScalar4DA.cxx
 castconvertScalar4DB.cxx
 castconvertMultiComponent.cxx
 castconvertMultiComponent2D2CA.cxx
 castconvertMultiComponent2D2CB.cxx
 castconvertMultiComponent2D3CA.cxx
 castconvertMultiComponent2D3CB.cxx
 castconvertMultiComponent2D4CA.cxx
 castconvertMultiComponent2D4CB.cxx
 castconvertMultiComponent3D2CA.cxx
 castconvertMultiComponent3D2CB.cxx
 castconvertMultiComponent3D3CA.cxx
 castconvertMultiComponent3D3CB.cxx
 castconvertMultiComponent3D4CA.cxx
 castconvertMultiComponent3D4CB.cxx
 castconvertMultiComponent4D2CA.cxx
 castconvertMultiComponent4D2CB.cxx
 castconvertMultiComponent4D3CA.cxx
 castconvertMultiComponent4D3CB.cxx
 castconvertMultiComponent4D4CA.cxx
 castconvertMultiComponent4D4CB.cxx
)

TARGET_LINK_LIBRARIES(crlConvertBetweenFileFormats ITKIO ITKCommon ITKBasicFilters)

INSTALL(TARGETS
crlConvertBetweenFileFormats
RUNTIME DESTINATION bin
)

# Testing:
IF( BUILD_TESTING )

  SET(IMAGE_PATH ${ConvertBetweenFileFormats_SOURCE_DIR}/Testing/Data)
  SET(IMAGE_PATH_OUTPUTS ${ConvertBetweenFileFormats_BINARY_DIR})
  # Convert from grayscale mhd to other image file formats
  ADD_TEST(ccvnt_to_png                   ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH}/image_in.mhd -out ${IMAGE_PATH_OUTPUTS}/image_out.png)
  ADD_TEST(ccvnt_to_jpg_unsigned_int      ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH}/image_in.mhd -out ${IMAGE_PATH_OUTPUTS}/image_out.jpg -opct unsigned_int)
  ADD_TEST(ccvnt_to_tiff                  ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH}/image_in.mhd -out ${IMAGE_PATH_OUTPUTS}/image_out.tiff)
  ADD_TEST(ccvnt_to_analyze75             ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH}/image_in.mhd -out ${IMAGE_PATH_OUTPUTS}/image_out.img)
  ADD_TEST(ccvnt_to_bmp_unsigned_char     ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH}/image_in.mhd -out ${IMAGE_PATH_OUTPUTS}/image_out.bmp -opct unsigned_char)
  ADD_TEST(ccvnt_to_nrrd_float            ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH}/image_in.mhd -out ${IMAGE_PATH_OUTPUTS}/image_out.nrrd -opct float)
  ADD_TEST(ccvnt_to_vtk_long              ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH}/image_in.mhd -out ${IMAGE_PATH_OUTPUTS}/image_out.vtk -opct long)
  ADD_TEST(ccvnt_to_gipl                  ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH}/image_in.mhd -out ${IMAGE_PATH_OUTPUTS}/image_out.gipl)
  ADD_TEST(ccvnt_to_nifti_test            ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH}/image_in.mhd -out ${IMAGE_PATH_OUTPUTS}/image_out.nii.gz)

  # Convert back to mhd
  ADD_TEST(ccvnt_from_png                 ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH_OUTPUTS}/image_out.png    -out ${IMAGE_PATH_OUTPUTS}/image_in1.mhd)
  ADD_TEST(ccvnt_from_jpg_unsigned_int    ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH_OUTPUTS}/image_out.jpg    -out ${IMAGE_PATH_OUTPUTS}/image_in2.mhd -opct unsigned_short)
  ADD_TEST(ccvnt_from_tiff                ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH_OUTPUTS}/image_out.tiff   -out ${IMAGE_PATH_OUTPUTS}/image_in3.mhd)
  ADD_TEST(ccvnt_from_analyze75           ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH_OUTPUTS}/image_out.hdr    -out ${IMAGE_PATH_OUTPUTS}/image_in4.mhd)
  ADD_TEST(ccvnt_from_bmp_unsigned_char   ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH_OUTPUTS}/image_out.bmp    -out ${IMAGE_PATH_OUTPUTS}/image_in5.mhd -opct unsigned_short)
  ADD_TEST(ccvnt_from_nrrd_float          ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH_OUTPUTS}/image_out.nrrd   -out ${IMAGE_PATH_OUTPUTS}/image_in6.mhd -opct unsigned_short)
  ADD_TEST(ccvnt_from_vtk_long            ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH_OUTPUTS}/image_out.vtk    -out ${IMAGE_PATH_OUTPUTS}/image_in7.mhd -opct unsigned_short)
  ADD_TEST(ccvnt_from_gipl                ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH_OUTPUTS}/image_out.gipl   -out ${IMAGE_PATH_OUTPUTS}/image_in8.mhd)
  ADD_TEST(ccvnt_from_nifti               ${ConvertBetweenFileFormats_BINARY_DIR}/ConvertBetweenFileFormats -in ${IMAGE_PATH_OUTPUTS}/image_out.nii.gz -out ${IMAGE_PATH_OUTPUTS}/image_in8.mhd)

  SET(IMAGE_COMPARE_TESTS ${CXX_TEST_PATH}/ImageCompareTests)

  # Compare with original png
  ADD_TEST(compare_png                    ${IMAGE_COMPARE_TESTS} --compare -in ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in1.mhd Dummy)
 #ADD_TEST(compare_jpg_unsigned_int       ${IMAGE_COMPARE_TESTS} -in ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in2.mhd)
  ADD_TEST(compare_tiff                   ${IMAGE_COMPARE_TESTS} --compare -in ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in3.mhd Dummy)
  ADD_TEST(compare_analyze75              ${IMAGE_COMPARE_TESTS} --compare -in ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in4.mhd Dummy)
  ADD_TEST(compare_bmp_unsigned_char      ${IMAGE_COMPARE_TESTS} --compare -in ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in5.mhd Dummy)
  ADD_TEST(compare_nrrd_float             ${IMAGE_COMPARE_TESTS} --compare -in ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in6.mhd Dummy)
  ADD_TEST(compare_vtk_long               ${IMAGE_COMPARE_TESTS} --compare -in ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in7.mhd Dummy)
  ADD_TEST(compare_gipl                   ${IMAGE_COMPARE_TESTS} --compare -in ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in8.mhd Dummy)
  ADD_TEST(ccvnt_compare_nifti            ${IMAGE_COMPARE_TESTS} --compare -in ${IMAGE_PATH}/image_in.mhd ${IMAGE_PATH_OUTPUTS}/image_in8.mhd Dummy)

  ADD_EXECUTABLE(ImageCompareTests ImageCompareTests.cxx)
  TARGET_LINK_LIBRARIES(ImageCompareTests ITKCommon ITKIO)
ENDIF( BUILD_TESTING )
