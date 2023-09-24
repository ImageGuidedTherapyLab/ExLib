/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkHessianToObjectnessMeasureImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkTestingMacros.h"


int
main(int argc, char * argv[])
{
  if (argc < 3)
  {
    std::cerr << "Missing parameters." << std::endl;
    std::cerr << "Usage: " << itkNameOfTestExecutableMacro(argv) << " inputImage"
              << " outputImage [ObjectDimension] [Bright/Dark] [alpha] [beta] [gamma] [radius] [scaleObjectness]" << std::endl;
    return EXIT_FAILURE;
  }

  // Define the dimension of the images
  constexpr unsigned char Dimension = 3;

  using PixelType = float;

  // Declare the types of the images
  using ImageType = itk::Image<PixelType, Dimension>;

  using FileReaderType = itk::ImageFileReader<ImageType>;

  // Declare the type of the recursive Gaussian filter
  using GaussianImageFilterType = itk::HessianRecursiveGaussianImageFilter<ImageType>;

  using HessianImageType = GaussianImageFilterType::OutputImageType;

  // Delcare the type of objectness measure image filter

  using ObjectnessFilterType = itk::HessianToObjectnessMeasureImageFilter<HessianImageType, ImageType>;

  FileReaderType::Pointer imageReader = FileReaderType::New();
  imageReader->SetFileName(argv[1]);

  ITK_TRY_EXPECT_NO_EXCEPTION(imageReader->Update());


  // Create a Gaussian filter
  GaussianImageFilterType::Pointer gaussianFilter = GaussianImageFilterType::New();

  // Create an objectness filter
  ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();

  ITK_EXERCISE_BASIC_OBJECT_METHODS(objectnessFilter, HessianToObjectnessMeasureImageFilter, ImageToImageFilter);


  // Connect the input images
  gaussianFilter->SetInput(imageReader->GetOutput());
  gaussianFilter->SetNormalizeAcrossScale(true);
  objectnessFilter->SetInput(gaussianFilter->GetOutput());

  unsigned int objectDimension = 1; 
  if (argc >  3)
  {
    objectDimension = std::stoi(argv[3]);
  }
  objectnessFilter->SetObjectDimension(objectDimension);

  bool brightObject = true;
  if (argc >  4)
  {
    brightObject = std::stoi(argv[4]);
  }
  objectnessFilter->SetBrightObject(brightObject);

  double alphaValue = 0.5;
  if (argc >  5)
  {
    alphaValue = std::stof(argv[5]);
  }
  objectnessFilter->SetAlpha(alphaValue);

  double betaValue = 0.5;
  if (argc >  6)
  {
    betaValue = std::stof(argv[6]);
  }
  objectnessFilter->SetBeta(betaValue);

  double gammaValue = 5.0;
  if (argc >  7)
  {
    gammaValue = std::stof(argv[7]);
  }
  objectnessFilter->SetGamma(gammaValue);

  
  double sigmaValue = 1.0;
  if (argc >  8)
  {
    sigmaValue = std::stof(argv[8]);
  }
  gaussianFilter->SetSigma(sigmaValue);

  // Set the filter properties
  bool scaleObjectnessMeasure = true;
  if (argc >  9)
  {
    scaleObjectnessMeasure = std::stoi(argv[9]);
  }
  objectnessFilter->SetScaleObjectnessMeasure(scaleObjectnessMeasure); // why?

  gaussianFilter->Print(std::cout);
  objectnessFilter->Print(std::cout);
  

  ITK_TRY_EXPECT_NO_EXCEPTION(objectnessFilter->Update());


  // Write the output image
  using FileWriterType = itk::ImageFileWriter<ImageType>;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName(argv[2]);
  // writer->UseCompressionOn();
  writer->SetInput(objectnessFilter->GetOutput());


  ITK_TRY_EXPECT_NO_EXCEPTION(writer->Update());


  std::cout << "Test finished." << std::endl;
  return EXIT_SUCCESS;
}
