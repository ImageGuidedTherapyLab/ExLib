
/* Facilitate the selection of a region of interest from the data. */

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImage.h"
#include <itkSymmetricSecondRankTensor.h>
#include <itkDiffusionTensor3D.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkSmartPointer.h>
#include <itkLightObject.h>

int main( int argc, char ** argv )
{
  // Verify the number of parameters in the command line
  if( argc != 9 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile  outputImageFile " << std::endl;
    std::cerr << " startX startY startZ sizeX sizeY sizeZ" << std::endl;
    return EXIT_FAILURE;
    }

  const char * inputFileName  = argv[1];
  const char * outputFileName = argv[2];

  static unsigned int const ImageDimension = 3;
  // typedef double   ComponentType;
  typedef float ComponentType;
  typedef   itk::DiffusionTensor3D<ComponentType>  PixelType;
  typedef   itk::Image< PixelType, ImageDimension >  InputImageType;
  typedef   itk::ImageFileReader<InputImageType> ImageReaderType;
  // This output will be the same as the input
  typedef   itk::Image< PixelType, ImageDimension >  OutputImageType;
  typedef   itk::ImageFileWriter<OutputImageType> ImageWriterType;

  ImageReaderType::Pointer  reader = ImageReaderType::New();
  reader->SetFileName( inputFileName );

  try {
    reader->GenerateOutputInformation();
  } catch (itk::ExceptionObject &excp)
  {
    std::cerr << "Failed to read from file " << inputFileName << std::endl;    std::cerr << "Exception error message is : " << excp << std::endl;
    return EXIT_FAILURE;
  }

  itk::ImageIOBase *imageIO = reader->GetImageIO();
  if (imageIO->GetComponentType() != itk::ImageIOBase::FLOAT) {
    std::cerr << "File component type does not match expected type FLOAT."
              << std::endl;
    std::cerr << "Use crlCastSymMatToDouble to change the type of the components."
              << std::endl;
    return EXIT_FAILURE;
  }

  // Actually read in the image now we know the component type is OK.
  reader->Update();

  typedef itk::RegionOfInterestImageFilter< InputImageType, 
                                            OutputImageType > FilterType;

  FilterType::Pointer filter = FilterType::New();

  OutputImageType::IndexType start;
  start[0] = atoi( argv[3] );
  start[1] = atoi( argv[4] );
  start[2] = atoi( argv[5] );


  OutputImageType::SizeType size;
  size[0] = atoi( argv[6] );
  size[1] = atoi( argv[7] );
  size[2] = atoi( argv[8] );

  OutputImageType::RegionType desiredRegion;
  desiredRegion.SetSize(  size  );
  desiredRegion.SetIndex( start );

  //  Then the region is passed to the filter using SetRegionOfInterest()
  filter->SetRegionOfInterest( desiredRegion );

  ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetFileName( outputFileName );

  filter->SetInput( reader->GetOutput() );
  writer->SetInput( filter->GetOutput() );

  filter->Update();

  // Transfer the measurement frame from the input data to the output data
  // Ensure the measurement frame is transferred to the output as well.
  itk::MetaDataDictionary &inDic = reader->GetOutput()->GetMetaDataDictionary();  std::vector<std::vector<double> > msrFrame(ImageDimension);
  for (unsigned int j = 0; j < ImageDimension; j++) {
    msrFrame[j].resize(ImageDimension);
  }
  itk::ExposeMetaData(inDic, "NRRD_measurement frame", msrFrame);

  itk::MetaDataDictionary &outDic =
                filter->GetOutput()->GetMetaDataDictionary();
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(outDic,
     std::string("NRRD_measurement frame"), msrFrame);

  try 
    { 
    writer->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 


  return EXIT_SUCCESS;
}

