//
// Create scalars for a surface model
//

#include <vtkPolyData.h>
#include <vtkProbeFilter.h>
#include <vtkImageData.h>
#include <vtkIdTypeArray.h>
#include <crlVtkMeshIO.h>

#include <itkOrientedImage.h>
#include <itkRGBPixel.h>
#include <itkImageFileReader.h>
#include <itkOrientImageFilter.h>
#include <itkImageToVTKImageFilter.h>

typedef unsigned char InternalPixelType;
typedef itk::OrientedImage<itk::RGBPixel<InternalPixelType>,3> ImageType;
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageToVTKImageFilter<ImageType> ITK2VTKConnectorFilterType;

int main(int argc, char *argv[])
{
  if (argc < 4 || argc > 5) {
    std::cout << "Usage: " << argv[0] << " surfaceModelFile imageFile outputModelFile [replaceBlack]" << std::endl;
    std::cout << "If replaceBlack is set, RGB values of (0,0,0) will be replaced by (200,200,200)." << std::endl;
    return EXIT_FAILURE;
  }

  // Read in the image
  ImageType::Pointer inputImage = ImageType::New();

  // Declare and initialize an image orienter filter
  itk::OrientImageFilter< ImageType, ImageType >::Pointer imageOrienter;
  imageOrienter = itk::OrientImageFilter< ImageType, ImageType >::New();

  // Now load the image off disk and initialize the image data.
  ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName( argv[2] );
  try {
    reader->Update();
  } catch (itk::ExceptionObject & err) {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return false;
  }
  imageOrienter->UseImageDirectionOn();
  imageOrienter->SetDesiredCoordinateOrientationToAxial();
  imageOrienter->SetInput( reader->GetOutput() );
  inputImage = imageOrienter->GetOutput();
  try {
    inputImage->Update();
  } catch (itk::ExceptionObject & err) {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return false;
  }

  // Disconnect the DataObject inputImage from the pipeline.
  inputImage->DisconnectPipeline();

  if (argc == 5)
    {
      itk::ImageRegionIterator<ImageType> it(inputImage,inputImage->GetLargestPossibleRegion());
      itk::RGBPixel<InternalPixelType> p;
      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
	  p = it.Get();
	  if (p[0] == 0 && p[1] == 0 && p[2] == 0)
	    {
	      p[0] = 200;
	      p[1] = 200;
	      p[2] = 200;
	      it.Set(p);
	    }
	}
    }

  ITK2VTKConnectorFilterType::Pointer ITK2VTKConnector = ITK2VTKConnectorFilterType::New();
  ITK2VTKConnector->SetInput( inputImage );

  // Read the surface
  vtkPolyData *surface = crlVtkMeshIO::ReadMesh(argv[1]);

  // Set up a pipeline to create scalars
  vtkProbeFilter *probe = vtkProbeFilter::New();
  probe->SetInputConnection(surface->GetProducerPort());
  probe->SetSource(ITK2VTKConnector->GetOutput());
  probe->Update();
  std::cout << probe->GetValidPoints()->GetNumberOfTuples() << std::endl;

  crlVtkMeshIO::WriteMesh(dynamic_cast<vtkPolyData*>(probe->GetOutput()),argv[3]);

  surface->Delete();
  probe->Delete();

  return EXIT_SUCCESS;

}
