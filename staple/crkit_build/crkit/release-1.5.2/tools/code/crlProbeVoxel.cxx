
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#include <itkVectorImage.h>

#include <iostream>

int GetImageType (std::string fileName,
                     itk::ImageIOBase::IOPixelType &pixelType,
                     itk::ImageIOBase::IOComponentType &componentType,
                     itk::ImageIOBase::Pointer &imageIO,
                     unsigned int &dimensionality )
{
  imageIO = itk::ImageIOFactory::CreateImageIO( fileName.c_str(), itk::ImageIOFactory::ReadMode );
  if ( imageIO.IsNull() )
    {
    std::cerr <<  "Can't make image IO for image "+fileName+" possibly an unsupported type." << std::endl;
    return EXIT_FAILURE;
    }
  imageIO->SetFileName( fileName );		// necessary?  ITK does it.
  imageIO->ReadImageInformation();

  pixelType = imageIO->GetPixelType();
  componentType = imageIO->GetComponentType();
  dimensionality = imageIO->GetNumberOfDimensions();
  return EXIT_SUCCESS;
}


template < class TImage >
class ProbeVoxelBase
{
public:
  typedef   TImage  ImageType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;

  ProbeVoxelBase(std::string filename, std::vector<unsigned int> coords)
  {
  typename ImageReaderType::Pointer r;

  r = ImageReaderType::New();
  r->SetFileName(filename);

  r->Update();	// should catch when instatiated
  typename ImageType::IndexType index;
  for (unsigned int i = 0; i < index.GetIndexDimension(); i++) 
    {
    index[i] = coords[i];
    }
  if ( r->GetOutput()->GetLargestPossibleRegion().IsInside(index) ) 
    {
    std::cout << "Voxel value : " << r->GetOutput()->GetPixel(index) << 
         std::endl;
    } 
  else
    {
    std::cout << "Index is out of range, you said: " << index << std::endl;
    std::cout << "but the image is only: " << r->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;
    }
  };
  
private:
  ProbeVoxelBase() {};
  ProbeVoxelBase(ProbeVoxelBase&) {};
  ProbeVoxelBase& operator=(ProbeVoxelBase&) {};

};



#define SWITCHP(p,c,d) \
  do { \
  if ( p == itk::ImageIOBase::SCALAR ) { \
  	ProbeVoxelBase<itk::Image<c,d> >(filename, coords); \
  } else {\
  	ProbeVoxelBase<itk::VectorImage<c,d> >(filename, coords);\
  }\
  }while(0)

#define SWITCHPC(p,c,d) \
  do{ \
  switch(c) { \
    case itk::ImageIOBase::UCHAR : \
    case itk::ImageIOBase::CHAR : \
    case itk::ImageIOBase::USHORT : \
    case itk::ImageIOBase::SHORT : \
    case itk::ImageIOBase::UINT :  \
    case itk::ImageIOBase::INT :   \
    case itk::ImageIOBase::LONG :  \
    	SWITCHP(p,long,d); break; \
    case itk::ImageIOBase::ULONG : \
    	SWITCHP(p, unsigned long, d); break; \
    case itk::ImageIOBase::FLOAT : \
        SWITCHP(p, float, d); break; \
    case itk::ImageIOBase::DOUBLE : \
        SWITCHP(p, double, d); break; \
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE : \
    default : break; \
  } \
  } while(0)

#define SWITCHPCD(p,c,d) \
  do { \
  switch(dim) { \
  case 2: SWITCHPC(p,c,2); break; \
  case 3: SWITCHPC(p,c,3); break; \
  case 4: SWITCHPC(p,c,4); break; \
  	} \
  } while(0)

 
int main(int argc, char *argv[])
{
  itk::ImageIOBase::IOPixelType	pixelType;
  itk::ImageIOBase::IOComponentType componentType;
  itk::ImageIOBase::Pointer imageIO;
  unsigned int dim;
  
  if ( argc < 2 )
    {
    std::cerr << "usage: " << argv[0] << " inimage dim1 [dim2..]" << std::endl;
    return(1);
    }
  
  std::string filename(argv[1]);

  if ( GetImageType( filename, pixelType, componentType, imageIO, dim ) )
    {
    return EXIT_FAILURE;
    }
  
  if ( argc != (signed int)( 2 + dim ) )
    {
    std::cerr << "usage: " << argv[0] << "inimage dim1 [dim2..]" << std::endl;
    std::cerr << "image is dimension " << dim << std::endl;
    std::cerr << "but you supplied " << argc - 2 << " coordinates" << std::endl;
    return EXIT_FAILURE;
    }
  
  std::vector<unsigned int> coords(dim);
  for ( unsigned int i = 0; i < dim; i++ )
    {
    if ( ! ( std::istringstream(argv[i+2]) >> coords[i])  )
      {
      std::cerr << "couldn't parse coordinate " << argv[i+2] << std::endl;
      return EXIT_FAILURE;
      }
    }
 
  try 
    {
    SWITCHPCD( pixelType, componentType, dim );	// macro makes all template parameters known at runtime.
    }
  catch( itk::ExceptionObject& err )
    {
    std::cerr << "error processing: " << err << std::endl;
    }
  return EXIT_SUCCESS;

#if 0
  unsigned int dim = 3;
  typedef   float PixelType;

  typedef   itk::Image<PixelType,ImageDimension>  ImageType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;

  typedef ImageType::RegionType RegionType;
  typedef ImageType::SizeType   SizeType;
  typedef ImageType::IndexType  IndexType;

  ImageReaderType::Pointer r;

  r = ImageReaderType::New();
  r->SetFileName(argv[1]);
  r->Update(); // read in the file , could catch errors here...

  ImageType::IndexType index;
  for (unsigned int i = 0; i < index.GetIndexDimension(); i++) {
    index[i] = atoi(argv[2 + i]);
  }

  std::cout << "Voxel value : " << r->GetOutput()->GetPixel(index) << 
     std::endl;

  exit(0); // success
#endif
}

