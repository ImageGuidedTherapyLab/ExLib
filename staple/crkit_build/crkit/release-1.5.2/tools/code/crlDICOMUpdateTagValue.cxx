//  This illustrates how to read a single DICOM slice and write it back
//  with some changed header information as another DICOM slice. Header
//  Key/Value pairs can be specified on the command line. The keys are defined
//  in the file
//
//  \code{Insight/Utilities/gdcm/Dicts/dicomV3.dic}

// We need to include the
// image reader, image writer, the image, the Meta data dictionary and its
// entries the Meta data objects and the GDCMImageIO. The Meta data dictionary
// is the data container that stores all the entries from the DICOM header once
// the DICOM image file is read into an ITK image.

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkGDCMImageIO.h"

#include <list>
#include <fstream>

int main(int argc, char* argv[])
{

  if( argc < 5 )
    {
      std::cerr << "Usage: " << argv[0] << 
                   " DicomImage OutputDicomImage Entry Value [Entry Value]*\n";
    return 1;
    }


// We declare the image type by selecting a particular pixel type and image
// dimension.
//

  typedef signed short InputPixelType;
  const unsigned int   Dimension = 2;
  typedef itk::Image< InputPixelType, Dimension > InputImageType;

// We instantiate the reader type by using the image type as template
// parameter. An instance of the reader is created and the file name to be read
// is taken from the command line arguments.

  typedef itk::ImageFileReader< InputImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

//
// The GDCMImageIO object is created in order to provide the services for
// reading and writing DICOM files. The newly created image IO class is
// connected to the reader.
//
 
  typedef itk::GDCMImageIO           ImageIOType;
  ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
  reader->SetImageIO( gdcmImageIO );

  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "exception in file reader " << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    return EXIT_FAILURE;
    }

// We take the Meta data dictionary from the image that the reader had loaded
// in memory.

  InputImageType::Pointer inputImage = reader->GetOutput();
  typedef itk::MetaDataDictionary   DictionaryType;
  DictionaryType & dictionary = inputImage->GetMetaDataDictionary();


// Now we access the entries in the Meta data dictionary, and for particular
// key values we assign a new content to the entry. This is done here by taking
// \{key,value\} pairs from the command line arguments. The relevant method is
// the EncapsulateMetaData that takes the dictionary and for a given key
// provided by \code{entryId}, replaces the current value with the content of
// the \code{value} variable. This is repeated for every potential pair present
// in the command line arguments.

  for (int i = 3; i < argc; i+=2)
    {
    std::string entryId( argv[i] );
    std::string value( argv[i+1] );
    itk::EncapsulateMetaData<std::string>( dictionary, entryId, value );
    }

// Now that the Dictionary has been updated, we proceed to save the image. This
// output image will have the modified data associated to its DICOM header.
//
// Using the image type, we instantiate a writer type and construct a writer.
// A short pipeline between the reader and the writer is connected. The
// filename to write is taken from the command line arguments. The image IO
// object is connected to the writer.


  typedef itk::ImageFileWriter< InputImageType >  Writer1Type;

  Writer1Type::Pointer writer1 = Writer1Type::New();

  writer1->SetInput( reader->GetOutput() );
  writer1->SetFileName( argv[2] );
  writer1->SetImageIO( gdcmImageIO );

  try
    {
    writer1->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "exception in file writer " << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;

}

