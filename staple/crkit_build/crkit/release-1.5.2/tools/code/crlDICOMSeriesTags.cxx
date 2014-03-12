
#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

// Read a DICOM series and print the main tags.

int main( int argc, char* argv[] )
{

  if( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0] << " DicomDirectory " << std::endl;
    return EXIT_FAILURE;
    }



//  We instantiate then the type to be used for storing the image once it is
//  read into memory.

  typedef signed short       PixelType;
  const unsigned int         Dimension = 3;
  
  typedef itk::Image< PixelType, Dimension >      ImageType;

// We use the image type for instantiating the series reader type and then we
// construct one object of this class.

  typedef itk::ImageSeriesReader< ImageType >     ReaderType;

  ReaderType::Pointer reader = ReaderType::New();

// A GDCMImageIO object is created and assigned to the reader.
//

  typedef itk::GDCMImageIO       ImageIOType;
    
  ImageIOType::Pointer dicomIO = ImageIOType::New();

  reader->SetImageIO( dicomIO );

// A GDCMSeriesFileNames is declared in order to generate the names of DICOM
// slices. We specify the directory with the \code{SetInputDirectory()} method
// and, in this case, take the directory name from the command line arguments.

  typedef itk::GDCMSeriesFileNames     NamesGeneratorType;
  
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  
  nameGenerator->SetInputDirectory( argv[1] );

// The list of files to read is obtained from the name generator by invoking
// the \code{GetInputFileNames()} method and receiving the results in a
// container of strings. The list of filenames is passed to the reader using
// the \code{SetFileNames()} method.

  typedef std::vector<std::string>    FileNamesContainer;
  FileNamesContainer fileNames = nameGenerator->GetInputFileNames();

  reader->SetFileNames( fileNames );


  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
    }

// ITK internally queries GDCM and obtain all the DICOM tags from the file
// headers. The tag values are stored in the \doxygen{MetaDataDictionary} that
// is a general purpose container for \{key,value\} pairs. The Meta data
// dictionary can be recovered from any ImageIO class by invoking the
// \code{GetMetaDataDictionary()} method.

  typedef itk::MetaDataDictionary   DictionaryType;

  const  DictionaryType & dictionary = dicomIO->GetMetaDataDictionary();

// In this example, we are only interested in the DICOM tags that can be
// represented as strings. We declare therefore a \doxygen{MetaDataObject} of
// string type in order to receive those particular values.

  typedef itk::MetaDataObject< std::string > MetaDataStringType;

// The Meta data dictionary is organized as a container with its corresponding
// iterators. We can therefore visit all its entries by first getting access to
// its \code{Begin()} and \code{End()} methods.
//
  DictionaryType::ConstIterator itr = dictionary.Begin();
  DictionaryType::ConstIterator end = dictionary.End();

// We are now ready for walking through the list of DICOM tags. For this
// purpose we use the iterators that we just declared. At every entry we
// attempt to convert it in to a string entry by using the \code{dynamic\_cast}
// based on RTTI information\footnote{Run Time Type Information}. The
// dictionary is organized like a \code{std::map} structure, we should use
// therefore the \code{first} and \code{second} members of every entry in order
// to get access to the \{key,value\} pairs.

  while( itr != end )
    {
    itk::MetaDataObjectBase::Pointer  entry = itr->second;

    MetaDataStringType::Pointer entryvalue = 
      dynamic_cast<MetaDataStringType *>( entry.GetPointer() ) ;

    if( entryvalue )
      {
      std::string tagkey   = itr->first;
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      std::cout << tagkey <<  " = " << tagvalue << std::endl;
      }

    ++itr;
    }
  
//  It is also possible to query for specific entries instead of reading all of
//  them as we did above. In this case, the user must provide the tag
//  identifier using the standard DICOM encoding. The identifier is stored in a
//  string and used as key on the dictionary. 

    std::string entryId = "0010|0010";

    DictionaryType::ConstIterator tagItr = dictionary.Find( entryId );

    if( tagItr == end )
      {
      std::cerr << "Tag " << entryId;
      std::cerr << " not found in the DICOM header" << std::endl;
      }

// Since the entry may or may not be of string type we must again use a
// \code{dynamic\_cast} in order to attempt to convert it to a string dictionary
// entry. If the conversion is successful, then we can print out its content.

    MetaDataStringType::ConstPointer entryvalue = 
      dynamic_cast<const MetaDataStringType *>( tagItr->second.GetPointer() );

    if( entryvalue )
      {
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      std::cout << "Patient's Name (" << entryId <<  ") ";
      std::cout << " is: " << tagvalue << std::endl;
      }

  return EXIT_SUCCESS;

}

