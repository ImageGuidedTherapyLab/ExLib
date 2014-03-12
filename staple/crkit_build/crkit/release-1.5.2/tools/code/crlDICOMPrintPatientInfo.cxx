
#include "itkImageSeriesReader.h"
#include "itkImageFileReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

#include "crlDICOMExtendDictionary.h"

int main( int argc, char* argv[] )
{

  if( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0] << " DicomFile " <<
              "[tag]*" << std::endl;
    return 1;
    }

  crlDICOMExtendDictionary();

  typedef itk::Image<signed short,3>              ImageType;

  // read a file rather than a series
  typedef itk::ImageFileReader< ImageType >       ReaderType;

  ReaderType::Pointer reader = ReaderType::New();

  typedef itk::GDCMImageIO ImageIOType;

  ImageIOType::Pointer dicomIO = ImageIOType::New();

  reader->SetImageIO( dicomIO );

  reader->SetFileName( argv[1] );

  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex << std::endl;
    return 1;
    }


  const unsigned int length = 2048;

  char patientName[  length ];
  char patientID[    length ];
  char patientSex[   length ];
  char patientAge[   length ];
  char studyDescription[    length ];
  char studyDate[    length ];
  char modality[     length ];
  char manufacturer[ length ];
  char institution[  length ];
  char model[        length ];

  dicomIO->GetPatientName(  patientName  );
  dicomIO->GetPatientID(    patientID    );
  dicomIO->GetPatientSex(   patientSex   );
  dicomIO->GetPatientAge(   patientAge   );
  dicomIO->GetStudyDescription(    studyDescription );
  dicomIO->GetStudyDate(    studyDate    );
  dicomIO->GetModality(     modality     );
  dicomIO->GetManufacturer( manufacturer );
  dicomIO->GetInstitution(  institution  );
  dicomIO->GetModel(        model        );

  std::cout << "Patient Name : " << patientName  << std::endl;
  std::cout << "Patient ID   : " << patientID    << std::endl;
  std::cout << "Patient Sex  : " << patientSex   << std::endl;
  std::cout << "Patient Age  : " << patientAge   << std::endl;
  std::cout << "Study Description : " << studyDescription << std::endl;
  std::cout << "Study Date   : " << studyDate    << std::endl;
  std::cout << "Modality     : " << modality     << std::endl;
  std::cout << "Manufacturer : " << manufacturer << std::endl;
  std::cout << "Institution  : " << institution  << std::endl;
  std::cout << "Model        : " << model        << std::endl;


  // ITK internally queries GDCM and obtain all the DICOM tags from the file
  // headers. The tag values are stored in the \doxygen{MetaDataDictionary} that
  // is a general purpose container for \{key,value\} pairs. The Meta data
  // dictionary can be recovered from any ImageIO class by invoking the
  // GetMetaDataDictionary() method.

  typedef itk::MetaDataDictionary   DictionaryType;
  const  DictionaryType & dictionary = dicomIO->GetMetaDataDictionary();

  // Other tags that we will always report
  std::vector< std::string > tagNames;
  std::vector< std::string > tagNamesDescriptions;
  // Now collect the value for tag: 0008 103e LO 1 Series Description
  tagNames.push_back("0008|103e");
  tagNamesDescriptions.push_back("Series Description ");

  tagNames.push_back("0020|0011");
  tagNamesDescriptions.push_back("Series Number ");

  tagNames.push_back("0018|0024");
  tagNamesDescriptions.push_back("Sequence Name ");

  tagNames.push_back("0018|0050");
  tagNamesDescriptions.push_back("Slice Thickness ");

  tagNames.push_back("0028|0010");
  tagNamesDescriptions.push_back("Number of Rows ");

  tagNames.push_back("0028|0011");
  tagNamesDescriptions.push_back("Number of Columns ");

  tagNames.push_back("0020|0032");
  tagNamesDescriptions.push_back("Image Position (Patient)");

  tagNames.push_back("0020|0037");
  tagNamesDescriptions.push_back("Image Orientation (Patient)");

  tagNames.push_back("0028|0030");
  tagNamesDescriptions.push_back("Pixel Spacing");

  DictionaryType::ConstIterator tagItr;
  DictionaryType::ConstIterator end = dictionary.End();
  typedef itk::MetaDataObject< std::string > MetaDataStringType;
  MetaDataStringType::ConstPointer entryvalue;
  for (unsigned int i = 0; i < tagNames.size(); i++) {
    tagItr = dictionary.Find( tagNames[i] );

    if( tagItr == end ) {
      std::cerr << "Value for tag " << tagNamesDescriptions[i] << " "
                     << tagNames[i];
      std::cerr << " not found in the DICOM header" << std::endl;
    } else {
      // Since the entry may or may not be of string type we must again use a
      // dynamic_cast in order to attempt to convert it to a string dictionary
      // entry. If the conversion is successful, then print out its content.

      entryvalue =
        dynamic_cast<const MetaDataStringType *>( tagItr->second.GetPointer() );
      if ( entryvalue )
        {
        std::string tagvalue = entryvalue->GetMetaDataObjectValue();
        std::cout << tagNamesDescriptions[i] << " " << tagNames[i] <<  ") ";
        std::cout << " : " << tagvalue << std::endl;
        }
      }
    }

  // Look for the value of any tags requested on the command line.
  for (int i = 2; i < argc; i++) {
    std::string entryId = std::string( argv[i] );
    tagItr = dictionary.Find( entryId );
    if( tagItr == end ) {
      std::cout << "tag " << entryId << " not found." << std::endl;
      continue;
    }
    entryvalue =
      dynamic_cast<const MetaDataStringType *>( tagItr->second.GetPointer() );
    if ( entryvalue ) {
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      std::cout << "tag (" << entryId <<  ") ";
      std::cout << " : " << tagvalue << std::endl;
    }
  }

  return 0;

}
