
#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageFileWriter.h"
#include "itkMetaImageIO.h"
#include "itkOrientedImage.h"

int main( int argc, char* argv[] )
{

  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0] << 
      " DicomDirectory outputFileName [outputFileName...]" << std::endl
      << std::endl
      << "if, for instance, three output files are specified, then "
      << std::endl 
      << "it is assumed that this there are 3 echoes in the DICOM series " 
      << std::endl
      << "and they will be written out separately with the first echo "
      << std::endl
      << "to the first file, etc."
      << std::endl;

    return EXIT_FAILURE;
    }

  typedef itk::OrientedImage<unsigned short,3>    ImageType;
  typedef itk::ImageSeriesReader< ImageType >     ReaderType;
  typedef itk::ImageFileWriter< ImageType >       WriterType;
  typedef itk::GDCMImageIO                        InputImageIOType;
//  typedef itk::MetaImageIO                        OutputImageIOType;


  int numOutputs = argc - 2;
    
  try
    {
    // scan the directory for DICOM filenames
    std::cout << "Scanning directory (" << 
      argv[1] << ") for DICOM files...";

    typedef itk::GDCMSeriesFileNames NamesGeneratorType;
    NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
    nameGenerator->SetUseSeriesDetails( true );
    if ( numOutputs != 1 )
      {
      // for multiple outputs restrict based on Echo Number
      std::cout << " (restricting based on echo number)...";
      nameGenerator->AddSeriesRestriction("0018|0086");
      }
    nameGenerator->SetInputDirectory( argv[1] );

    typedef std::vector<std::string >   SeriesIdContainer;
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();

    std::cout << " found " << seriesUID.size() << " series'" << std::endl;

    
#if 0
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while( seriesItr != seriesEnd )
      {
      std::cout << seriesItr->c_str() << std::endl;
      seriesItr++;
      }


  
    typedef std::vector<std::string> fileNamesContainer;
    fileNamesContainer fileNames;
    fileNames = nameGenerator->GetInputFileNames();
    
    const int numFileNames = fileNames.size();

    std::cout << " found " << numFileNames << " files" << std::endl;
#endif



    // make sure that the number of found series matches the number of
    // output files specified
    if ( seriesUID.size() != numOutputs )
      {
      std::cerr << "You specified " << numOutputs << " outputs, but " <<
        " we found " << seriesUID.size() << " series in the directory. " <<
        std::endl << 
        "You must specify one output for each." << std::endl;
      return EXIT_FAILURE;

      }


    // try read/write for each output
    ReaderType::Pointer reader = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();
    InputImageIOType::Pointer dicomIO = InputImageIOType::New();
    reader->SetImageIO( dicomIO );
    writer->SetInput( reader->GetOutput() );

    for ( int thisOutput = 0; thisOutput < numOutputs; thisOutput++ )
      {
      typedef std::vector<std::string> fileNamesContainer;
      fileNamesContainer fileNames;
      fileNames = nameGenerator->GetFileNames( seriesUID[thisOutput] );

      std::cout << "read/write output " << argv[ 2+thisOutput ] << "...";

      reader->SetFileNames( fileNames );
      writer->SetFileName( argv[2 + thisOutput ]);
      writer->Update();

      std::cout << std::endl;
      }
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cerr << ex << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
