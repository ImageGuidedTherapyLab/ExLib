
// This reads a volume represented by a DICOM series and writes it out
// in a new format.
//

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"

#include <tclap/CmdLine.h>


int main( int argc, char* argv[] )
{

  int ignoreDuplicateData = 0;
  std::string *dicomDirectory = new std::string("");
  std::string *outputFileName = new std::string("");
  std::string *seriesName = new std::string("");

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',"1.0");
      TCLAP::UnlabeledValueArg<std::string> dicomDirectoryArg("dicomDirectory","DICOM Directory",true,"","directory name",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputFileNameArg("outputFileName","Output File Name",true,"","file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> seriesNameArg("seriesName","DICOM Series Name",false,"","series name",cmd);
      TCLAP::ValueArg<int> ignoreDuplicateDataArg("i","ignoreDuplicateData","Reject apparently duplicated file names",false,0,"boolean indicator",cmd);

      cmd.parse(argc,argv);

      if (dicomDirectoryArg.isSet()) dicomDirectory = new std::string(dicomDirectoryArg.getValue());
      if (outputFileNameArg.isSet()) outputFileName = new std::string(outputFileNameArg.getValue());
      if (seriesNameArg.isSet()) seriesName = new std::string(seriesNameArg.getValue());
      if (ignoreDuplicateDataArg.isSet()) ignoreDuplicateData = ignoreDuplicateDataArg.getValue();
    }
  catch (TCLAP::ArgException &e)
    {
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
      exit(1);
    }


  typedef float PixelType;
  const unsigned int      Dimension = 3;

  typedef itk::Image< PixelType, Dimension >         ImageType;


// We use the image type for instantiating the type of the series reader and
// for constructing one object of its type.

  typedef itk::ImageSeriesReader< ImageType >        ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

// A GDCMImageIO object is created and connected to the reader. This object is
// the one that is aware of the internal intricacies of the DICOM format.
//

  typedef itk::GDCMImageIO       ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();

  reader->SetImageIO( dicomIO );


// Now we face one of the main challenges of the process of reading a DICOM
// series. That is, to identify from a given directory the set of filenames
// that belong together to the same volumetric image. Fortunately for us, GDCM
// offers functionalities for solving this problem and we just need to invoke
// those functionalities through an ITK class that encapsulates a communication
// with GDCM classes. This ITK object is the GDCMSeriesFileNames. Conveniently
// for us, we only need to pass to this class the name of the directory where
// the DICOM slices are stored. This is done with the \code{SetDirectory()}
// method. The GDCMSeriesFileNames object will explore the directory and will
// generate a sequence of filenames for DICOM files for one study/series.
// In this example, we also call the \code{SetUseSeriesDetails(true)} function
// that tells the GDCMSereiesFileNames object to use additional DICOM
// information to distinguish unique volumes within the directory.  This is
// useful, for example, if a DICOM device assigns the same SeriesID to
// a scout scan and its 3D volume; by using additional DICOM information
// the scout scan will not be included as part of the 3D volume.  Note that
// \code{SetUseSeriesDetails(true)} must be called prior to calling
// \code{SetDirectory()}.
//

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails( true );

  nameGenerator->SetDirectory( dicomDirectory->c_str() );

  try
    {
    std::cout << std::endl << "The directory: " << std::endl;
    std::cout << std::endl << (*dicomDirectory) << std::endl << std::endl;
    std::cout << "Contains the following DICOM Series: ";
    std::cout << std::endl << std::endl;



// The GDCMSeriesFileNames object first identifies the list of DICOM series
// that are present in the given directory. We receive that list in a reference
// to a container of strings and then we can do things like printing out all
// the series identifiers that the generator had found. Since the process of
// finding the series identifiers can potentially throw exceptions, it is
// wise to put this code inside a try/catch block.
//

    typedef std::vector< std::string >    SeriesIdContainer;

    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();

    if (seriesUID.size() <= 0) {
      std::cout << "No series were identified." << std::endl;
      return EXIT_FAILURE;
    }

    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while( seriesItr != seriesEnd )
      {
      std::cout << seriesItr->c_str() << std::endl;
      seriesItr++;
      }

//
// Given that it is common to find multiple DICOM series in the same directory,
// we must tell the GDCM classes what specific series do we want to read. In
// this example we do this by checking first if the user has provided a series
// identifier in the command line arguments. If no series identifier has been
// passed, then we simply use the first series found during the exploration of
// the directory.
//

    std::string seriesIdentifier;

    if((*seriesName) != "" ) // If no optional series identifier
      {
      seriesIdentifier = seriesName->c_str();
      }
    else
      {
      seriesIdentifier = seriesUID.begin()->c_str();
      }


    std::cout << std::endl << std::endl;
    std::cout << "Now reading series: " << std::endl << std::endl;
    std::cout << seriesIdentifier << std::endl;
    std::cout << std::endl << std::endl;

//
// We pass the series identifier to the name generator and ask for all the
// filenames associated to that series. This list is returned in a container of
// strings by the \code{GetFileNames()} method.
//
// \index{itk::GDCMSeriesFileNames!GetFileNames()}
//

    typedef std::vector< std::string >   FileNamesContainer;
    FileNamesContainer fileNames;

    fileNames = nameGenerator->GetFileNames( seriesIdentifier );

    // if the data is pushed to the scanner twice, it can get the same data
    // but stored with a leading 00x_ file name.
    // We want to exclude such files from further consideration.

    if (ignoreDuplicateData) {
      FileNamesContainer newFileNames;
      for (unsigned int i = 0; i < fileNames.size(); i++) {
        if (
             (fileNames[i].find("/001_") == std::string::npos) &&
             (fileNames[i].find("/002_") == std::string::npos) &&
             (fileNames[i].find("/003_") == std::string::npos)
           ) { // keep the file name
          newFileNames.push_back( fileNames[i] );
        }
      }
      fileNames = newFileNames;
    }

//
// The list of filenames can now be passed to the \doxygen{ImageSeriesReader}
// using the \code{SetFileNames()} method.
//
//  \index{itk::ImageSeriesReader!SetFileNames()}

    reader->SetFileNames( fileNames );

//
// Finally we can trigger the reading process by invoking the \code{Update()}
// method in the reader. This call as usual is placed inside a \code{try/catch}
// block.
//

    try
      {
      reader->Update();
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
      return 1;
      }


//
// At this point, we have a volumetric image in memory that we can access by
// invoking the \code{GetOutput()} method of the reader.
//


//
// We proceed now to save the volumetric image in another file, as specified by
// the user in the command line arguments of this program. Thanks to the
// ImageIO factory mechanism, only the filename extension is needed to identify
// the file format in this case.
//

    typedef itk::ImageFileWriter< ImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();

    writer->SetFileName( outputFileName->c_str() );

    writer->SetInput( reader->GetOutput() );
    writer->UseCompressionOn( );

    std::cout  << "Writing the image as " << std::endl << std::endl;
    std::cout  << (*outputFileName) << std::endl << std::endl;


//
// The process of writing the image is initiated by invoking the
// \code{Update()} method of the writer.
//

    try
      {
      writer->Update();
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
      }
    }
  catch (itk::ExceptionObject &ex)
    {
      std::cout << ex << std::endl;
    return EXIT_FAILURE;
    }


// Note that in addition to writing the volumetric image to a file we could
// have used it as the input for any 3D processing pipeline. Keep in mind that
// DICOM is simply a file format and a network protocol. Once the image data
// has been loaded into memory, it behaves as any other volumetric dataset that
// you could have loaded from any other file format.
//


  return EXIT_SUCCESS;

}
