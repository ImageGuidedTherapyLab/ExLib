
// This reads a volume represented by a DICOM series and writes it out
// in a new format.
//

#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkImageSeriesReader.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <set>
#include <ctype.h>

#include "configuration.h"
#include "tclap/CmdLine.h"
#include "crlConvertDICOMMOSAIC.txx"

#include "crlDICOMExtendDictionary.h"

int main( int argc, char* argv[] )
{
  crlDICOMExtendDictionary();
  // TCLAP command line argument parsing

  std::string *dicomDirectory = new std::string(".");
  std::string *outputPrefix = new std::string("out");
  std::string *seriesID = new std::string(""); // optional
  std::vector<std::string> restrictionTags;
  int ignoreDuplicateData = 1;
  std::string *outputSuffix = new std::string(".nrrd");
  int processNonMOSAICData = 1;
  int recursiveProcessing = 1;
  int fMRIProcessing = 0;
  std::string mosaicOrientation = "axial";
  int color = 0;

  try {
    // Define the command line object, and insert a message
    // that describes the program.
    TCLAP::CmdLine cmd("Computational Radiology Laboratory program to convert from DICOM format. Some useful tags to restrict the set of images included in a volume are: Series Number 0020|0011, Series Date 0008|0021, Series Time 0008|0031, Series Description 0008|103e, Echo Number 0018|0086, Instance Number 0020|0013, Rescale Intercept 0028|1052 and Rescale Slope 0028|1053", ' ',
		       CRKIT_VERSION_STRING );

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> dicomArg("d", "dicomDir",
					  "Input DICOM Directory Name", true,
					  "/dev/null", "file name",
					  cmd);

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> outputPrefixArg("p", "outputPrefix",
						 "Output Prefix to File Name", false,
						 "out", "file name",
						 cmd);

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> outputSuffixArg("s", "outputSuffix",
						 "Output Suffix to File Name", false,
						 ".nrrd", "file name extension",
						 cmd);

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<std::string> seriesIDArg("n", "seriesName",
					     "Identifier for a particular series", false,
					     "", "DICOM series identifier string",
					     cmd);

    // Define a value argument and add it to the command line.
    TCLAP::ValueArg<unsigned int> ignoreDuplicateArg("i", "ignoreDuplicateData",
						     "Reject apparently duplicated file names", false,
						     1, "Boolean indicator",
						     cmd);

    TCLAP::ValueArg<unsigned int> processNonMOSAICDataArg("o", "processOtherData",
							  "Process and convert other (i.e. non-MOSAIC) data",false,
							  1,"Boolean indicator", cmd);

    TCLAP::ValueArg<unsigned int> recursiveProcessingArg("c", "recursiveConversion",
							 "Convert input directory and all subdirectories",false,
							 1,"Boolean indicator",cmd);

    TCLAP::ValueArg<unsigned int> fMRIProcessingArg("f","fMRI","Convert 4D fMRI volumes into a set of Analyze files (for analysis with SPM)",false,0,"Boolean indicator",cmd);

    // Define a value argument and add it to the command line.
    TCLAP::MultiArg<std::string> restrictSeries("r",
						"restrictSeries",
						"Tag ID group|element to restrict file names",
						false,
						"Use Tag ID group|element to restrict file names - can be used more than once on the command line to add further restrictions.",
						cmd);

    TCLAP::ValueArg<std::string> mosaicOrientationArg("m","mosaicOrientation","Orientation of the MOSAIC images",false,"axial","image orientation",cmd);
    TCLAP::ValueArg<unsigned int> colorImageArg("l","color","Convert RGB DICOM images",false,0,"Boolean indicator",cmd);

    // Now parse the command line and assign variable values as required.
    cmd.parse( argc, argv );

    if (dicomArg.isSet()) { dicomDirectory = new std::string( dicomArg.getValue()); };
    if (outputPrefixArg.isSet()) { outputPrefix = new std::string( outputPrefixArg.getValue()); };
    if (outputSuffixArg.isSet()) { outputSuffix = new std::string( outputSuffixArg.getValue() ); };
    if (seriesIDArg.isSet()) { seriesID = new std::string( seriesIDArg.getValue()) ; };
    if (ignoreDuplicateArg.isSet()) {
      ignoreDuplicateData = ignoreDuplicateArg.getValue();
    };
    if (processNonMOSAICDataArg.isSet())
      processNonMOSAICData = processNonMOSAICDataArg.getValue();

    if (recursiveProcessingArg.isSet())
      recursiveProcessing = recursiveProcessingArg.getValue();

    if (fMRIProcessingArg.isSet())
      fMRIProcessing = fMRIProcessingArg.getValue();

    if (mosaicOrientationArg.isSet())
      mosaicOrientation = mosaicOrientationArg.getValue();

    if (colorImageArg.isSet())
      color = colorImageArg.getValue();

    if (restrictSeries.isSet()) {
      restrictionTags = restrictSeries.getValue();
    };
  } catch (TCLAP::ArgException &e)
    {
      std::cerr << "error: " << e.error() << " for arg " <<
	e.argId() << std::endl;
      exit(1);
    }


  // The purpose of this section is to identify an image from which we
  // can learn the appropriate output type for all of the data.
  //   For the time being, we are going to assume we can
  // learn this from a single file.

  typedef unsigned char PixelType;
  const unsigned int      ImageDimension = 3;

  typedef itk::Image< PixelType, ImageDimension >         ImageType;
  typedef itk::ImageSeriesReader< ImageType > SeriesReaderType;
  SeriesReaderType::Pointer seriesReader = SeriesReaderType::New();

  typedef   itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer volumeReader = ImageReaderType::New();

  typedef itk::GDCMImageIO ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  seriesReader->SetImageIO( dicomIO );
  volumeReader->SetImageIO( dicomIO );

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator =
    NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails(true);
  if (recursiveProcessing != 0)
    nameGenerator->RecursiveOn();
  else
    nameGenerator->RecursiveOff();

  for (std::vector<std::string>::iterator it = restrictionTags.begin(); it != restrictionTags.end(); it++) {
    std::cout << "Setting a series restriction with tag "
              << (*it) << std::endl;
    nameGenerator->AddSeriesRestriction( *it );
  }
  nameGenerator->SetDirectory( *dicomDirectory );

  try {
    typedef std::vector< std::string > SeriesIdContainer;
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    if (seriesUID.size() <= 0) {
      std::cout << "No DICOM series were identified in directory " <<
	(*dicomDirectory) << std::endl;
      return EXIT_FAILURE;
    } else {
      std::cout << "Found " << seriesUID.size() << " series." << std::endl;
    }

    // Iterate over all the series, seeing if we want to read each one.
    // If we want to read it, read it in and push it to disk using the
    // appropriate output type.
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

    std::string seriesIdentifier;

    // vectors to hold DWI data
    std::vector<DWIData> dwiData;

    while (seriesItr != seriesEnd) {
      if (seriesID->empty()) { // No optional series ID
        seriesIdentifier = seriesItr->c_str();
        seriesItr++;
      } else {
        seriesIdentifier = (*seriesID);
        seriesItr = seriesEnd;
      }
      // Now read in the series first file name and discover its type.
      typedef std::vector< std::string > FileNamesContainer;
      FileNamesContainer fileNames;
      nameGenerator->SetDirectory( *dicomDirectory );
      std::cout << "Finding names for " << seriesIdentifier << std::endl;
      fileNames = nameGenerator->GetFileNames( seriesIdentifier );
      itk::ImageIOBase *imageIO;
      // Might want to exclude duplicates here...
      if (fileNames.size() == 1) {
        try {
          volumeReader->SetFileName( fileNames[0] );
          volumeReader->GenerateOutputInformation();
          imageIO = volumeReader->GetImageIO();
        } catch (itk::ExceptionObject &excp)
	  {
	    std::cerr << "Failed to read from file " << fileNames[0] << std::endl;
	    std::cerr << "Exception error message is : " << excp << std::endl;
	    return EXIT_FAILURE;
	  }
      } else if (fileNames.size() > 1) {
        try {
          seriesReader->SetFileNames( fileNames );
          seriesReader->GenerateOutputInformation();
          imageIO = seriesReader->GetImageIO();
        } catch (itk::ExceptionObject &excp)
	  {
	    std::cerr << "Failed to read from files starting with "
		      << fileNames[0] << std::endl;
	    std::cerr << "Exception error message is : " << excp << std::endl;
	    return EXIT_FAILURE;
	  }
      } else {
        std::cerr << "Failed to find any files. " << std::endl;
        return EXIT_FAILURE;
      }

      crlConvertDICOMMOSAICBase::Pointer crlConvert;
      if (color)
	crlConvert = crlConvertDICOMMOSAIC<ImageDimension, itk::RGBPixel<unsigned char> >::New();
      else
	// Now attempt to read in the series and generating the most
	// appropriate closest matching output type.
	switch (imageIO->GetComponentType()) {
	case itk::ImageIOBase::DOUBLE : {
	  // Component type of data on file is DOUBLE, so use filter templated
	  // over double for this operation.
	  crlConvert = crlConvertDICOMMOSAIC<ImageDimension,
	      double>::New();
	};
	  break;
	case itk::ImageIOBase::FLOAT : {
	  // Component type of data on file is FLOAT, so use filter templated
	  // over float for this operation.
	  crlConvert = crlConvertDICOMMOSAIC<ImageDimension,
	      float>::New();
	};
	  break;
	case itk::ImageIOBase::UCHAR : {
	  crlConvert = crlConvertDICOMMOSAIC<ImageDimension,
	      unsigned char>::New();
	}; break;
	case itk::ImageIOBase::CHAR : {
	  crlConvert = crlConvertDICOMMOSAIC<ImageDimension,
	      signed char>::New();
	}; break;
	case itk::ImageIOBase::USHORT : {
	  crlConvert = crlConvertDICOMMOSAIC<ImageDimension,
	      unsigned short>::New();
	}; break;
	case itk::ImageIOBase::SHORT : {
	  crlConvert = crlConvertDICOMMOSAIC<ImageDimension,
	      signed short>::New();
	}; break;
	case itk::ImageIOBase::UINT : {
	  crlConvert = crlConvertDICOMMOSAIC<ImageDimension,
	      unsigned int>::New();
	}; break;
	case itk::ImageIOBase::INT : {
	  crlConvert = crlConvertDICOMMOSAIC<ImageDimension,
	      signed int>::New();
	}; break;
	case itk::ImageIOBase::ULONG : {
	  crlConvert = crlConvertDICOMMOSAIC<ImageDimension,
	      unsigned long>::New();
	}; break;
	case itk::ImageIOBase::LONG : {
	  crlConvert = crlConvertDICOMMOSAIC<ImageDimension,
	      signed long>::New();
	}; break;
	case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE :
	default : {
	  std::cerr << "The was a problem determining the file component type."
		    << std::endl;
	  return EXIT_FAILURE;
	}
	}
      crlConvert->SetInputDICOMDirectoryName( *dicomDirectory );
      crlConvert->SetOutputPrefix( *outputPrefix );
      crlConvert->SetOutputSuffix( *outputSuffix );
      crlConvert->SetSeriesID( seriesIdentifier );
      if (restrictionTags.size() != 0) {
        crlConvert->AddRestrictionTags( restrictionTags );
      }
      crlConvert->SetIgnoreDuplicateData( (ignoreDuplicateData != 0) );
      crlConvert->SetProcessNonMOSAICData((processNonMOSAICData != 0));
      crlConvert->SetRecursiveProcessing(recursiveProcessing != 0);
      crlConvert->SetFMRIProcessing(fMRIProcessing != 0);
      crlConvert->SetMOSAICOrientation(mosaicOrientation);


      int ret = crlConvert->Execute();
      if (ret != 0) {
        std::cerr << "crlConvertDICOMMOSAIC returned error code "
		  << std::endl;
        return EXIT_FAILURE;
      }

      if (crlConvert->dwiImage && (crlConvert->dwiDataVector[0].filename != ""))
	{
	  std::vector<DWIData>::iterator dit;
	  for (dit = crlConvert->dwiDataVector.begin(); dit != crlConvert->dwiDataVector.end(); dit++)
	    dwiData.push_back(*dit);
	}
    } // end of loop iterating over DICOM series in this directory


    // find number of unique series IDs in DWI data
    std::vector<DWIData>::iterator it;
    std::set<std::string> seriesNumbers;

    for (it = dwiData.begin(); it != dwiData.end(); it++)
      seriesNumbers.insert(it->series_number);
    std::cout << "Found " << seriesNumbers.size() << " DTI series." << std::endl;

    // create NHDR file for each series ID
    for (std::set<std::string>::iterator seriesIt = seriesNumbers.begin(); seriesIt != seriesNumbers.end(); seriesIt++)
      {
	std::string currentNumber = *seriesIt;
	std::vector<DWIData> currentdwiData;

	// copy relevant DWI data
	for (it = dwiData.begin(); it != dwiData.end(); it++)
	  if (it->series_number == currentNumber)
	    currentdwiData.push_back(*it);

	if (currentdwiData.size() < 2)
	  // if there is only one data set in a series, we are most likely dealing with fMRI data - skip it
	  {
	    std::cout << "Skipping series " << currentNumber << " because it contains only one data set." << std::endl;
	    continue;
	  }

	// write NHDR file
	unsigned int max_b_value = 0;

	// find highest b value
	for (it = currentdwiData.begin(); it != currentdwiData.end(); it++)
	  if (it->b_value > max_b_value) max_b_value = it->b_value;

	// adjust gradients for lower b values
	for (it = currentdwiData.begin(); it != currentdwiData.end(); it++)
	  if (it->b_value < max_b_value)
	    {
	      double factor = sqrt((double)it->b_value/(double)max_b_value);
	      double gradients[3];
	      std::stringstream ss;
	      ss.str(it->g_directions);
	      ss >> gradients[0];
	      ss >> gradients[1];
	      ss >> gradients[2];
	      gradients[0] = gradients[0] * factor;
	      gradients[1] = gradients[1] * factor;
	      gradients[2] = gradients[2] * factor;
	      std::stringstream ss2;
	      ss2 << gradients[0] << " " << gradients[1] << " " << gradients[2];
	      it->g_directions = ss2.str();
	    }

	if (currentdwiData.begin() != currentdwiData.end())
	  {
	    std::ofstream nhdrFile(std::string(*outputPrefix+"diffusion"+currentNumber+".nhdr").c_str());
	    if (!nhdrFile)
	      {
		std::cout << "Failed writing diffusion header file." << std::endl;
		return EXIT_FAILURE;
	      }
	    nhdrFile << "NRRD0005" << std::endl;
	    nhdrFile << "# Complete NRRD file format specification at:" << std::endl;
	    nhdrFile << "# http://teem.sourceforge.net/nrrd/format.html" << std::endl;
	    nhdrFile << "content: epi-dti-asset" << std::endl;
	    nhdrFile << "dimension: 4" << std::endl;
	    nhdrFile << "sizes: " << currentdwiData.begin()->sizes << " " << currentdwiData.size() << std::endl;
	    nhdrFile << "type: unsigned short" << std::endl;
	    nhdrFile << "encoding: raw" << std::endl;
	    nhdrFile << "endian: little" << std::endl;
	    nhdrFile << "byte skip: -1" << std::endl;
	    nhdrFile << "space: left-posterior-superior" << std::endl;
	    nhdrFile << "centers: cell cell cell none" << std::endl;
	    //nhdrFile << "space directions: " << currentdwiData.begin()->space_directions << " none" << std::endl;
	    // do not use the space directions of the first volume, they might be wrong (in case of several b=0 volumes)
	    nhdrFile << "space directions: " << currentdwiData[currentdwiData.size()-1].space_directions << " none" << std::endl;
	    nhdrFile << "space units: \" mm\" \"mm\" \"mm\"" << std::endl;
	    if ( currentdwiData[currentdwiData.size()-1].space_origin!="" )
	        nhdrFile << "space origin: " << currentdwiData[currentdwiData.size()-1].space_origin << std::endl;
	    nhdrFile << "kinds: space space space list" << std::endl;
	    nhdrFile << "measurement frame: (1,0,0) (0,1,0) (0,0,1)" << std::endl;
	    nhdrFile << "modality:=DWMRI" << std::endl;
	    nhdrFile << "DWMRI_b-value:=" << max_b_value << std::endl;
	    unsigned int g = 0;
	    for (it = currentdwiData.begin(); it != currentdwiData.end(); it++)
	      {
		nhdrFile << "DWMRI_gradient_" << std::right << std::setw(4) << std::setfill('0') << g << ":= " << currentdwiData.at(g).g_directions << std::endl;
		g++;
	      }
	    nhdrFile << "data file: LIST" << std::endl;
	    for (it = currentdwiData.begin(); it != currentdwiData.end(); it++)
	      nhdrFile << it->filename << std::endl;

	    nhdrFile.close();
	  }
      } // end of iteration over series IDs

  } catch (itk::ExceptionObject &ex) {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}
