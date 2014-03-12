
#ifndef _CRL_CONVERTDICOMMOSAIC
#define _CRL_CONVERTDICOMMOSAIC 1

// This reads a volume represented by a DICOM series and writes it out
// in a new format, based on the file extension.

#include "crlConvertDICOMMOSAIC.h"

#include <iostream>
#include <fstream>
#include <sys/stat.h>

#include <itkExtractImageFilter.h>
#include <itkTileImageFilter.h>
#include <itkOrientImageFilter.h>
#include <itkDiffusionTensor3D.h>

#include "gdcmDictSet.h"
#include "gdcmDict.h"
#include "gdcmFile.h"
#include "gdcmDocEntry.h"
#include "gdcmBinEntry.h"
#include "gdcmDictEntry.h"
#include "gdcmGlobal.h"

#ifdef WIN32
#include <direct.h>
#endif

std::string basename(std::string path)
{
  /* Find the last path element: */
  size_t pos = path.find_last_of('/');
  
  /* Make sure we do the right thing for, /file/name/ */
  if (pos == std::string::npos) return path;
  if (pos == (path.size() - 1)) {
    return basename(path.substr(0, path.size()-1));
  }
  
  return path.substr(pos+1);
  
}

int sign(double v)
{
  return ((v > 0) - (v < 0));
}


template <const unsigned int TImageDimension, class TComponentType >
int crlConvertDICOMMOSAIC< TImageDimension, TComponentType>
::Execute()
{
  // The image type is largely defined by the template parameeters;
  typedef itk::Image< TComponentType, TImageDimension >         ImageType;
  typedef itk::Image< TComponentType, TImageDimension-1 > ImageTypeSlice;
  
  typedef itk::ImageSeriesReader< ImageType >        SeriesReaderType;
  typename SeriesReaderType::Pointer reader = SeriesReaderType::New();
  
  typedef itk::ImageFileReader< ImageType > VolumeReaderType;
  typename VolumeReaderType::Pointer volumereader = VolumeReaderType::New();
  
  typedef itk::GDCMImageIO       ImageIOType;
  typename ImageIOType::Pointer dicomIO = ImageIOType::New();
  
  reader->SetImageIO( dicomIO );
  volumereader->SetImageIO( dicomIO );
  dicomIO->SetLoadPrivateTags(true);
  
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  typename NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  
  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->SetLoadPrivateTags(true);
  if (m_recursiveProcessing)
    nameGenerator->RecursiveOn();
  else
    nameGenerator->RecursiveOff();
  if (!m_RestrictionTag.empty()) {
    nameGenerator->AddSeriesRestriction( m_RestrictionTag );
  }
  for (std::vector<std::string>::iterator it = m_RestrictionTags.begin(); it != m_RestrictionTags.end(); it++)
    {
      nameGenerator->AddSeriesRestriction(*it);
    }
  
  nameGenerator->SetDirectory( m_dicomDir );
  
  std::string outputFilePrefix( m_outputPrefix );
  
  try
    {
      //std::cout << std::endl << "The directory: " << std::endl;
      //std::cout << std::endl << m_dicomDir;
      //std::cout << " contains the following DICOM Series: ";
      //std::cout << std::endl << std::endl;
      
      typedef std::vector< std::string >    SeriesIdContainer;
      
      const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
      
      if (seriesUID.size() <= 0) {
	std::cout << "No series were identified." << std::endl;
	return EXIT_FAILURE;
      }
      
      SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
      SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
      /*while( seriesItr != seriesEnd )
	{
	std::cout << seriesItr->c_str() << std::endl;
	seriesItr++;
	}*/
      
      seriesItr = seriesUID.begin();
      seriesEnd = seriesUID.end();
      while( seriesItr != seriesEnd ) {
	
	std::string seriesIdentifier;
	
	if( !m_SeriesID.empty() ) // If no optional series identifier
	  {
	    seriesIdentifier = m_SeriesID;
	    seriesItr = seriesEnd;
	  }
	else
	  {
	    seriesIdentifier = seriesItr->c_str();
	    seriesItr++;
	  }
	
	std::cout << std::endl ;
	std::cout << "Now reading series: " << std::endl << std::endl;
	std::cout << "   " << seriesIdentifier << std::endl;
	std::cout << std::endl ;
	
	// We pass the series identifier to the name generator and ask for all the
	// filenames associated to that series. This list is returned in a container of
	// strings by the GetFileNames() method.
	
	typedef std::vector< std::string >   FileNamesContainer;
	FileNamesContainer fileNames;
	
	fileNames = nameGenerator->GetFileNames( seriesIdentifier );
	
	// if the data is pushed to the scanner twice, it can get the same data
	// but stored with a leading 00x_ file name.
	// We want to exclude such files from further consideration.
	
	if (m_IgnoreDuplicateData) {
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
	
	// Declare the writer here so we can adjust the reader input as appropriate
	typedef itk::ImageFileWriter< ImageType > WriterType;
	typename WriterType::Pointer writer = WriterType::New();
	
	bool writeThisFile = m_processNonMOSAICData;
	
	// FIX FIX - CHECK FOR 2 NEX DT-MRI SCAN WITH ONE VOLUME SLICE
	if (fileNames.size() == 1) {
	  volumereader->SetFileName( fileNames[0].c_str() );
	  try {
	    volumereader->Update();
	  } catch (itk::ExceptionObject &ex) {
	    std::cout << ex << std::endl;
	    return 1;
	  }
	  writer->SetInput( volumereader->GetOutput() );
	  
	} else if (fileNames.size() > 1) {
	  reader->SetFileNames( fileNames );
	  try {
	    reader->Update();
	  } catch (itk::ExceptionObject &ex) {
	    std::cout << ex << std::endl;
	    return 1;
	  }
	  writer->SetInput( reader->GetOutput() );
	  
	} else {
	  std::cerr << "Didn't find enough file names." << std::endl;
	}
	
	// Let's interrogate the series information
	typedef itk::MetaDataDictionary   DictionaryType;
	const  DictionaryType & dictionary = dicomIO->GetMetaDataDictionary();
	std::vector< std::string > tagNames(5);
	std::vector< std::string > tagNamesDescriptions(5);
	std::vector< std::string> tagValues(5);
	
	tagNames[0] = "0008|103e";
	tagNamesDescriptions[0] = "Series Description ";
	tagNames[1] = "0018|0024";
	tagNamesDescriptions[1] = "Sequence Name ";
	tagNames[2] = "0019|100a";
	tagNamesDescriptions[2] = "Number of Images in Mosaic ";
	tagNames[3] = "0020|0011";
	tagNamesDescriptions[3] = "Series Number";
	tagNames[4] = "0008|0031";
	tagNamesDescriptions[4] = "Series Time";
	
	
	DictionaryType::ConstIterator tagItr;// = dictionary.Find(seriesDescriptionTag);
	DictionaryType::ConstIterator end = dictionary.End();
	typedef itk::MetaDataObject< std::string > MetaDataStringType;
	MetaDataStringType::ConstPointer entryvalue;
	
	for (unsigned int i = 0; i < tagNames.size(); i++) {
	  tagItr = dictionary.Find( tagNames[i] );
	  if( tagItr == end )
	    {
	      std::cerr << "   " << tagNamesDescriptions[i] ;
	      std::cerr << " not found in the DICOM header" << std::endl;
	      tagValues[i] = tagNamesDescriptions[i];
	      tagValues[i] += "Not_Found";
	    } else {
	    // Since the entry may or may not be of string type we must again use a
	    // dynamic_cast in order to attempt to convert it to a string dictionary
	    // entry. If the conversion is successful, then print out its content.
	    
	    entryvalue =
	      dynamic_cast<const MetaDataStringType *>( tagItr->second.GetPointer());
	    if ( entryvalue )
	      {
		tagValues[i] = entryvalue->GetMetaDataObjectValue();
		std::cout << "   " << tagNamesDescriptions[i] << " (" << tagNames[i] <<  ") ";
		std::cout << " : " << tagValues[i] << std::endl;
	      }
	    else
	      // if the conversion is not successful, print out an annoyed grunt
	      std::cout << "GRUNT!" << std::endl;
	  }
	}
	
	// Convert spaces in the tag value to underscore to make me happy
	// Also, convert path separator characters to underscore.
	for (unsigned int i = 0; i < tagValues.size(); i++) {
	  for (unsigned int j = 0; j < tagValues[i].size(); j++) {
	    if ( isspace(tagValues[i][j]) ||
		 (tagValues[i][j] == '/') ||
		 (tagValues[i][j] == '\\') ||
		 (tagValues[i][j] == ';') ||
		 ( (i==0 || i==1) && (tagValues[i][j] == ',') ) ||
		 (tagValues[i][j] == ':') ||
		 (tagValues[i][j] == '*') ||
		 (tagValues[i][j] == '"') ||
		 (tagValues[i][j] == '\0') ||
		 (tagValues[i][j] == '\'')
		 ) {
	      tagValues[i].replace(j,1,"_");
	    }
	  }
	}
	
	std::string outName = m_outputPrefix;
        outName += tagValues[0];
        outName += tagValues[1];
        outName += seriesIdentifier;
        outName += m_outputSuffix;
	
	writer->SetFileName( outName.c_str() );
	std::cout  << "   Writing the image as " << std::endl << std::endl;
	std::cout  << "   " << outName << std::endl << std::endl;
	
	//---------------------------------------------------
	// If 0x0019 was not found, try to read the CSA header
	//---------------------------------------------------
	bool readFromCSA = false;
	std::string CSAHeader;
	
	if ((tagValues[2].find("Not_Found",0) != std::string::npos)||(atoi(tagValues[2].c_str()) == 0))
	  {
	    //------------------------------------------------
	    // Read the Siemens CSA header in 'CSAHeader' to get
	    // the number of mosaics if '0019|100a' is not in the dicom
	    // The whole CSA header is loaded in CSAHeader.
	    //
	    // The CSA header (0029|1010) is a binary tag.
	    // It seems that tagItr=dictionary.Find("0029|1010")
	    // doesn't work (because it is binary?) so we have
	    // to read it with gdcm
	    //------------------------------------------------
	    
	    // Load the file with gdcm
	    gdcm::File *header0 = new gdcm::File;
	    gdcm::BinEntry* binEntry;
	    header0->SetMaxSizeLoadEntry(163840);
	    header0->SetFileName( fileNames[0] );
	    header0->SetLoadMode( gdcm::LD_ALL );
	    header0->Load();
	    
	    // Parse the tags to find '0029|1010'
	    gdcm::DocEntry* docEntry = header0->GetFirstEntry();
	    while(docEntry)
	      {
		if ( docEntry->GetKey() == "0029|1010"  )
		  {
		    std::cout << "   Read the CSA header" << std::endl;
		    
		    // Get the binary entry
		    binEntry = dynamic_cast<gdcm::BinEntry*> ( docEntry );
		    int binLength = binEntry->GetFullLength();
		    
		    CSAHeader.resize( binLength );
		    uint8_t * p = binEntry->GetBinArea();
		    
		    for (int k = 0; k < binLength; k++)
		      CSAHeader[k] = *(p+k);
		    break;
		  }
		docEntry = header0->GetNextEntry();
	      }
	    
	    // If we found 0029|1010... fix the value of tagValues[2]
	    // (the b-value/grad dir will be extracted after, if "0019|1027"
	    // is not found)
	    if ( CSAHeader.size()!=0 )
	      {
		size_t pos = CSAHeader.find( "NumberOfImagesInMosaic" );
		if ( pos != std::string::npos)
		  {
		    int itemLength = *(CSAHeader.c_str()+pos+72+16);
		    tagValues[2] = CSAHeader.substr( pos+84+16, itemLength );
		    readFromCSA = true;
		    // note: +16 = size of the CSA2 header
		    //       +84 = +64+20 = end of that tag = offset of the value
		    // cf http://nipy.sourceforge.net/dipy/theory/siemens_csa.html
		  }
		else
		  std::cout << "   !! WARNING: CSA header found but not NumberOfImagesInMosaic" << std::endl;
	      }
	  }
	
	if ((tagValues[2].find("Not_Found",0) == std::string::npos) && (atoi(tagValues[2].c_str()) != 0))
	  {
	    std::cout << "We are converting a MOSAIC image." << std::endl;
	    // We are dealing with a MOSAIC image
	    writeThisFile = true;
	    
	    unsigned int numImages = atoi(tagValues[2].c_str());
	    unsigned int rowImages = (unsigned int)ceil(sqrt((double)numImages));
	    
	    typedef itk::ExtractImageFilter<ImageType,ImageTypeSlice> ExtractFilterType;
	    typename ExtractFilterType::Pointer ExtractFilter = ExtractFilterType::New();
	    
	    typedef itk::TileImageFilter<ImageTypeSlice,ImageType> TileFilterType;
	    typename TileFilterType::Pointer TileFilter = TileFilterType::New();
	    
	    typename TileFilterType::LayoutArrayType layout;
	    layout[0] = 1;
	    layout[1] = 1;
	    layout[2] = 0;
	    TileFilter->SetLayout(layout);
	    
	    typename ImageType::ConstPointer mosaicImage;
	    mosaicImage  = writer->GetInput();
	    
	    writer->SetInput(TileFilter->GetOutput());
	    this->dwiImage=false;
	    
	    
	    // Retrieve B matrix from DICOM header
	    tagItr=dictionary.Find("0019|1027");
	    if (tagItr!=end)
	      {
		typedef itk::MetaDataObject<std::string> StringMatrixType;
		StringMatrixType::ConstPointer matrixObject;
		
		matrixObject = dynamic_cast<const StringMatrixType*>(tagItr->second.GetPointer());
		if (matrixObject)
		  {
		    double *matrix = (double*)(matrixObject->GetMetaDataObjectValue().c_str());
		    std::cout << "   B matrix: bxx=" << matrix[0] << " bxy=" << matrix[1] << " bxz=" << matrix[2] << " byy=" << matrix[3] << " byz=" << matrix[4] << " bzz=" << matrix[5] << std::endl;
		    
		    tagItr = dictionary.Find("0019|100c");
		    matrixObject = dynamic_cast<const StringMatrixType*>(tagItr->second.GetPointer());
		    std::cout << "   B value: " << matrixObject->GetMetaDataObjectValue() << std::endl;
		    
		    // store some information
		    this->dwiImage = true;
		    dwiData.filename = basename(outName);
		    dwiData.b_value = atoi(matrixObject->GetMetaDataObjectValue().c_str());
		    dwiData.series_number = tagValues[4]+tagValues[3];
		    
		    itk::DiffusionTensor3D<double> t;
		    
		    for (int i = 0; i < 6; i++)
		      t[i] = matrix[i];
		    
		    itk::DiffusionTensor3D<double>::EigenValuesArrayType eigenValues;
		    itk::DiffusionTensor3D<double>::EigenVectorsMatrixType eigenVectors;
		    t.ComputeEigenAnalysis(eigenValues,eigenVectors);
		    // eigenVectors[2] contains the principal eigenvector;
		    
		    //double gradientDirections[3];
		    //int bsign[3];
		    
		    //bsign[0] = sign(sign(matrix[0])+0.01);
		    //bsign[1] = sign(sign(matrix[1])+0.01);
		    //bsign[2] = sign(sign(matrix[2])+0.01);
		    
		    //gradientDirections[0] = bsign[0] * sqrt(matrix[0] / dwiData.b_value);
		    //gradientDirections[1] = bsign[1] * sqrt(matrix[3] / dwiData.b_value);
		    //gradientDirections[2] = bsign[2] * sqrt(matrix[5] / dwiData.b_value);
		    
		    std::stringstream ss;
		    //ss << gradientDirections[0] << " " << gradientDirections[1] << " " << gradientDirections[2];
		    ss << eigenVectors[2][0] << " " << eigenVectors[2][1] << " " << eigenVectors[2][2];
		    dwiData.g_directions = ss.str();
		    std::cout << "   Computed gradient directions: " << dwiData.g_directions << std::endl;
		  } // if (matrixObject)
		else
		  std::cout << "   B-Matrix found, but DICOM attribute format is unknown." << std::endl;
	      }
	    
	    //----------------------------------------------
	    // If "0019|1027" was not found, but we have the CSA
	    // header, try to extract the B value and the gradient vectors from the CSA!
	    //----------------------------------------------
	    else if ( readFromCSA )
	      {
		// store some informations
		dwiData.filename = basename(outName);
		dwiData.series_number = tagValues[4]+tagValues[3];
		this->dwiImage = true;
		
		// Extract the b-value
		size_t pos = CSAHeader.find( "B_value" );
		if ( pos != std::string::npos)
		  {
		    int itemLength = *(CSAHeader.c_str()+pos+72+16);
		    dwiData.b_value = atoi(CSAHeader.substr( pos+84+16, itemLength ).c_str());
		    std::cout << "   B value (from CSA): " << dwiData.b_value << std::endl;
		  }
		else
		  {
		    std::cout << "   !! WARNING: B not found in the CSA. Set to 0." << std::endl;
		    dwiData.b_value = 0;
		  }
		
		if ( dwiData.b_value==0 )
		  {
		    dwiData.g_directions="0 0 0";
		  }
		else
		  {
		    // Extract the diffusion gradients
		    pos = CSAHeader.find( "DiffusionGradientDirection" );
		    if ( pos != std::string::npos)
		      {
			int itemLength = *(CSAHeader.c_str()+pos+72+16);
			std::string v1 = CSAHeader.substr( pos+84+16, itemLength-1 );
			std::string v2 = CSAHeader.substr( pos+84+16+28, itemLength-1 );
			std::string v3 = CSAHeader.substr( pos+84+16+28+28, itemLength-1 );
			
			// !!!! not sure at all about the +28  and +28+28 !!!!!!
			
			std::stringstream ss;
			ss << v1.c_str() << " " << v2.c_str() << " " << v3.c_str();
			dwiData.g_directions = ss.str();
			std::cout << "   Gradient directions (from CSA): " << dwiData.g_directions << std::endl;
		      }
		    else
		      {
			std::cout << "   !! WARNING: Gradient directions not found in the CSA. Set to 0 0 0." << std::endl;
			dwiData.g_directions="0 0 0";
		      }
		  }
	      }
	    
	    
	    
	    // Retrieve gradient directions from DICOM header
	    /*	    tagItr=dictionary.Find("0019|100e");
	      if (tagItr!=end)
	      {
	      typedef itk::MetaDataObject< std::string > StringDirectionType;
	      StringDirectionType::ConstPointer gradientDirectionsObject;
	      
	      gradientDirectionsObject = dynamic_cast< const StringDirectionType* >(tagItr->second.GetPointer());
	      if (gradientDirectionsObject)
	      {
	      const char* gradientDirectionsString = gradientDirectionsObject->GetMetaDataObjectValue().c_str();
	      double* gradientDirections = (double*)gradientDirectionsString;
	      std::cout << "Gradient directions: " << gradientDirections[0] << " " << gradientDirections[1] << " " << gradientDirections[2] << std::endl;
	      tagItr=dictionary.Find("0019|100c");
	      gradientDirectionsObject = dynamic_cast< const StringDirectionType* >(tagItr->second.GetPointer());
	      std::cout << "B value: " << gradientDirectionsObject->GetMetaDataObjectValue() << std::endl;
	      this->dwiImage=true;
	      
	      // store filename, b value, gradient directions, series id
	      dwiData.filename = basename(outName);
	      dwiData.b_value = atoi(gradientDirectionsObject->GetMetaDataObjectValue().c_str());
	      std::stringstream ss;
	      ss << gradientDirections[0] << " " << gradientDirections[1] << " " << gradientDirections[2];
	      dwiData.g_directions = ss.str();
	      dwiData.series_number = tagValues[4]+tagValues[3];
	      }
	      else
	      std::cout << "Gradient directions found, but DICOM attribute format is unknown." << std::endl;
	      }*/
	    else
	      {
		// if there is no B matrix, can we assume that this is a b=0 image?
		// let's give it a try!
		std::cout << "   No B matrix found, but the image is MOSAIC --> b = 0" << std::endl;
		this->dwiImage=true;
		dwiData.filename = basename(outName);
		dwiData.b_value = 0;
		dwiData.g_directions="0 0 0";
		dwiData.series_number = tagValues[4]+tagValues[3];
	      }
	    
	    // Compute size of each tile of the MOSAIC image
	    typename ImageType::SizeType inputSize;
	    inputSize = mosaicImage->GetLargestPossibleRegion().GetSize();
	    unsigned int tileX = inputSize[0] / rowImages;
	    unsigned int tileY = inputSize[1] / rowImages;
	    bool fMRI = false;
	    
	    if (inputSize[2] > 1 && m_fMRIProcessing)
	      {
		// we have a 4D fMRI volume here
		std::cout << "   4D fMRI volume detected" << std::endl;
		fMRI = true;
	      }
	    
	    std::stringstream ss;
	    if (m_mosaicOrientation == "coronal")
	      ss << tileX <<" " << numImages << " " << tileY;
	    else if (m_mosaicOrientation == "sagittal")
	      ss << numImages << " " << tileX << " " << tileY;
	    else // axial
	      ss << tileX << " " << tileY << " " << numImages;
	    
	    dwiData.sizes = ss.str();
	    
	    // get slice thickness - use spacing between slices instead of slice thickness (0018|0050)
	    tagItr = dictionary.Find("0018|0088");
	    
	    entryvalue = dynamic_cast<const MetaDataStringType *>( tagItr->second.GetPointer());
	    if (entryvalue)
	      {
		double newSpacing[3];
		newSpacing[0] = mosaicImage->GetSpacing()[0];
		newSpacing[1] = mosaicImage->GetSpacing()[1];
		newSpacing[2] = atof(entryvalue->GetMetaDataObjectValue().c_str());
		reader->GetOutput()->SetSpacing(newSpacing);
	      }
	    
	    ss.str("");
	    //ss << "(" << mosaicImage->GetSpacing()[0] << ", 0.0, 0.0) ";
	    //ss << "(0.0, " << mosaicImage->GetSpacing()[1] << ", 0.0) ";
	    //ss << "(0.0, 0.0, " << mosaicImage->GetSpacing()[2] << ")";
	    ss << "(" << mosaicImage->GetDirection()[0][0]*mosaicImage->GetSpacing()[0] << ", " << mosaicImage->GetDirection()[1][0]*mosaicImage->GetSpacing()[0] << ", " << mosaicImage->GetDirection()[2][0]*mosaicImage->GetSpacing()[0] << ") ";
	    ss << "(" << mosaicImage->GetDirection()[0][1]*mosaicImage->GetSpacing()[1] << ", " << mosaicImage->GetDirection()[1][1]*mosaicImage->GetSpacing()[1] << ", " << mosaicImage->GetDirection()[2][1]*mosaicImage->GetSpacing()[1] << ") ";
	    ss << "(" << mosaicImage->GetDirection()[0][2]*mosaicImage->GetSpacing()[2] << ", " << mosaicImage->GetDirection()[1][2]*mosaicImage->GetSpacing()[2] << ", " << mosaicImage->GetDirection()[2][2]*mosaicImage->GetSpacing()[2] << ")";
	    dwiData.space_directions = ss.str();
	    
	    
	    // read the origin
	    ss.str("");
	    tagItr=dictionary.Find("0020|0032");
	    if (tagItr!=end)
	      {
		typedef itk::MetaDataObject<std::string> StringVectorType;
		StringVectorType::ConstPointer vectObject;
		
		vectObject = dynamic_cast<const StringVectorType*>(tagItr->second.GetPointer());
		if (vectObject)
		  {
		    // replace the '\' by a ',' to make a valid vector
		    std::string origin = (vectObject->GetMetaDataObjectValue().c_str());
		    std::replace( origin.begin(), origin.end(), '\\', ',' );
		    
		    ss << "(" << origin << ")";
		    std::cout << "   Origin: "<<origin<<std::endl;
		  }
	      }
	    dwiData.space_origin = ss.str();
	    
	    typename ImageType::RegionType ExtractRegion;
	    ExtractRegion.SetSize(0,tileX);
	    ExtractRegion.SetSize(1,tileY);
	    
	    ExtractFilter->SetInput(mosaicImage);
	    typename ImageTypeSlice::Pointer currentSlice;
	    
	    // Now extract all mosaic tiles and put them together into a new volume
	    for (unsigned int k = 0; k < inputSize[2]; k++)
	      {
		for (unsigned int j = 0;j < rowImages*tileY-1; j+=tileY)
		  for (unsigned int i = 0;i < rowImages*tileX-1; i+=tileX)
		    {
		      // only process this tile if there actually is an image
		      // (in the last row there may be empty spots)
		      if ((i/tileX + j/tileY*rowImages) < numImages)
			{
			  ExtractRegion.SetIndex(0,i);
			  ExtractRegion.SetIndex(1,j);
			  ExtractRegion.SetIndex(2,k);
			  ExtractFilter->SetExtractionRegion(ExtractRegion);
			  ExtractFilter->UpdateLargestPossibleRegion();
			  currentSlice = ExtractFilter->GetOutput();
			  currentSlice->DisconnectPipeline();
			  TileFilter->PushBackInput(currentSlice);
			  currentSlice = NULL;
			}
		    }
		
		TileFilter->UpdateLargestPossibleRegion();
		
		TileFilter->GetOutput()->SetOrigin(mosaicImage->GetOrigin());
		TileFilter->GetOutput()->SetSpacing(mosaicImage->GetSpacing());
		TileFilter->GetOutput()->SetDirection(mosaicImage->GetDirection());
		
		if (fMRI)
		  {
		    // create NIFTI format files for SPM
		    std::string name = outName.substr(0,outName.rfind(".")-1);
		    typename WriterType::Pointer fMRIWriter = WriterType::New();
		    typename itk::OrientImageFilter<ImageType,ImageType>::Pointer orienter = itk::OrientImageFilter<ImageType,ImageType>::New();
		    orienter->SetInput(TileFilter->GetOutput());
		    orienter->UseImageDirectionOn();
		    orienter->SetDesiredCoordinateOrientationToAxial();
		    orienter->Update();
		    std::stringstream fs;
		    fs << name << "/" << "fmri_" << std::setw(3) << std::setfill('0') << k+1 << ".nii";
#ifdef WIN32
		    _mkdir(name.c_str());
#else
		    mkdir(name.c_str(),0777);
#endif
		    fMRIWriter->SetFileName(fs.str().c_str());
		    fMRIWriter->SetInput(orienter->GetOutput());
		    
		    // change direction
		    typename ImageType::DirectionType direction;
		    direction.SetIdentity();
		    orienter->GetOutput()->SetDirection(direction);
		    
		    // change slice thickness - use spacing between slices for a correct value
		    tagItr = dictionary.Find("0018|0088");
		    
		    entryvalue = dynamic_cast<const MetaDataStringType *>( tagItr->second.GetPointer());
		    if (entryvalue)
		      {
			double newSpacing[3];
			newSpacing[0] = orienter->GetOutput()->GetSpacing()[0];
			newSpacing[1] = orienter->GetOutput()->GetSpacing()[1];
			newSpacing[2] = atof(entryvalue->GetMetaDataObjectValue().c_str());
			orienter->GetOutput()->SetSpacing(newSpacing);
		      }
		    
		    fMRIWriter->Write();
		    
		    while (TileFilter->GetNumberOfInputs() > 0)
		      TileFilter->PopBackInput();
		    writeThisFile = false;
		    this->dwiImage = false;
		  }
		else if (inputSize[2] > 1)
		  {
		    // no fMRI, but a 4D DWI dataset
		    /*typename itk::OrientImageFilter<ImageType,ImageType>::Pointer orienter = itk::OrientImageFilter<ImageType,ImageType>::New();
		      orienter->SetInput(TileFilter->GetOutput());
		      orienter->UseImageDirectionOn();
		      orienter->SetDesiredCoordinateOrientationToAxial();
		      orienter->Update();*/
		    
		    std::string name = outName.substr(0,outName.rfind(".")-1);
		    typename WriterType::Pointer dwiWriter = WriterType::New();
		    //writer->SetInput(orienter->GetOutput());
		    std::stringstream fs;
		    fs << name << k << ".nrrd";
		    std::cout << "Re-setting the file name to be : " << fs.str().c_str() << std::endl;
		    writer->SetFileName(fs.str().c_str());
		    dwiData.filename = basename(fs.str());
		    writer->UseCompressionOff();
		    
		    // change slice thickness - use spacing between slices for a correct value
		    tagItr = dictionary.Find("0018|0088");
		    
		    entryvalue = dynamic_cast<const MetaDataStringType *>( tagItr->second.GetPointer());
		    if (entryvalue)
		      {
			double newSpacing[3];
			newSpacing[0] = TileFilter->GetOutput()->GetSpacing()[0];
			newSpacing[1] = TileFilter->GetOutput()->GetSpacing()[1];
			newSpacing[2] = atof(entryvalue->GetMetaDataObjectValue().c_str());
			TileFilter->GetOutput()->SetSpacing(newSpacing);
		      }
		    
		    // we need to create a copy of the output of the tile filter, so that our spacing and origin
		    // changed don't get changed by the rest of the ITK pipeline
		    typename ImageType::Pointer volume = TileFilter->GetOutput();
		    volume->DisconnectPipeline();
		    writer->SetInput(volume);
		    writer->Write();
		    while (TileFilter->GetNumberOfInputs() > 0)
		      TileFilter->PopBackInput();
		    
		    dwiDataVector.push_back(dwiData);
		    writeThisFile = false;
		  }
	      }
	  }
	else
	  {
	    // not a MOSAIC image, but it could still be DTI
	    // let's check for gradient directions
	    tagItr=dictionary.Find("0019|100e");
	    if (tagItr!=end)
	      {
		typedef itk::MetaDataObject<std::string> StringVectorType;
		StringVectorType::ConstPointer vectObject;
		
		vectObject = dynamic_cast<const StringVectorType*>(tagItr->second.GetPointer());
		if (vectObject)
		  {
		    double *vect = (double*)(vectObject->GetMetaDataObjectValue().c_str());
		    std::stringstream ss;
		    ss << vect[0] << " " << vect[1] << " " << vect[2];
		    std::cout << "   Gradient directions: " << ss.str() << std::endl;
		    dwiData.g_directions=ss.str();
		    tagItr = dictionary.Find("0019|100c");
		    vectObject = dynamic_cast<const StringVectorType*>(tagItr->second.GetPointer());
		    std::cout << "   B value: " << vectObject->GetMetaDataObjectValue() << std::endl;
		    this->dwiImage = true;
		    dwiData.filename = basename(outName);
		    dwiData.b_value = atoi(vectObject->GetMetaDataObjectValue().c_str());
		    dwiData.series_number = tagValues[4]+tagValues[3];
		    ss.str("");
		    typename ImageType::ConstPointer image = writer->GetInput();
		    ss << "(" << image->GetDirection()[0][0]*image->GetSpacing()[0] << ", " << image->GetDirection()[1][0]*image->GetSpacing()[0] << ", " << image->GetDirection()[2][0]*image->GetSpacing()[0] << ") ";
		    ss << "(" << image->GetDirection()[0][1]*image->GetSpacing()[1] << ", " << image->GetDirection()[1][1]*image->GetSpacing()[1] << ", " << image->GetDirection()[2][1]*image->GetSpacing()[1] << ") ";
		    ss << "(" << image->GetDirection()[0][2]*image->GetSpacing()[2] << ", " << image->GetDirection()[1][2]*image->GetSpacing()[2] << ", " << image->GetDirection()[2][2]*image->GetSpacing()[2] << ")";
		    dwiData.space_directions = ss.str();
		    ss.str("");
		    ss << image->GetLargestPossibleRegion().GetSize()[0] << " " << image->GetLargestPossibleRegion().GetSize()[1] << " " << image->GetLargestPossibleRegion().GetSize()[2];
		    dwiData.sizes = ss.str();
		    writeThisFile = true;
		    
		    // Read the origin
		    ss.str("");
		    tagItr=dictionary.Find("0020|0032");
		    if (tagItr!=end)
		      {
			vectObject = dynamic_cast<const StringVectorType*>(tagItr->second.GetPointer());
			if (vectObject)
			  {
			    // replace the '\' by a ',' to make a valid vector
			    std::string origin = (vectObject->GetMetaDataObjectValue().c_str());
			    std::replace( origin.begin(), origin.end(), '\\', ',' );
			    ss << "(" << origin << ")";
			  }
		      }
		    dwiData.space_origin = ss.str();
		    
		  }
	      }
	  }
	
	try
	  {
	    if (!this->dwiImage) {
	      writer->UseCompressionOff();
            }
	    else {
	      writer->UseCompressionOff();
            }
	    if (writeThisFile)
	      {
		writer->Update();
		dwiDataVector.push_back(dwiData);
	      }
	    else {
	      std::cout << "   Not writing this file." << std::endl;
            }
	  }
	catch (itk::ExceptionObject &ex)
	  {
	    std::cout << ex << std::endl;
            std::cout << "Exception in first exception handler : writing the file failed." << std::endl;
	    return EXIT_FAILURE;
	  }
      }
    }
  catch (itk::ExceptionObject &ex)
    {
      std::cout << ex << std::endl;
      std::cout << "Exception in second exception handler : writing the file failed." << std::endl;
      return EXIT_FAILURE;
    }
  
  std::cout << std::endl << std::endl;
  return EXIT_SUCCESS;
  
};

#endif
