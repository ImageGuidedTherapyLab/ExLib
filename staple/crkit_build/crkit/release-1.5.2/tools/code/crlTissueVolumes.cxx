#include <iostream>
#include "tclap/CmdLine.h"
#include "configuration.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"

typedef float PixelType;
typedef itk::Image<PixelType,3> ImageType;
typedef itk::ImageFileReader<ImageType> ImageReaderType;
typedef itk::ImageRegionConstIterator<ImageType> ImageConstIterator;


int main( int argc, char* argv[] )
{
  const unsigned int maxLabel = 4096;       // 0..maxLabel.  Note that for -l we assume digits only 0..4096
  

  TCLAP::CmdLine cmd("Computational Radiology Laboratory", ' ', 
        CRKIT_VERSION_STRING  
  );
  TCLAP::ValueArg<std::string> optLabelsArg( "l", "labels", 
    "A string of integer labels to restrict our interest to", false, "", "string",cmd);
  TCLAP::UnlabeledMultiArg<std::string> optFilesArg("inputs","input file name(s)",true,"string",cmd) ;
  try
    {
    cmd.parse(argc,argv);
    }
  catch (TCLAP::ArgException& e)
    {
    std::cerr << "command line error: " << e.error() << " for arg " <<
      e.argId() << std::endl;
    return 1;
    }

  ///////// 
  ///// we have a bool vector from 0..maxLabel telling whether to use a
  ///// label.  We parse a command line parameter here '-l' which takes
  ///// string of integer label names.
  /////////
  std::vector<bool> labels(maxLabel+1,true);
  if ( optLabelsArg.isSet() )
    {
    labels.assign(maxLabel+1,false);
    std::string ss = optLabelsArg.getValue();
    const char *s = optLabelsArg.getValue().c_str();
    char *end;
    unsigned int l = 0;
    while ( (static_cast<unsigned char>(*s) != '\0') &&
            isspace(static_cast<unsigned char>(*s)) ) s++;
    do {
      l = static_cast<unsigned int>( strtol(s, &end, 10));
      if (s == end) break; // we are done converting
      
      if ( l > maxLabel || l < 0 )
        {
          std::cerr << "bad label " << l << std::endl;
        return 1;
        }
      labels[l] = true;
      s = end;  // if not done, then start again where we left off
      while ( (static_cast<unsigned char>(*s) != '\0') &&
            isspace(static_cast<unsigned char>(*s)) ) s++;
    } while (static_cast<unsigned char>(*s) != '\0');

    }

  /* if verbose output is requested:
  std::cerr << "using labels: ";
  for ( unsigned int i = 0; i <= maxLabel; i++ )
    {
    if ( labels[i] ) std::cerr << i << " ";
    }
  std::cerr << std::endl;
  */


  ///
  // load each image, generate stats, and quit
  ///
  const std::vector<std::string>& inputs = optFilesArg.getValue();

  for ( unsigned int subj  = 0; subj < inputs.size(); subj++ )
    {
    // std::cerr << "reading " << inputs[subj] << std::endl;

    // read image
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName( inputs[subj] );
    try 
      {
      reader->Update();
      }
    catch (itk::ExceptionObject& e )
      {
      std::cerr << "error reading input " << inputs[subj] << ": " << e << std::endl;
      }

    ImageType::Pointer image = reader->GetOutput();

    // compute number of voxels in image 
    const ImageType::SizeType& size = image->GetBufferedRegion().GetSize();
    unsigned long totalVoxels = 1;
    for ( unsigned int i = 0; i < size.GetSizeDimension(); i++ )
      {
      totalVoxels *= size[i];
      }

    // computer voxel volume and total image volume
    const ImageType::SpacingType& spacing = image->GetSpacing();
    double voxelVolume = 1.0;
    for ( unsigned int i = 0; i < spacing.Dimension; i++ )
      {
      voxelVolume *= spacing[i];
      }

    double totalVolume = voxelVolume * totalVoxels;

    std::map<unsigned int, unsigned long> count;        // missing keys start out as zero
    ImageConstIterator iter( image, image->GetBufferedRegion() );

    // generate tissue statistics -- we count all labels <= maxLabel, even if
    // they're excluded by labels[].  We exclude printing them.
    for ( iter.GoToBegin(); !iter.IsAtEnd(); ++iter ) 
      {
      unsigned int v = static_cast<unsigned int>( iter.Value() );
      if ( v <= maxLabel )
        {
        count[v]++;
        }
      else
        {
        std::cerr << "warning: ignoring value: " << v << " in image" << std::endl;
        }
      }

    const char sep = ',';               // the separator character

  // print title first time around
    if ( subj == 0 )
      {
      std::cout << "filename" << sep << "0" << sep << "x!=0";
      for ( unsigned int i = 1; i <= maxLabel; i++ )
        if ( labels[i] )
          {
          std::cout << sep << i;
          }
      std::cout << std::endl;
      }

    // output statistics
    std::cout.precision(2);
    std::cout << inputs[subj] << sep;
    std::cout << std::fixed << count[0]*voxelVolume << sep << totalVolume - ( count[0]*voxelVolume );
    for ( unsigned int i = 1; i <= maxLabel; i++ )
      {
      if ( labels[i] )
        {
        std::cout << sep << count[i]*voxelVolume;
        }
      }
    std::cout << std::endl << std::flush;


    }	// for subject...

  return 0;
}

