#include <iostream>
#include <fstream>
#include <sstream>

#include <itkSignedDanielssonDistanceMapImageFilter.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#ifdef WIN32
#include <cstdlib>
#include <ctime>

#define random rand
#endif

int main(int argc, const char *argv[])
{

  static int const ImageDimension = 3;
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ChannelImageType;
  typedef itk::VectorImage<PixelType, ImageDimension> AtlasImageType;
  typedef itk::ImageFileReader<ChannelImageType> ChannelReaderType;
  typedef itk::ImageFileReader<AtlasImageType> AtlasReaderType;

  if ( (argc < 7 ) )
    {
      std::cerr <<
	"Usage: crlatlassample vatlas.nrrd distance[+limit] output-mapfile target_nvoxels labels"
		<< std::endl << " channel1.nrrd [channel2.nrrd...]" << std::endl <<
	"Takes a multi-component atlas (probability) volume and writes out a map file"
		<< std::endl <<
	"with random samples generated from the channel volumes for each"
		<< std::endl <<
	"label in the atlas."
		<< std::endl << std::endl <<
	"target_nvoxels is a desired number of samples from each class"
		<< std::endl <<
	"labels is either the word 'all' (no quotes) or a quoted list of labels"
		<< std::endl <<
	"distance is the distance from class boundaries (ML) in the atlas to stay away from and, optionally, " 
		<< std::endl << 
	" limit states how far beyond that distance to go FOR CLASS 0 only!"
		<< std::endl
		<< "HEY!!!!!!!!!!!!!!! READ THIS---> only works for two class 0 and non-zero right now"
		<< std::endl
		<< "also, it's possible for this to hang if there aren't enough of a given class to sample"
		<< std::endl;

      return 1;
    }

  const int numchannels = argc - 6;
  const char* const atlasfile = argv[1];

  // parsing the distance in argv2[]
  std::istringstream temp(argv[2]);
  float distance = 0.;
  float limit = 0.; 		// flag value: 0.0 is "not set"
  temp >> distance;
  if ( !temp.eof() && temp.peek() == '+' )
    { 	// optional limit distance specified
    temp.ignore();
    temp >> limit;
    }
  

  const char* const mapfile = argv[3];
  const int ndesired = atoi( argv[4] );
  const char* const labelstr = argv[5];
  const char** const chfiles = argv + 6;

  std::cout << "output will be written to " << mapfile << std::endl;
  std::cout << "desired number of samples per class is " << ndesired << std::endl;
  std::cout << "distance from maximum-likelihood boundaries: " << distance << std::endl;
  if ( limit != 0.0 )
    std::cout << "... limit distance from ML boundaries for class 0 is " << limit << std::endl;


  // load all of the relevant volumes and check that their sizes match
  AtlasReaderType::Pointer atlasReader = AtlasReaderType::New();
  AtlasImageType::Pointer atlas;

  std::vector<ChannelReaderType::Pointer> channelReaders(numchannels);
  for (int n = 0; n < numchannels; n++)
    channelReaders[n] = ChannelReaderType::New();

  std::vector<ChannelImageType::Pointer> channel(numchannels);


  /// load the atlas

  std::cout << "reading atlas volume " << atlasfile << "...";

  atlasReader->SetFileName(atlasfile);
  try
    {
      atlasReader->Update();
    }
  catch (...)
    {
      std::cerr << "failure reading atlas volume: " << atlasfile << std::endl;
      return 1;
    }

  std::cout << std::endl;
  atlas = atlasReader->GetOutput();

  /// load the channels

  for ( int ch = 0; ch < numchannels; ch++ )
    {
      std::cout << "reading channel " << ch+1 << ": " << chfiles[ch] << "...";
      channelReaders[ch]->SetFileName(chfiles[ch]);
      try
	{
	  channelReaders[ch]->Update();
	}
      catch (...)
	{
	  std::cerr << "failure reading input channel: " << chfiles[ch] <<
	    std::endl;
	  return 1;
	}
      channel[ch] = channelReaders[ch]->GetOutput();
      std::cout << std::endl;
    }


  /// check the sizes and that each volume matches

  unsigned int size[3];
  for (int i = 0; i < 3; i++) size[i] = atlas->GetLargestPossibleRegion().GetSize()[i];

  std::cout << "atlas shape is " << size[0] << " x " << size[1] << " x " << size[2] << std::endl;

  for ( int ch = 0; ch < numchannels; ch++ )
    {
      if ( size[0] != channel[ch]->GetLargestPossibleRegion().GetSize()[0] ||
	   size[1] != channel[ch]->GetLargestPossibleRegion().GetSize()[1] ||
	   size[2] != channel[ch]->GetLargestPossibleRegion().GetSize()[2] )
	{
	  std::cerr << "channel " << ch+1 << ", file " << chfiles[ch] <<
	    " is not the same shape as the atlas.  Exiting." << std::endl;
	  return 1;
	}

      if ( channel[ch]->GetNumberOfComponentsPerPixel() != 1 )
	{
	  std::cerr << "WARNING: channel " << ch+1 << " has more than one component." <<
	    std::endl << "We're ignoring past the first component of each." << std::endl;
	}
    }

  std::cout << "each of the channels matches" << std::endl;

  /// generate a maximum-likelihood segmentation of the atlas
  ChannelImageType::Pointer mlseg = ChannelImageType::New();
  mlseg->CopyInformation(channel[0]);
  mlseg->SetRegions(channel[0]->GetLargestPossibleRegion());
  mlseg->SetNumberOfComponentsPerPixel(1);
  mlseg->Allocate();

  for ( int z = 0; z < size[2]; z++ )
    for ( int y = 0; y < size[1]; y++ )
      for ( int x = 0; x < size[0]; x++ )
        {
	  AtlasImageType::IndexType index;
	  index[0] = x; index[1] = y; index[2] = z;
	  itk::VariableLengthVector<float> v = atlas->GetPixel(index);
	  float max = v[0];
	  int imax = 0;
	  for (int comp = 1; comp < atlas->GetNumberOfComponentsPerPixel(); comp++)
	    {
	      if ( v[comp] > max )
		{
		  max = v[comp];
		  imax = comp;
		}
	    }
	  mlseg->SetPixel(index,imax);       // PLUS ONE TO MAKE THE DISTANCE MAP HAPPY
        }

  /// generate a distance transform of the ML segmentation
  typedef itk::SignedDanielssonDistanceMapImageFilter<ChannelImageType,ChannelImageType> FilterType;

  FilterType::Pointer filter = FilterType::New();
  std::cout << "(using image spacing for distance map)" << std::endl;
  filter->UseImageSpacingOn();
  filter->SetInput( mlseg );
  filter->Update();

  ChannelImageType::Pointer dmap = filter->GetDistanceMap();

#if 0
  /// write stuff for debugging
  typedef itk::ImageFileWriter< ChannelImageType > ChannelWriterType;
  ChannelWriterType::Pointer writer = ChannelWriterType::New();

  try
    {
      writer->SetInput(mlseg);
      writer->SetFileName("mlseg.nrrd");
      writer->Update();
      writer->SetInput(dmap);
      writer->SetFileName("dmap.nrrd");
      writer->Update();
    }
  catch( itk::ExceptionObject exp)
    {
      std::cerr << "itk exception caught! " << exp << std::endl;
      exit(1);
    }

#endif

  /// sample the sample and write the map file


  //// open the output file
  //// yeah, this is so much better than C (not)
  std::ofstream map;
  map.exceptions( std::ofstream::badbit | std::ofstream::failbit
		  | std::ofstream::eofbit );
  try { map.open( mapfile ); }
  catch ( std::ios_base::failure& e )
    {
      std::cerr << "failure opening map file " << mapfile << " for output " << std::endl;
      return 1;
    }

  std::vector<bool> luse(atlas->GetNumberOfComponentsPerPixel(), false);

  std::cout << "using labels: ";
  if ( strcmp( labelstr, "all" ) != 0 )         // specified label restrictions
    {
      const char* head;
      const char* tail;
      head = tail = labelstr;
      while ( *head )
	{
	  while ( *tail && *tail++ != ' ' )         // find space or end of string
	    ;
	  if ( head != tail )
	    {
	      int label = atoi( head );
	      luse[label] = true;
	      std::cout << label << " ";
	    }
	  head = tail;
	}
    }
  else
    {
      for (int i = 0; i < atlas->GetNumberOfComponentsPerPixel(); i++)
	{
	  luse[i] = true;
	  std::cout << i << " ";
	}
    }
  std::cout << std::endl;


  ChannelImageType::Pointer samples = ChannelImageType::New();
  samples->CopyInformation(channel[0]);
  samples->SetRegions(channel[0]->GetLargestPossibleRegion());
  samples->SetNumberOfComponentsPerPixel(1);
  samples->Allocate();

#ifdef WIN32
	time_t t;
	time(&t);
	srand(t);
#else
	srandom( getpid() );
#endif

  std::vector<int> lsampled(atlas->GetNumberOfComponentsPerPixel(), 0);

  bool done = false;

  while ( ! done )
    {
      // generate random coordinate
      long int x = static_cast<long int>(size[0]-1) * random() / RAND_MAX;
      long int y = static_cast<long int>(size[1]-1) * random() / RAND_MAX;
      long int z = static_cast<long int>(size[2]-1) * random() / RAND_MAX;

      // we subtract one here because we add one to the label when
      // creating the mlseg since the distance map ignores zero
      // (background), which we're possibly interested in
      ChannelImageType::IndexType pixelIndex;
      pixelIndex[0] = x; pixelIndex[1] = y; pixelIndex[2] = z;
      int label = static_cast<int>(mlseg->GetPixel(pixelIndex));

      // is that label being scanned
      if ( luse[label] && lsampled[label] < ndesired )
        {
	  // are we above distance threshold
	  float dist = fabs(dmap->GetPixel( pixelIndex ));
	  /*std::cerr << "label " << label << ": " << x << "," << y << ","
		    << z << ": ";
		    std::cerr << " distance " << dist;*/
	  // if label is zero and limit is not precisely 0.0 (therefore
	  // set on the cmdline), then we check it.
	  if ( dist > distance && ( label > 0 || limit == 0. || dist <= limit+distance ) )
	    {
	      samples->SetPixel(pixelIndex,label+1);
	      // okay, print a sample
	      map << label << std::endl;
	      for ( int ch = 0; ch < numchannels; ch++ )
		{
		  map << static_cast<int>( channel[ch]->GetPixel(pixelIndex)) << " ";
		}
	      map << ":row-" << y << "-col-" << x << "-sliceindex-" << z << std::endl;
	      // increse count for this label
	      lsampled[label]++;
	    }
	  //std::cerr << std::endl;

        }

      // check to see if we're done
      done = true;
      for (int i = 0; i < atlas->GetNumberOfComponentsPerPixel(); i++)
        {
	  if ( luse[i] && lsampled[i] < ndesired )
	    {
	      done = false;
	      break;
	    }
        }
    }

  //samples.writeVolume("samplemap.crl");

  map.close();

  return 0;
}
