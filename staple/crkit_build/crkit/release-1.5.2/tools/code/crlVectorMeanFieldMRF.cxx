
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkVariableLengthVector.h"
#include "itkVectorImage.h"
#include "itkConstantBoundaryCondition.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkMultiThreader.h"
#include "itkFastMutexLock.h"
#include "itkBarrier.h"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>

double* initializeInteractionMatrix( const char* argv2, const unsigned int vectorSize );
ITK_THREAD_RETURN_TYPE singleMethod(void* data);

typedef double ProbabilityType;
typedef itk::VariableLengthVector<ProbabilityType > VectorType;

static int const ImageDimension = 3; 

typedef itk::VectorImage< ProbabilityType, ImageDimension > ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
typedef itk::ImageFileWriter< ImageType > WriterType;

//typedef itk::FastMutexLock MutexType;
//typedef itk::Barrier BarrierType;


struct ThreadUserData
{
  ImageType::Pointer curnEstimate;
  ImageType::Pointer prevEstimate;
  ImageType::SizeType inSize;
  double* interactionmatrix;
  double mrfboost;
  unsigned int vectorSize;
};

int main( int argc, char ** argv ) 
{
  bool optThreading = true;
  
  if ( argc != 7 && argc != 8 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " VectorPDFInputImage interactionExampleSegmentation|'auto' "
      << " convergenceTolerance mrfStrength maxIterationCount outputImageFile [numthreads]"
      << std::endl;
    std::cerr << "This is a work in progress." << std::endl;
      return -1;
    }

  int numThreads = 1;	
  if ( argc == 8 )
    {
    std::istringstream s(argv[7]);

    if ( ! ( s >> numThreads ) )
      {
      std::cerr << "can't parse thread argument" << std::endl;
      return 1;
      }

    }
  std::cout << "requesting " << numThreads << " threads." << std::endl;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  double mrfboost = std::atof( argv[4] );
  signed int maxIterationCount = std::atoi( argv[5] );
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[6] );

  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }


  ImageType::Pointer prevEstimate;
  ImageType::Pointer curnEstimate;
  ImageType::Pointer tmpImagePointer;

  prevEstimate = reader->GetOutput();
  prevEstimate->DisconnectPipeline();
  
  const unsigned int vectorSize = prevEstimate->GetVectorLength();
  std::cerr << "probability vector size is " << vectorSize << std::endl;
  
  curnEstimate = ImageType::New();
  curnEstimate->SetRegions( prevEstimate->GetLargestPossibleRegion() );
  curnEstimate->CopyInformation( prevEstimate );
  curnEstimate->SetVectorLength( vectorSize );
  curnEstimate->Allocate();
  


  ImageType::SizeType inSize;
  inSize = prevEstimate->GetLargestPossibleRegion().GetSize();
  
  unsigned long nvoxels = 1;
  for ( signed int d = 0; d < ImageDimension; d++ )
    {
    nvoxels *= inSize[d];
    }
#ifndef WIN32
	unsigned long nvoxelsthread = static_cast<unsigned long>(vcl_ceil(nvoxels / numThreads ));
#else
	unsigned long nvoxelsthread = static_cast<unsigned long>(vcl_ceil((double)nvoxels / (double)numThreads ));
#endif
  std::cerr << "numvoxels=" << nvoxels << std::endl;
  std::cerr << "numvoxelsthread=" << nvoxelsthread << std::endl;
  
  double* interactionmatrix = initializeInteractionMatrix( argv[2], vectorSize );
  
  if (!interactionmatrix)
    {
    std::cerr << "failed initializing the interaction matrix.  Exiting" << std::endl;
    return EXIT_FAILURE;
    }
  

  // Repeatedly iterate over the calculation until we run out of iterations.
  signed int thisIter = 1;

  
  do 
    {
    std::cerr << "Iteration #" << thisIter << std::endl;
    
    ThreadUserData tud;
    tud.mrfboost = mrfboost;
    tud.interactionmatrix = interactionmatrix;
    tud.curnEstimate = curnEstimate;
    tud.prevEstimate = prevEstimate;
    tud.vectorSize = vectorSize;
    tud.inSize = inSize;
    
    if ( optThreading )
      {
      typedef itk::MultiThreader MultiThreadType;
      MultiThreadType::Pointer threader = MultiThreadType::New();
      if (numThreads) threader->SetNumberOfThreads(numThreads);
      threader->SetSingleMethod(singleMethod,&tud);
      threader->SingleMethodExecute();     
      }
    else
      {
      throw "not implemented";
      }

    tmpImagePointer = curnEstimate;
    curnEstimate = prevEstimate;
    prevEstimate = tmpImagePointer;

  } while (thisIter++ < maxIterationCount);
  
  // Now prevEstimate holds the most up to date estimate of the mean field (because we swapped)

  writer->SetInput( prevEstimate );
  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }

  return EXIT_SUCCESS;

}


ITK_THREAD_RETURN_TYPE singleMethod(void* data)
{
  typedef itk::MultiThreader::ThreadInfoStruct ThreadInfoType;
  
  ThreadInfoType* tis = static_cast<ThreadInfoType*>(data);
  
  std::cerr << "Thread info: thread id: " << tis->ThreadID << std::endl;
//  std::cerr << "number of threads: " << tis->NumberOfThreads << std::endl;
  
  ThreadUserData *tudp = static_cast<ThreadUserData*>(tis->UserData);
  
//  std::cerr << "mrfboost is " << tudp->mrfboost << std::endl;
  
  
  //unsigned long voxelcount = 0;
  IteratorType inIt( tudp->prevEstimate, tudp->prevEstimate->GetLargestPossibleRegion());
  
  inIt.GoToBegin();
  for ( unsigned int tmp = tis->ThreadID; !inIt.IsAtEnd() && tmp; --tmp )
    ++inIt;
  
  for ( ; !inIt.IsAtEnd();  )
    {
    ImageType::IndexType inIdx = inIt.GetIndex();
    ImageType::IndexType neighborIdx = inIdx;
    VectorType inVoxel = inIt.Get();
    VectorType outVoxel = tudp->curnEstimate->GetPixel(inIdx);
    for (unsigned int i = 0; i < outVoxel.Size(); i++)
      outVoxel[i] = inVoxel[i];
    
    // Hard coded for a particular neighborhood size, and for 3D.
    for (int rk = -1; rk <= 1; rk++) 
	{
	for (int rj = -1; rj <= 1; rj++) 
	  {
	  for (int ri = -1; ri <= 1; ri++) 
	    {
	    if ( (inIdx[2] == 0) && (rk == -1) ) continue;
	    if ( (inIdx[2] == (signed)(tudp->inSize[2]-1)) && (rk == +1) ) continue;
	    if ( (inIdx[1] == 0) && (rj == -1) ) continue;
	    if ( (inIdx[1] == (signed)(tudp->inSize[1]-1)) && (rj == +1) ) continue;
	    if ( (inIdx[0] == 0) && (ri == -1) ) continue;
	    if ( (inIdx[0] == (signed)(tudp->inSize[0]-1)) && (ri == +1) ) continue;
	    neighborIdx[0] = inIdx[0] + ri;
	    neighborIdx[1] = inIdx[1] + rj;
	    neighborIdx[2] = inIdx[2] + rk;
	    // mrfboost plays the role of inverse temperature, and modifies
	    // the strength of the spatial prior.
	    VectorType neighbor = tudp->prevEstimate->GetPixel(neighborIdx);
	    for ( unsigned int s = 0; s < tudp->vectorSize; s++ )
	      for ( unsigned int n = 0; n < tudp->vectorSize; n++ )
      		{
      		outVoxel[s] *= exp( log( tudp->interactionmatrix[s*tudp->vectorSize + n]) * tudp->mrfboost * neighbor[n ]);
      		}
//           tmp *= exp( log(interactionmatrix[s*numcomp + n])*mrfboost * 
//                  prevEstimate->GetPixel(neighborIdx) );
	    }
	  }
	}
    
    double sum = 0.0;
    for ( unsigned int s = 0; s < tudp->vectorSize; s++ )
	{
	sum += outVoxel[s];
	}
    if ( sum <= 0.0 )
	{
	sum = 1.0;
	std::cerr << "Found a sum of " << sum << " at " 
	                << inIt.GetIndex() << std::endl;
	}
    outVoxel /= sum;
//    tudp->curnEstimate->SetPixel(inIdx,tmp);
    
    for ( unsigned int tmp = tis->NumberOfThreads; !inIt.IsAtEnd() && tmp; --tmp )
      ++inIt;
    
    }  // for each voxel
  
  
  return ITK_THREAD_RETURN_VALUE;
}

// http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
void tokenizeString( const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " " )
{
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (std::string::npos != pos || std::string::npos != lastPos)
    {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
    }
}


double* readInteractionMatrixFile( const char* argv2, unsigned int& matsize )
{
  std::ifstream infile( argv2 );
  std::string line;

  if ( !infile.good() )
    {
    std::cerr << "error opening input file " << argv2 << std::endl;
    return 0;
    }

  typedef std::vector< std::vector<double> > VVListType;
  VVListType mmatrix;
  matsize = 0;
  
  while (!std::getline(infile, line).eof())
    {
    std::vector<std::string > tokens;
    tokenizeString( line, tokens, " \t,");

    if ( tokens.size() == 0)
      continue;
    
    std::vector<double > ptokens;
    for ( std::vector<std::string >::iterator i = tokens.begin(); i != tokens.end(); ++i )
      {
      double temp;
      if ( std::istringstream( *i ) >> temp )
	ptokens.push_back( temp );
      }
    
    if ( matsize == 0 ) matsize = ptokens.size();
    
    if ( ptokens.size() > 0 )		// if line wasn't blank
      {
      if ( ptokens.size() == matsize )	// we're trying to make sure that each row is same length
	mmatrix.push_back( ptokens );
      else
	{
	std::cerr << "last row was size: " << matsize << " and this row is " << ptokens.size() << std::endl;
	std::cerr << "this is bad.  we're quitting now." << std::endl;
	std::cerr << "the offending row was: " << std::endl;
	std::cerr << line << std::endl;
	return 0;
	}
      }

    }
  
  if ( mmatrix.size() != matsize )
    {
    std::cerr << "read matrix was NOT SQUARE.  We read " << mmatrix.size() << " rows and " << matsize << " cols.  We're quitting." << std::endl;
    return 0;
    }

  std::cerr << "interaction matrix file " << argv2 << " parsed as: " << std::endl;

  double* interactionMatrix = new double[matsize*matsize];
  
  unsigned int row = 0;
  for (  VVListType::iterator i = mmatrix.begin() ; i != mmatrix.end(); ++i )
    {
    std::copy(i->begin(), i->end(), std::ostream_iterator<double>(std::cerr, ", ")); 
    std::cerr << std::endl;
    for ( unsigned int col = 0; col < matsize; col++ )
      interactionMatrix[row*matsize+col] = (*i)[col];
    ++row;
    }

  infile.close();
  return interactionMatrix;
}

double* initializeInteractionMatrix( const char* argv2, const unsigned int vectorSize )
{
  double* interactionmatrix = new double[vectorSize*vectorSize];
  for ( unsigned int i = 0; i < vectorSize*vectorSize; i++ )
    interactionmatrix[i] = 0.0;
  
  if (!strcmp(argv2,"auto"))
    {
    std::cout << "Automatic initialiation of the interaction matrix " <<
       " with size " << 
                (vectorSize * 
                 vectorSize) << std::endl;

    for (unsigned int i = 0; i < vectorSize; i++) 
      	{
    	for (unsigned int j = 0; j < vectorSize; j++) 
    	  {
    	  if (i == j) 
    	    {
    	    interactionmatrix[i*vectorSize + j] = 0.95;
    	    } 
    	  else 
    	    {
    	    interactionmatrix[i*vectorSize + j] = 
                                    0.05/(vectorSize-1);
    	    }
    	  std::cout << "intmat["<<i<<","<<j<<"] " << 
    	  interactionmatrix[i*vectorSize + j] << std::endl;
    	  }
      	}
    }
  else if (!strcmp(argv2+strlen(argv2)-4, ".txt"))
    {
    unsigned int matsize;
    interactionmatrix = readInteractionMatrixFile(argv2, matsize);
    if (!interactionmatrix) return 0;
    if ( matsize != vectorSize )
      {
      std::cerr << "matrix size read from file was " << matsize << "x" << matsize << " but vector size from PDF is " << vectorSize << std::endl;
      return 0;
      }
    }
  else
    {
    std::cerr << "Attempting initialization of interaction matrix from segmentation" << std::endl;
    // okay, so we don't recognize any special words in argv2, so we assume that it's a filename of an example segmentation to parse.
    typedef unsigned char SegPixelType;
    typedef itk::Image<SegPixelType, 3> SegImageType;
    typedef itk::ImageFileReader<SegImageType > SegImageReaderType;
    
    SegImageReaderType::Pointer reader = SegImageReaderType::New();
    reader->SetFileName(argv2);
    try
      {
      reader->Update();
      }
    catch (itk::ExceptionObject& err )
      {
      std::cerr << "failed reading example segmentation file: " << argv2 << ", maybe you want to specify 'auto'" << std::endl;
      delete [] interactionmatrix;
      return 0;
      }
    
    // check image here for bad labels first...
    typedef itk::ImageRegionIteratorWithIndex< SegImageType > TempIteratorType;
    
    TempIteratorType it1( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
    
    for ( it1.GoToBegin(); !it1.IsAtEnd(); ++it1 )
      {
      if ( it1.Value() >= vectorSize )
	{
	std::cerr << "bad label " << it1.Value() << " in example segmentation at index " << std::endl << it1.GetIndex() << std::endl;
	std::cerr << "changed to zero." << std::endl;
	it1.Value() = 0;
	}
      }
    
    
    // ...then walk it with confidence.
    typedef itk::ConstantBoundaryCondition< SegImageType > BoundaryType;
    typedef itk::ConstNeighborhoodIterator< SegImageType, BoundaryType > SegIteratorType;
    SegIteratorType::RadiusType radius;
    radius.Fill(1U);	// one voxel out in all directions
    SegIteratorType it( radius, reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
    BoundaryType bounds;
    bounds.SetConstant(0U);
    it.OverrideBoundaryCondition(&bounds);
    
    SegPixelType center;
    SegPixelType neighbor;
    
    // calculate the size of the neighborhood.
    SegIteratorType::SizeType faces = it.GetSize();
    SegIteratorType::SizeValueType nfaces = 1;
    double count = 0.0;
    for ( unsigned int i = 0; i < faces.GetSizeDimension(); ++i )
      nfaces *= faces[i];
    std::cerr << "number of faces in the neighborhood is " << nfaces << std::endl;
    for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      center = it.GetCenterPixel();
      for ( unsigned int i = 0; i < nfaces; ++i )
	{
	neighbor = it.GetPixel(i);
	interactionmatrix[center*vectorSize+neighbor] += 1.0;
	interactionmatrix[neighbor*vectorSize+center] += 1.0;
	count += 2.0;
	}
      interactionmatrix[center*vectorSize+center] -= 2.0;		// we'll have visited the center once as a neighbor.
      count-=2.0;
      }
    
    for ( unsigned int i = 0; i < vectorSize; i++ )
      {
      for ( unsigned int j = 0; j < vectorSize; j++ )
        {
        if ( interactionmatrix[i*vectorSize+j] == 0.0 )
          {
          interactionmatrix[i*vectorSize+j] = 1.0;
          count += 1.0;
          }
        }
      }
    
    for ( unsigned int i = 0; i < vectorSize; i++ )
      {
      for ( unsigned int j = 0; j < vectorSize; j++ )
        {
        interactionmatrix[i*vectorSize+j] /= count;
        }
      }
    }
  

  
  
#if 0
  // renormalize matrix
  for ( unsigned int i = 0; i < vectorSize; i++)
    {
    double sum = 0.0;			// we want to normalize each row to sum to one.
    for ( unsigned int j = 0; j < vectorSize; j++ )
	{
	sum += interactionmatrix[i*vectorSize+j];
	}
    if ( sum == 0.0 )			// but if the row sums to zero, we initialize it sensibly
	{
	for ( unsigned int j = 0; j < vectorSize; j++ )
	  {
	  if ( j == i )
	    interactionmatrix[i*vectorSize+j] = 0.9999;
	  else
	    interactionmatrix[i*vectorSize+j] = (1.0 - 0.9999) / static_cast<double>(vectorSize);
	  }
	}
    else
	{
	for ( unsigned int j = 0; j < vectorSize; j++ )		// so the row didn't sum to zero, but we check for zeros as elements
	  {
	  if (interactionmatrix[i*vectorSize+j] == 0.0)		// we found an isolated (not the entire row) zero element
	    {
	    interactionmatrix[i*vectorSize+j] = 0.000001;	// so we make it a small non-zero number.
	    sum += 0.000001;					// and increase the sum accordingly
	    }
	  }
	for ( unsigned int j = 0; j < vectorSize; j++ )
	  {
	  interactionmatrix[i*vectorSize+j] /= sum;		// finally, we can normalize by the (possibly) adjusted sum
	  }
	}
    }
#endif
  
  std::cerr << "interaction matrix computed from: " << argv2 << ": " << std::endl;
  for ( unsigned int i = 0; i < vectorSize; i++ )
    {
    for ( unsigned int j = 0; j < vectorSize; j++ )
      {
      std::cerr << interactionmatrix[i*vectorSize+j] << "  ";
      }
    std::cerr << std::endl;
    }
  
  return interactionmatrix;
}

