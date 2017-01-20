/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

//  Software Guide : BeginCommandLineArgs
//    INPUTS:  {FivePointsDilated.png}
//    OUTPUTS: {DanielssonDistanceMapImageFilterOutput1.png}
//    OUTPUTS: {DanielssonDistanceMapImageFilterOutput2.png}
//    ARGUMENTS: {DanielssonDistanceMapImageFilterOutput3.mhd}
//  Software Guide : EndCommandLineArgs

// Software Guide : BeginLatex
//
// This example illustrates the use of the
// \doxygen{DanielssonDistanceMapImageFilter}.  This filter generates a
// distance map from the input image using the algorithm developed by
// Danielsson \cite{Danielsson1980}. As secondary outputs, a Voronoi
// partition of the input elements is produced, as well as a vector image
// with the components of the distance vector to the closest point. The input
// to the map is assumed to be a set of points on the input image. The label
// of each group of pixels is assigned by the
// \doxygen{ConnectedComponentImageFilter}.
//
// \index{itk::Danielsson\-Distance\-Map\-Image\-Filter!Instantiation}
// \index{itk::Danielsson\-Distance\-Map\-Image\-Filter!Header}
//
// The first step required to use this filter is to include its header file.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkBresenhamLine.h"
// Software Guide : EndCodeSnippet

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputVectorMapImageFile ";
    std::cerr << " inputProfileImage";
    std::cerr << " outputLineIntegralFile ";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  //  Software Guide : BeginLatex
  //
  //  Then we must decide what pixel types to use for the input and output
  //  images. Since the output will contain distances measured in pixels, the
  //  pixel type should be able to represent at least the width of the image,
  //  or said in $N$-dimensional terms, the maximum extension along all the dimensions.
  //  The input, output (distance map), and voronoi partition image types are
  //  now defined using their respective pixel type and dimension.
  //
  //  Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  typedef  unsigned char                    InputPixelType;
  typedef           float                   FloatPixelType;
  typedef  unsigned char                    VoronoiPixelType;
  typedef itk::Image< InputPixelType,  3 >  InputImageType;
  typedef itk::Image< FloatPixelType,  3 >  FloatImageType;
  typedef itk::Image< VoronoiPixelType, 3 > VoronoiImageType;
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  The filter type can be instantiated using the input and output image
  //  types defined above. A filter object is created with the \code{New()}
  //  method.
  //
  //  \index{itk::Danielsson\-Distance\-Map\-Image\-Filter!instantiation}
  //  \index{itk::Danielsson\-Distance\-Map\-Image\-Filter!New()}
  //  \index{itk::Danielsson\-Distance\-Map\-Image\-Filter!Pointer}
  //
  //  Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  typedef itk::DanielssonDistanceMapImageFilter<
               InputImageType, FloatImageType, VoronoiImageType >  FilterType;
  // Software Guide : EndCodeSnippet

  //
  // Reader and Writer types are instantiated.
  //
  typedef itk::ImageFileReader< InputImageType  >  ReaderType;

  // Software Guide : BeginCodeSnippet
  typedef FilterType::VectorImageType   OffsetImageType;
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  We can use this type for instantiating an \doxygen{ImageFileWriter} type
  //  and creating an object of this class in the following lines.
  //
  //  Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  typedef itk::ImageFileReader< OffsetImageType >  ReaderOffsetType;
  ReaderOffsetType::Pointer offsetReader = ReaderOffsetType::New();
  // Software Guide : EndCodeSnippet

  offsetReader->SetFileName( argv[1]  );


  //  Software Guide : BeginLatex
  //
  //  Execution of the writer is triggered by the invocation of the
  //  \code{Update()} method. Since this method can potentially throw
  //  exceptions it must be placed in a \code{try/catch} block.
  //
  //  Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  try
    {
    offsetReader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr <<     exp    << std::endl;
    }
  // Software Guide : EndCodeSnippet
  //
  
  //  Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  OffsetImageType::Pointer distanceComponents  =  offsetReader->GetOutput();


  // Read Profile image
  typedef itk::ImageFileReader< FloatImageType >  ProfileReaderType;
  ProfileReaderType::Pointer profilereader = ProfileReaderType::New();
  profilereader->SetFileName( argv[2] );
  profilereader->Update();
  FloatImageType::Pointer     ProfileImage = profilereader->GetOutput();

  // Duplicate Profile image with zeros
  typedef itk::ImageDuplicator< FloatImageType > DuplicatorType;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(ProfileImage );
  duplicator->Update();
  FloatImageType::Pointer LineProfileImage = duplicator->GetOutput();

  //  Software Guide : BeginLatex
  //
  //  Note that only the \doxygen{MetaImageIO} class supports reading and
  //  writing images of pixel type \doxygen{Offset}.
  //
  //  Software Guide : EndLatex

  std::cout << "ComputeLineProfile Start" << std::endl;

  typename FloatImageType::RegionType region  = LineProfileImage->GetRequestedRegion();

  itk::ImageRegionIteratorWithIndex< FloatImageType   >  ot(LineProfileImage,    region);
  itk::ImageRegionIteratorWithIndex< OffsetImageType  >  ct(distanceComponents,  region);

  int count = 0;
  ot.GoToBegin();
  ct.GoToBegin();
  while ( !ot.IsAtEnd() )
    {
    itk::BresenhamLine<3> bhline;
    FilterType::IndexType pixel0 = ct.GetIndex();
    FilterType::IndexType pixel1 = ct.GetIndex() + ct.Get();
 
    std::vector< FilterType::IndexType  > pixels = bhline.BuildLine(pixel0, pixel1);
 
    std::cout  << count << ": " ;
    float lineintegral = 0.0;
    for(unsigned int i = 0; i < pixels.size(); i++)
      {
      lineintegral = lineintegral + ProfileImage->GetPixel(pixels[i]);
      }
    lineintegral = lineintegral /pixels.size();

    if ( region.IsInside(pixel0 ) )
      {
      ot.Set( lineintegral  );
      }

 

    // OffsetType distanceVector = ct.Get();
    // double     distance = 0.0;
    // if ( m_UseImageSpacing )
    //   {
    //   for ( unsigned int i = 0; i < InputImageDimension; i++ )
    //     {
    //     double component = distanceVector[i] * static_cast< double >( m_InputSpacingCache[i] );
    //     distance += component * component;
    //     }
    //   }
    // else
    //   {
    //   for ( unsigned int i = 0; i < InputImageDimension; i++ )
    //     {
    //     distance += distanceVector[i] * distanceVector[i];
    //     }
    //   }

    ++ot;
    ++ct; ++count;
    }
    
  std::cout<< std::endl  << "ComputeLineProfile End" << std::endl;


  // save line profile image
  typedef itk::ImageFileWriter< FloatImageType >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( LineProfileImage  );
  writer->Update();


  return EXIT_SUCCESS;
}
