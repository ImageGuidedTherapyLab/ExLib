#ifndef CRLRELAXATIONSCHEME_H_
#define CRLRELAXATIONSCHEME_H_

#define EXECUTE1 1

#include "itkMacro.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkImageRegionIterator.h"
#include "itkVariableLengthVector.h"
#include "itkMatrix.h"
#include <map>
#include <cassert>

struct classcomp {
  bool operator() (const char& lhs, const char& rhs) const
  {return lhs<rhs;}
};


namespace crl
{

struct Triplet
{
  struct comp
  {
    bool operator() ( const Triplet& a, const Triplet& b ) const
    {
      if ( a.first < b.first )
        return true;
      else if ( a.first > b.first )
        return false;
      else if ( a.second < b.second )
        return true;
      else if ( a.second > b.second )
        return false;
      else if ( a.third <= b.third )
        return true;
      else
        return false;
    };
  };

  void Set(unsigned int _first, unsigned int _second, unsigned int _third )
  {
  unsigned int tmp;
  first = _first;
  second = _second;
  third = _third;
  if ( second < first || third < first )
    {
    if ( second < third )
      {
      tmp = first;
      first = second;
      second = tmp;
      }
    else
      {
      tmp = first;
      first = third;
      third = tmp;
      }
    }
  if ( third < second )
    {
    tmp = second;
    second = third;
    third = tmp;
    }
  }

  unsigned int first;
  unsigned int second;
  unsigned int third;
};

struct Pair
{
  struct comp
  {
    bool operator() ( const Pair& a, const Pair& b ) const
      {
      if ( a.first < b.first )
        return true;
      else if ( a.first > b.first )
        return false;
      else if ( a.second <= b.second )
        return true;
      else
        return false;
      };
  };
  unsigned int first;
  unsigned int second;
};


class RelaxationScheme
{
public:
  // use the following static function (defined at the end of the file) in order to interrogate the input image and allocate the
  // proper templated algorithm based on the image dimension.  You can then use the pointer to the base class generically.
  //
  // So you do something like this:
  // crl::RelaxationScheme* diff = crl::RelaxationScheme::ProbeFile( inputfile );
  // diff->ReadFile(inputfile);
  // diff->Execute();
  // diff->WriteFile(outputfile);
  // delete diff;
  static RelaxationScheme* ProbeFile( const std::string& inputfile );


  virtual void ReadFile(const std::string&) = 0;
  virtual void ReadExampleImage(const std::string&) = 0;
  virtual void Execute() = 0;
  virtual void Execute1() = 0;
  virtual void WriteFile(const std::string&, const std::string& = "") = 0;
  virtual void SetIterations(unsigned int& i)  = 0;
};

template <const unsigned int ImageDimension = 3>
class RelaxationSchemeSpecial : public RelaxationScheme
{
public:
  typedef float ElementType;
  typedef typename itk::VariableLengthVector< ElementType >	VectorType;
  typedef typename itk::VariableLengthVector< double >          AccVectorType;
  typedef typename itk::VectorImage<ElementType, ImageDimension >	ImageType;
  typedef unsigned char	LabelType;
  typedef typename itk::Image< LabelType, ImageDimension > LabelImageType;

  // recommended timestep <= ( 0.5 / 2.0^ImageDimension )
  // for 2D: 0.125
  // for 3D: 0.0625
  RelaxationSchemeSpecial() : m_InputImage(0),m_OutputImage(0),m_DebugImage(0),m_Iterations(5) {};

  void ReadFile(const std::string& inputfile)
    {
    typedef itk::ImageFileReader<ImageType > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(inputfile);
    reader->Update();
    m_InputImage = reader->GetOutput();
    m_InputImage->DisconnectPipeline();
    };

  void ReadExampleImage( const std::string& examplefile )
    {
    typedef itk::ImageFileReader<LabelImageType > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName( examplefile ) ;
    reader->Update();
    m_ExampleImage = reader->GetOutput();
    m_ExampleImage->DisconnectPipeline();
    }

  // Execute1() -- this is a different strategy, where an example image is provided and the update rule is based on statistics
  // for all labels
  void Execute1()
    {
#if EXECUTE1
    if ( !m_ExampleImage )
      {
      throw std::runtime_error("BUG: example image not read but Execute1 called");
      }

#ifndef WIN32
    #warning "this is not a finished product"
#endif
    unsigned int maxlabel = 8;
    unsigned int nlabels = maxlabel+1;
    unsigned int radius = 1;
    double *triplets = new double[ nlabels*nlabels*nlabels ];
    double *pairs = new double [ nlabels*nlabels ];
    for ( unsigned int i = 0; i < nlabels*nlabels*nlabels; i++ )
      triplets[i] = 0.0;
    for ( unsigned int i = 0; i < nlabels*nlabels; i++ )
      pairs[i] = 0.0;

    // generate statistics based on example image

    typedef typename itk::ConstantBoundaryCondition< LabelImageType > BoundaryType;
    typedef typename itk::ConstNeighborhoodIterator< LabelImageType > NeighborhoodIteratorType;
    typename NeighborhoodIteratorType::RadiusType nRadius;
    nRadius.Fill(radius);
    NeighborhoodIteratorType nit( nRadius, m_ExampleImage, m_ExampleImage->GetBufferedRegion() );
    BoundaryType bounds;
    bounds.SetConstant(0U);
    nit.OverrideBoundaryCondition(&bounds);

    unsigned int nlength = 1;
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      nlength *= nit.GetSize(d);
      }
    std::cerr << "nlength = " << nlength << std::endl;

    std::cerr << "Generating neighborhood statistics..." << std::endl;

    unsigned long count = 0;
    for ( nit.GoToBegin(); !nit.IsAtEnd(); ++nit )
      {
      ++count;
      if ( count % 100000 == 0 )
	std::cerr << "count=" << count << std::endl;

      typename LabelImageType::PixelType thisi = nit.GetCenterPixel();
      unsigned int i = nlength/2; 	// purposeful integer division

      // find all possible pairs != center and update a triplet to reflect.
      for ( unsigned int j = 0; j < nlength-1; j++ )
	{
	if ( j == i ) continue;		// skip j==i
	typename LabelImageType::PixelType thisj = nit.GetPixel(j);

	for ( unsigned int k = j+1; k < nlength; k++ )
	  {
	  if ( k == i ) continue; 	// skip k==i
	  typename LabelImageType::PixelType thisk = nit.GetPixel(k);
	  triplets[thisi+nlabels*(thisj+nlabels*thisk)] += 1.0;
	  pairs[thisj+nlabels*thisk] += 1.0;
	  }
	}
      }


      // print the map for debugging purposes and "fold" the diagonal matrices
      // i.e. p(x,y) <- p(x,y)+p(y,x) if x != y
      //      p(y,x) <- 0             if x != y
      std::cerr << "Printing triplets and folding diagonal matrix..." << std::endl;

	for (unsigned int lj = 0; lj < nlabels; lj++ )
	  for (unsigned int lk = lj; lk < nlabels; lk++ )
	    {
	    if ( lj != lk )
	      {
	      pairs[lj+nlabels*lk] += pairs[lk+nlabels*lj];
	      pairs[lk+nlabels*lj] = 0.0;
	      }
	    }
      for (unsigned int lj = 0; lj < nlabels; lj++ )
	for (unsigned int lk = lj; lk < nlabels; lk++ )
	  for (unsigned int li = 0; li < nlabels; li++ )
	    {
	    if ( lj != lk )
	      {
	      triplets[li+nlabels*(lj+nlabels*lk)] +=
		triplets[li+nlabels*(lk+nlabels*lj)];
	      triplets[li+nlabels*(lk+nlabels*lj)] = 0.0;
	      }
	    if ( pairs[lj+nlabels*lk] > 0.0 )
	      triplets[li+nlabels*(lj+nlabels*lk)] /= pairs[lj+nlabels*lk];
	    double val = triplets[li+nlabels*(lj+nlabels*lk)];
	    if ( val != 0.0 )
	      std::cerr << "p( " << li << "|" << lj << "," << lk << ")=" << val << std::endl;
	    }

      std::cerr << "hacking partial volume targets" << std::endl;
      // p(l7 | l4,l5)
      double weight = triplets[7+nlabels*(4+nlabels*5)];
      triplets[7+nlabels*(4+nlabels*5)] /= 10000;
      weight -= triplets[7+nlabels*(4+nlabels*5)];
      for ( unsigned int li = 0; li < nlabels; li++ )
        if ( li != 7 )
          triplets[li+nlabels*(4+nlabels*5)] += weight / (nlabels - 1);
      std::cerr << "p(7|4,5)=" << triplets[7+nlabels*(4+nlabels*5)] << std::endl;

      weight = triplets[7+nlabels*(5+nlabels*8)];
      triplets[7+nlabels*(5+nlabels*8)] /= 10000;
      weight -= triplets[7+nlabels*(5+nlabels*8)];
      for ( unsigned int li = 0; li < nlabels; li++ )
        if ( li != 7 )
          triplets[li+nlabels*(5+nlabels*8)] += weight / (nlabels - 1);
      std::cerr << "p(7|5,8)=" << triplets[7+nlabels*(5+nlabels*8)] << std::endl;

    // allocate output image
    m_OutputImage = ImageType::New();
    m_OutputImage->CopyInformation( m_InputImage );
    m_OutputImage->SetRegions( m_InputImage->GetBufferedRegion() );
    m_OutputImage->SetVectorLength( m_InputImage->GetVectorLength() );
    m_OutputImage->Allocate();

    // allocate mask image
    typename LabelImageType::Pointer maskImage = LabelImageType::New();
    maskImage->CopyInformation( m_InputImage );
    maskImage->SetRegions( m_InputImage->GetBufferedRegion() );
    maskImage->Allocate();
    maskImage->FillBuffer(1U);



    // normalize the input image so each vector sums to one.
    // and populate the mask image
    typedef itk::ImageRegionIterator< ImageType > InputIteratorType;
    InputIteratorType it( m_InputImage, m_InputImage->GetBufferedRegion() );
    typedef itk::ImageRegionIterator< LabelImageType > LabelIteratorType;
    LabelIteratorType mit( maskImage, maskImage->GetBufferedRegion() );

    for ( it.GoToBegin(),mit.GoToBegin(); !it.IsAtEnd(); ++it,++mit )
      {
      VectorType v = it.Get();
      double sum  = 0.0;
      for ( unsigned int i = 0; i < v.GetSize(); i++ )
        sum += v[i];
      v /= sum;
      it.Set(v);                // probably not necessary the way VectorImage works, but whatever.
      if ( v[0] > 0.99999 )
        mit.Set(0U);            // if we're sure that it's background, save us the pain of optimizing it
      }




    while ( m_Iterations-- > 0 )
      {
      std::cerr << "Iterations remaining: " << m_Iterations+1 << std::endl;
      // set-up neighborhood iterator
      typedef typename itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
      typename NeighborhoodIteratorType::RadiusType radius;
      radius.Fill(1);
      NeighborhoodIteratorType nit( radius, m_InputImage, m_InputImage->GetBufferedRegion() );
      unsigned int nlength = 1;
      for ( unsigned int d = 0; d < ImageDimension; d++ )
	{
	nlength *= nit.GetSize(d);
	}


      // set-up output iterator
      typedef typename itk::ImageRegionIterator< ImageType > OutputIteratorType;
      OutputIteratorType oit( m_OutputImage, m_OutputImage->GetBufferedRegion() );


      nit.GoToBegin();
      oit.GoToBegin();
      mit.GoToBegin();

      unsigned long progress = 0;

      while ( ! ( nit.IsAtEnd() || oit.IsAtEnd() || mit.IsAtEnd() )  )		// should be the same
	{
        if ( ++progress % 100000 == 0 )
          {
          std::cerr << "progress = " << progress << std::endl;
          }
        if ( mit.Value() != 0U )
          {


          // compute delta for the center voxel
          AccVectorType delta;
          delta.SetSize(nit.GetPixel(0).GetSize() );
          assert( nit.GetPixel(0).GetSize() == nlabels );
          delta.Fill(0.);

          VectorType vli =  nit.GetCenterPixel();

          unsigned int i = nlength / 2;

          double delta_count = 0.0;

          for ( unsigned int j = 0; j < nlength-1; j++ )
            {
            if ( j == i ) continue;
            VectorType vlj = nit.GetPixel(j);

            for ( unsigned int k = j+1; k < nlength; k++ )
              {
              if ( k == i ) continue;
              VectorType vlk = nit.GetPixel(k);

              delta_count += 1.0;

              // we're looking at P(l_i|l_j,l_k) where l_i is the label at voxel i
              //
              for (unsigned int lj = 0; lj < nlabels; lj++ )
                for (unsigned int lk = lj; lk < nlabels; lk++ )
                  {
                  for ( unsigned int li = 0; li < nlabels; li++ )
                    {
                    delta[li] += ( triplets[ li+nlabels*(lj+nlabels*lk) ] * vlj[lj] * vlk[lk]  );
                    }
                  }
              }
            }



          // update the output voxel with the delta
          VectorType vo = oit.Get();
          double sum = 0.0;
          for ( unsigned int li = 0; li < nlabels; li++ )
            {
            vo[li] = vli[li] * ( 1. + 1000000*delta[li] / delta_count );
//            vo[li] = delta[li] / delta_count ;
            sum += vo[li];
            }
          for ( unsigned int li = 0; li < nlabels; li++ )
            {
            vo[li] /= sum;
            }
          oit.Set(vo);
          }

	++nit;
	++oit;
        ++mit;
	}

      if ( m_Iterations > 0 )
        {
        std::cerr << "swap output image and input image" << std::endl;
        // swap m_OutputImage and m_InputImage
        typename ImageType::Pointer temp;
        temp = m_OutputImage;
        m_OutputImage = m_InputImage;
        m_InputImage = temp;
        }
      }

#endif // EXECUTE1

    }

  void Execute()
    {
    // allocate debug image
    m_DebugImage = LabelImageType::New();
    m_DebugImage->CopyInformation( m_InputImage );
    m_DebugImage->SetRegions( m_InputImage->GetBufferedRegion() );
    m_DebugImage->Allocate();
    m_DebugImage->FillBuffer(0U);

    // allocate output image
    m_OutputImage = ImageType::New();
    m_OutputImage->CopyInformation( m_InputImage );
    m_OutputImage->SetRegions( m_InputImage->GetBufferedRegion() );
    m_OutputImage->SetVectorLength( m_InputImage->GetVectorLength() );
    m_OutputImage->Allocate();

    // normalize the input image so each vector sums to one.
    typedef itk::ImageRegionIterator< ImageType > InputIteratorType;
    InputIteratorType it( m_InputImage, m_InputImage->GetBufferedRegion() );
    for (it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      VectorType v = it.Get();
      double sum  = 0.0;
      for ( unsigned int i = 0; i < v.GetSize(); i++ )
        sum += v[i];
      v /= sum;
      it.Set(v);                // probably not necessary the way VectorImage works, but whatever.
      }


    while ( m_Iterations-- > 0 )
      {
      std::cerr << "Iterations remaining: " << m_Iterations+1 << std::endl;
      // set-up neighborhood iterator
      typedef typename itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
      typename NeighborhoodIteratorType::RadiusType radius;
      radius.Fill(1);
      NeighborhoodIteratorType nit( radius, m_InputImage, m_InputImage->GetBufferedRegion() );
      unsigned int nlength = 1;
      for ( unsigned int d = 0; d < ImageDimension; d++ )
	{
	nlength *= nit.GetSize(d);
	}

      // set-up output iterator
      typedef typename itk::ImageRegionIterator< ImageType > OutputIteratorType;
      OutputIteratorType oit( m_OutputImage, m_OutputImage->GetBufferedRegion() );


      // set-up debug iterator
      typedef typename itk::ImageRegionIterator< LabelImageType > LabelIteratorType;
      LabelIteratorType dit( m_DebugImage, m_DebugImage->GetBufferedRegion() );

      nit.GoToBegin();
      dit.GoToBegin();
      oit.GoToBegin();


      VectorType cum;
      cum.SetSize(nit.GetPixel(0).GetSize());

#if 0
      unsigned long progress = 0;
#endif

      while ( ! ( nit.IsAtEnd() || oit.IsAtEnd() || dit.IsAtEnd() )  )		// should be the same
	{

#if 0
	if ( ++progress % 100000 == 0 )
	  {
	  std::cerr << "progress = " << progress << std::endl;
	  }
#endif

	const VectorType& c = nit.GetCenterPixel();

	// count up neighborhood statistics
	cum.Fill(0.);

	for ( unsigned int i = 0; i < nlength; i++ )
	  {
	  cum += nit.GetPixel(i);
	  }

	cum -= c;
	cum /= ( nlength - 1 );

        VectorType v = oit.Get();
        for ( unsigned int i = 0; i < v.Size(); i++ )
          v[i] = c[i];

#if 0
	if ( cum[4] > 0.2 && cum[5] > 0.2 || cum[4] > 0.2 && cum[0] > 0.2 )
	  {
          dit.Value() += 1U;
          v[7] = 0.001 * v[7];
          v[8] = 0.001 * v[8];
#else
	if ( cum[4] > 0.2 && cum[5] > 0.2 || cum[4] > 0.2 && cum[0] > 0.2 ||
              cum[8] > 0.2 && cum[5] > 0.2 )
	  {
          dit.Value() += 1U;
          double penalty = 0.1 * v[7];
          v[7] -= penalty;
          if ( cum[8] > 0.2 && cum[5] > 0.2 )
            {
            // basal ganglia gray matter
            v[8] += 0.5*penalty;
            v[5] += 0.5*penalty;
            }
          else
            {
            // cortical gray matter
            v[4] += 0.25*penalty;
            v[5] += 0.75*penalty;
            }
#endif

          // renormalize
          double sum = 0.0;
          for ( unsigned int i = 0; i < v.Size(); i++ )
            {
            sum += v[i];
            }
          v /= sum;
          }

	++nit;
	++dit;
	++oit;
	}

      if ( m_Iterations > 0 )
        {
        std::cerr << "swap output image and input image" << std::endl;
        // swap m_OutputImage and m_InputImage
        typename ImageType::Pointer temp;
        temp = m_OutputImage;
        m_OutputImage = m_InputImage;
        m_InputImage = temp;
        }
      }
    };

  void WriteFile(const std::string& outputfile, const std::string& debugfile = "")
    {

    if ( debugfile != "" )
      {
      std::cerr << "writing debug file " << debugfile << std::endl;
      typedef itk::ImageFileWriter<LabelImageType > LabelWriterType;
      typename LabelWriterType::Pointer writer = LabelWriterType::New();
      m_DebugImage->Print(std::cerr);
      writer->SetInput( m_DebugImage );
      writer->SetFileName( debugfile );
      writer->Update();
      }

    typedef itk::ImageFileWriter<ImageType > ImageWriterType;
    typename ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput( m_OutputImage );
    writer->SetFileName( outputfile );
    writer->Update();

#if 0
    if ( VectorDimension != 1 )
      {
      typedef itk::ImageFileWriter<ImageType > ImageWriterType;
      typename ImageWriterType::Pointer writer = ImageWriterType::New();
      writer->SetInput( m_OutputImage );
      writer->SetFileName( outputfile );
      writer->Update();
      }
    else
      {
      typedef itk::Image<ElementType, ImageDimension > ScalarImageType;
      typedef itk::ImageFileWriter< ScalarImageType > ImageWriterType;
      typedef itk::ImageRegionConstIterator< ImageType > VectorIteratorType;
      typedef itk::ImageRegionIterator<ScalarImageType > ScalarIteratorType;

      typename ScalarImageType::Pointer scalar = ScalarImageType::New();
      scalar->CopyInformation( m_OutputImage );
      scalar->SetRegions( m_OutputImage->GetLargestPossibleRegion() );
      scalar->Allocate();

      VectorIteratorType viter(m_OutputImage, m_OutputImage->GetBufferedRegion() );
      ScalarIteratorType siter(scalar, scalar->GetBufferedRegion() );

      viter.GoToBegin();
      siter.GoToBegin();
      while ( !viter.IsAtEnd() && !siter.IsAtEnd() )
	{
	siter.Set(viter.Value()[0]);
	++siter;
	++viter;
	}

      typename ImageWriterType::Pointer writer = ImageWriterType::New();
      writer->SetInput( scalar );
      writer->SetFileName(outputfile);
      writer->Update();
      }
#endif
    };

  void SetIterations(unsigned int& i) { m_Iterations = i; };

  void GetImageDimension() { return ImageDimension; };

  virtual ~RelaxationSchemeSpecial() {};

private:
  typename LabelImageType::Pointer m_ExampleImage;
  typename ImageType::Pointer m_InputImage;
  typename ImageType::Pointer m_OutputImage;
  typename LabelImageType::Pointer m_DebugImage;
  unsigned int m_Iterations;
};

RelaxationScheme* RelaxationScheme::ProbeFile(const std::string& inputfile)
  {
  typedef itk::ImageFileReader< itk::VectorImage<float> > ImageReaderType;
  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(inputfile);
  reader->GenerateOutputInformation();
  unsigned int imageDimension = reader->GetImageIO()->GetNumberOfDimensions();
  unsigned int vectorDimension = reader->GetImageIO()->GetNumberOfComponents();

  std::cout << "Image Dimension = " << imageDimension << ", Vector Size = " << vectorDimension << std::endl;

  if      ( imageDimension == 2)
    return new RelaxationSchemeSpecial<2>();
  else if ( imageDimension == 3)
    return new RelaxationSchemeSpecial<3>();
  else
    throw std::runtime_error("we only handle 2D or 3D images");

  return 0;	// NEVER REACHED
  }

};


#endif /*CRLRELAXATIONSCHEME_H_*/
