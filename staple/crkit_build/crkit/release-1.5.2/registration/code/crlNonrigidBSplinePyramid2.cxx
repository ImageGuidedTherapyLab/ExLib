#include "crlNonrigidBSplinePyramidClass.h"
#include "configuration.h"
#include <tclap/CmdLine.h>
#include "itkOrientImageFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkScaleSkewVersor3DTransform.h"


typedef itk::ScaleSkewVersor3DTransform< double > ScaleSkewVersor3DTransformType;

class Exception : public std::runtime_error
{
public:
Exception(std::string err) : std::runtime_error(err) {};
Exception(const char* err) : std::runtime_error(err) {};
};
  
ScaleSkewVersor3DTransformType::Pointer
readBulkTransform( std::string& file )
{
  std::cout << "reading bulk transform from file " << file << std::endl;

  typedef itk::TransformFileReader TransformReader;
  TransformReader::Pointer trsfReader = TransformReader::New();
  trsfReader->SetFileName(file);
  trsfReader->Update();

  //         itk::GenericExceptionMacro(std::string("failed reading transform from file: ")+file);

  typedef itk::TransformFileReader::TransformListType * TransformListType;
  TransformListType transforms = trsfReader->GetTransformList();
  std::cout << "Number of transforms = " << transforms->size() << std::endl;


  ScaleSkewVersor3DTransformType::Pointer trsf_read =
                   ScaleSkewVersor3DTransformType::New();

  if (transforms->size() == 1) 
    {
    TransformReader::TransformListType::const_iterator trsfit = transforms->begin();

    if (!strcmp((*trsfit)->GetNameOfClass(),"ScaleSkewVersor3DTransform")) 
      {
      trsf_read = static_cast<ScaleSkewVersor3DTransformType*>((*trsfit).GetPointer());
      std::cout << "Just loaded a scale skew versor 3d transform." << std::endl;
      trsf_read->Print(std::cout);
      } 
    else 
      {
      std::ostringstream sa;
      sa << "For Bulk transform file " << file << ", need a ScaleSkewVersor3DTransform instead of " << (*trsfit)->GetNameOfClass();
      throw Exception( sa.str() );
      }
    } 
  else 
    {
    std::ostringstream sa;
    sa << "Bulk transform file " << file << " has multiple( " << transforms->size() << ") transforms." <<
    	"Can't initialize.";
    throw Exception(sa.str());
    }
  
  return trsf_read;
}
      
      
      

int main( int argc, char* argv[] )
{
  try 
    {
    // TCLAP command line parsing
    TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',CRKIT_VERSION_STRING );

    TCLAP::UnlabeledValueArg<std::string> fixedImageFileArg("fixedImage", "Filename of the fixed (target) image.", true, "", "fixed image", cmd);
    TCLAP::UnlabeledValueArg<std::string> movingImageFileArg("movingImage", "Filename of the moving image.", true, "", "moving image", cmd);
    TCLAP::UnlabeledValueArg<std::string> outputImageFileArg("outputImage", "Filename of the output image.", true, "", "output image", cmd);
    TCLAP::UnlabeledValueArg<std::string> outputTransformFileArg("outputTransform","Filename transform (spline coefficients) will be written to.", true, "", "output transform", cmd);
    TCLAP::UnlabeledValueArg<std::string> outputDefFileArg("outputDeform","Filename to write a deformation field to.", true, "", "output deform field", cmd);

    TCLAP::ValueArg<unsigned int> pyramidLevelArg("p","pyramidLevels","The number of levels of a multi-resolution pyramid to run", false, 2, "pyramid levels", cmd);
    TCLAP::ValueArg<unsigned int> binsArg("b","bins","The number of bins to use when histogramming for the MI metric.", false, 20, "bins", cmd);
    TCLAP::ValueArg<unsigned int> samplesArg("s","samples","The number points to sample randomly when computing the metric.  Zero means 'all'", 
		false, 0, "samples", cmd);
    TCLAP::ValueArg<unsigned int> iterationsArg("i","iters","The maximum number of iterations of the Powell optimizer per pyramid level.", 
		false, 200, "iterations", cmd);
    TCLAP::ValueArg<unsigned int> controlPointsArg("c","nctrlpts","Number of B-spline control points to use in each dimension", false, 5, "control points", cmd);
    TCLAP::ValueArg<float> maxStepLengthArg("u","maxStepLength","Optimizer maximum (upper) step length.", false, 10, "max steplength", cmd);
    TCLAP::ValueArg<float> minStepLengthArg("l","minStepLength","Optimizer minimum (lower) step length.", false, 0.01, "min steplength", cmd);
    TCLAP::ValueArg<float> relaxFactorArg("r","relaxationFactor","Optimizer relaxation factor.", false, 0.7, "relaxation factor", cmd);

    TCLAP::ValueArg<std::string> optInitTransformFileArg("t","initTrans",
	"Filename specifying an initial ScaleSkewVersor (affine) bulk transform to read (e.g. from crlAffineMIRegistration)",
	false, "", "initial affine transform", cmd);
    TCLAP::ValueArg<std::string> optFixedMaskFileArg("m","fixedmask",
	"Image specifying a mask of the fixed image domain.  Locations where image is precisely zero will NOT be sampled for the metric", false, "", "fixed mask image", cmd);
    TCLAP::SwitchArg optFixedMaskFlagArg("M","fixedmaskflag",
	"Flag specifying that the fixed image has already been masked.  Locations where _it_ is precisely zero will NOT be sampled for the metric",
	 cmd);
    TCLAP::SwitchArg optNormalizedMI("n", "BROKENnmiflag", "BROKENFlag specifying that normalized mutual information should be used rather than Mattes' MI.  Default is Mattes' MI.",
	cmd);
    cmd.parse(argc,argv);
    
    
    if (optFixedMaskFlagArg.isSet() && optFixedMaskFileArg.isSet() )
      throw Exception("you can't set -m AND -M");
    
    // read and orient all inputs
    static const unsigned int ImageDimension = 3;
    typedef float PixelType;
    typedef itk::OrientedImage< PixelType, ImageDimension > ImageType;
    typedef crl::NonrigidBSplinePyramid< ImageType, ImageDimension > RegistrationType;
    
    typedef itk::ImageFileReader<ImageType> ImageReaderType;
    typedef itk::ImageFileWriter<ImageType> ImageWriterType;
    
    ImageReaderType::Pointer fixedReader = ImageReaderType::New();
    ImageReaderType::Pointer movingReader = ImageReaderType::New();
    
    fixedReader->SetFileName(fixedImageFileArg.getValue());
    movingReader->SetFileName(movingImageFileArg.getValue());
    
    typedef itk::OrientImageFilter<ImageType,ImageType> OrientType;
    OrientType::Pointer fixedOrient = OrientType::New();
    OrientType::Pointer movingOrient = OrientType::New();
    fixedOrient->SetInput( fixedReader->GetOutput() );
    movingOrient->SetInput( movingReader->GetOutput() );
    fixedOrient->UseImageDirectionOn();
    movingOrient->UseImageDirectionOn();
    fixedOrient->SetDesiredCoordinateOrientationToAxial();
    movingOrient->SetDesiredCoordinateOrientationToAxial();
    
    fixedOrient->Update();	// errors will be caught below
    movingOrient->Update();	
    
    
    
    // set up the process and run
    typedef crl::NonrigidBSplinePyramid<ImageType,3> RegType;
    RegType::Pointer registration = RegType::New();
    registration->SetFixedImage(fixedOrient->GetOutput());
    registration->SetMovingImage(movingOrient->GetOutput());
    if (optNormalizedMI.isSet()) registration->NormalizedMIOn();
    if (pyramidLevelArg.isSet()) registration->SetNumberPyramidLevels(pyramidLevelArg.getValue());
    if (controlPointsArg.isSet()) registration->SetControlPoints(controlPointsArg.getValue());
    if (maxStepLengthArg.isSet()) registration->SetMaxStepLength(maxStepLengthArg.getValue());
    if (minStepLengthArg.isSet()) registration->SetMinStepLength(minStepLengthArg.getValue());
    if (relaxFactorArg.isSet()) registration->SetRelaxationFactor(relaxFactorArg.getValue());
    if (binsArg.isSet()) registration->SetNumberBins(binsArg.getValue());
    if (samplesArg.isSet()) registration->SetNumberSamples(samplesArg.getValue());
    if (iterationsArg.isSet()) registration->SetNumberIterations(iterationsArg.getValue());
    
    // optionally read in an Affine Transform
    if ( optInitTransformFileArg.isSet() )
      {
      registration->SetInitialBulkTransform( readBulkTransform(optInitTransformFileArg.getValue()));
      }

    
    // optionally set a mask
    RegistrationType::ImageMaskObjectType::Pointer mask = 0;
    typedef RegistrationType::ImageMaskObjectType::ImageType MaskImageType;
    
    typedef itk::CastImageFilter<ImageType,MaskImageType> MaskCasterType;
    MaskCasterType::Pointer maskCaster = MaskCasterType::New();
    
    if ( optFixedMaskFlagArg.isSet() ) // -M
      {
      maskCaster->SetInput( fixedOrient->GetOutput() );
      maskCaster->Update();
      mask = RegistrationType::ImageMaskObjectType::New();
      mask->SetImage( maskCaster->GetOutput() ); // the image *is* the mask with -M
      }
    else if ( optFixedMaskFileArg.isSet() ) // -m file
      {
      ImageReaderType::Pointer maskreader = ImageReaderType::New();
      maskreader->SetFileName( optFixedMaskFileArg.getValue() );
      OrientType::Pointer maskOrient = OrientType::New();
      maskOrient->SetInput( maskreader->GetOutput() );
      maskOrient->UseImageDirectionOn();
      maskOrient->SetDesiredCoordinateOrientationToAxial();
      maskCaster->SetInput( maskOrient->GetOutput() );
      maskCaster->Update();
      
      mask = RegistrationType::ImageMaskObjectType::New();
      mask->SetImage( maskCaster->GetOutput() );
      }
    if ( mask.IsNotNull() )
      {
      std::cerr << "setting a mask for the fixed image" << std::endl;
      registration->SetMaskObject( mask );
      }
    

    //
    // RUN THE REGISTRATION
    //
    std::cerr << "running the filter" << std::endl;
    registration->Update();
    std::cerr << "done! running the filter" << std::endl;
    //////////////////////////////////////////////////////////////
  
    // get the results
    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetInput(registration->GetOutput());
    writer->SetFileName(outputImageFileArg.getValue());
    writer->UseCompressionOn();
    writer->Update();
    

    typedef itk::ImageFileWriter< RegType::DeformationFieldType > DefWriterType;
    DefWriterType::Pointer defwriter = DefWriterType::New();
    defwriter->SetInput(registration->GetDeformationField());
    defwriter->SetFileName(outputDefFileArg.getValue() );
    defwriter->UseCompressionOn();
    defwriter->Update();

    
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer transwriter = TransformWriterType::New();
    transwriter->AddTransform(registration->GetFinalTransform());
    transwriter->SetFileName( outputTransformFileArg.getValue() );
    transwriter->Update();
    
    
    } 
  catch (TCLAP::ArgException& e)
    {
    std::cerr << "error parsing " << e.argId() << ": " << e.error() << std::endl;
    exit(1);
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "itk error: " << e << std::endl;
    exit(1);
    }
  catch (Exception& e)
    {
    std::cerr << "Error: " << e.what() << std::endl;
    exit(1);
    }
  
  exit(0);
}

  

