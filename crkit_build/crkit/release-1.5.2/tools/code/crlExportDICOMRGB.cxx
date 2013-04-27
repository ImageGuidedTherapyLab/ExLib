#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageSeriesWriter.h>
#include <itkImageSeriesReader.h>
#include <itkNumericSeriesFileNames.h>
#include <itkMetaDataObject.h>
#include <itkRGBPixel.h>

#include <tclap/CmdLine.h>

#include "configuration.h"


typedef unsigned char PixelType;
const unsigned int ImageDimension = 3;
const unsigned int OutputDimension = 2;

typedef itk::OrientedImage<itk::RGBPixel<PixelType>,ImageDimension> ImageType;
typedef itk::ImageFileReader<ImageType> ReaderType;

typedef itk::OrientedImage<itk::RGBPixel<PixelType>,OutputDimension> Image2DType;
typedef itk::ImageSeriesWriter<ImageType,Image2DType> SeriesWriterType;

typedef itk::ImageSeriesReader<ImageType> TemplateReaderType;

itk::GDCMImageIO::Pointer gdcmIO = itk::GDCMImageIO::New();
itk::GDCMImageIO::Pointer templateIO = itk::GDCMImageIO::New();

int main( int argc, char *argv[] )
{
  std::string inputFile;
  std::string outputPrefix;
  std::string dataTemplate = "";
  std::string seriesDesc = "PROCESSED DATA";
  std::string mrnPrefix = "";

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ',CRKIT_VERSION_STRING );
      TCLAP::ValueArg<std::string> seriesDescArg("d","description","Series Description",false,"","series description",cmd);
      TCLAP::UnlabeledValueArg<std::string> inputFileArg("inputFile","Input File Name",true,"","input file name",cmd);
      TCLAP::UnlabeledValueArg<std::string> outputPrefixArg("outputPrefix","Output Prefix",true,"","output prefix",cmd);
      TCLAP::UnlabeledValueArg<std::string> dataTemplateArg("dataTemplate","DICOM Data Template Directory",false,"","data template",cmd);
      TCLAP::ValueArg<std::string> mrnPrefixArg("p","prefix","MRN prefix",false,"","MRN prefix",cmd);

      cmd.parse(argc,argv);

      if (inputFileArg.isSet()) inputFile = inputFileArg.getValue();
      if (outputPrefixArg.isSet()) outputPrefix = outputPrefixArg.getValue();
      if (dataTemplateArg.isSet()) dataTemplate = dataTemplateArg.getValue();
      if (seriesDescArg.isSet()) seriesDesc = seriesDescArg.getValue();
      if (mrnPrefixArg.isSet()) mrnPrefix = mrnPrefixArg.getValue();
    }
  catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(EXIT_FAILURE);
    }

  ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName(inputFile);
  try
    {
      reader->Update();
    }
  catch (itk::ExceptionObject &e)
    {
      std::cerr << "Could not read input file: " << e << std::endl;
      return EXIT_FAILURE;
    }

  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();

  seriesWriter->SetInput(reader->GetOutput());
  seriesWriter->SetImageIO(gdcmIO);

  TemplateReaderType::Pointer templateReader = TemplateReaderType::New();
  itk::GDCMSeriesFileNames::Pointer templateNames = itk::GDCMSeriesFileNames::New();

  if (dataTemplate != "")
    {
      templateReader->SetImageIO(templateIO);
      templateNames->SetInputDirectory(dataTemplate);
      templateReader->SetFileNames(templateNames->GetInputFileNames());
      try
	{
	  templateReader->Update();
	}
      catch (itk::ExceptionObject &e)
	{
	  std::cerr << "Could not read DICOM data template: " << e << std::endl;
	  return EXIT_FAILURE;
	}
    }

  itk::MetaDataDictionary &dict = gdcmIO->GetMetaDataDictionary();

  if (dataTemplate != "")
    {
      itk::MetaDataDictionary &templatedict = templateIO->GetMetaDataDictionary();
      std::string tmp;
      // patient name
      itk::ExposeMetaData(templatedict,"0010|0010",tmp);
      itk::EncapsulateMetaData<std::string>(dict,"0010|0010",tmp);
      // patient ID
      itk::ExposeMetaData(templatedict,"0010|0020",tmp);
      itk::EncapsulateMetaData<std::string>(dict,"0010|0020", mrnPrefix+tmp);
      // patient's DOB
      itk::ExposeMetaData(templatedict,"0010|0030",tmp);
      itk::EncapsulateMetaData<std::string>(dict,"0010|0030",tmp);
      // patient's sex
      itk::ExposeMetaData(templatedict,"0010|0040",tmp);
      itk::EncapsulateMetaData<std::string>(dict,"0010|0040",tmp);
      // patient's weight
      itk::ExposeMetaData(templatedict,"0010|1030",tmp);
      itk::EncapsulateMetaData<std::string>(dict,"0010|1030",tmp);
      // ethnic group
      itk::ExposeMetaData(templatedict,"0010|2160",tmp);
      itk::EncapsulateMetaData<std::string>(dict,"0010|2160",tmp);
      // referring physician
      itk::ExposeMetaData(templatedict,"0008|0090",tmp);
      itk::EncapsulateMetaData<std::string>(dict,"0008|0090",tmp);
      // performing physician
      itk::ExposeMetaData(templatedict,"0008|1050",tmp);
      itk::EncapsulateMetaData<std::string>(dict,"0008|1050",tmp);
      // physician of record
      itk::ExposeMetaData(templatedict,"0008|1048",tmp);
      itk::EncapsulateMetaData<std::string>(dict,"0008|1048",tmp);
      // accession number
      itk::ExposeMetaData(templatedict,"0008|0050",tmp);
      itk::EncapsulateMetaData<std::string>(dict,"0008|0050",tmp);

    }
  // modality
  itk::EncapsulateMetaData<std::string>(dict,"0008|0060","MR");
  // image type
  itk::EncapsulateMetaData<std::string>(dict,"0008|0008","DERIVED\\SECONDARY");
  // manufacturer
  itk::EncapsulateMetaData<std::string>(dict,"0008|0070","CRL");
  // institution name
  itk::EncapsulateMetaData<std::string>(dict,"0008|0080","Children's Hospital Boston");
  // institution address
  itk::EncapsulateMetaData<std::string>(dict,"0008|0081","300 Longwood Ave, Boston, MA 02115");
  // conversion type
  itk::EncapsulateMetaData<std::string>(dict,"0008|0064","DV");
  // study description
  itk::EncapsulateMetaData<std::string>(dict,"0008|1030","Computational Radiology Lab");
  // series description
  itk::EncapsulateMetaData<std::string>(dict,"0008|103e",seriesDesc);
  // sequence name
  itk::EncapsulateMetaData<std::string>(dict,"0018|0024","COMPUTATIONAL");
  // patient position
  itk::EncapsulateMetaData<std::string>(dict,"0018|5100","HFS");
  // slice thickness
  std::stringstream ss;
  ss << reader->GetOutput()->GetSpacing()[2];
  itk::EncapsulateMetaData<std::string>(dict,"0018|0050",ss.str());

  itk::NumericSeriesFileNames::Pointer names = itk::NumericSeriesFileNames::New();
  outputPrefix = outputPrefix+"im%04d.dcm";
  names->SetSeriesFormat(outputPrefix.c_str());
  names->SetStartIndex(reader->GetOutput()->GetLargestPossibleRegion().GetIndex()[2]);
  names->SetEndIndex(reader->GetOutput()->GetLargestPossibleRegion().GetIndex()[2]+reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2]-1);

  seriesWriter->SetFileNames(names->GetFileNames());

  try
    {
      seriesWriter->Write();
    }
  catch (itk::ExceptionObject &exp)
    {
      std::cout << "Could not write DICOM series." << std::endl;
      std::cout << exp << std::endl;
      return EXIT_FAILURE;
    }
}
