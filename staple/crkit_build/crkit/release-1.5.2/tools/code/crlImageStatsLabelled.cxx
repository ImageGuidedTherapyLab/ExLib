
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkLabelStatisticsImageFilter.h>

#include <iostream>
// Enable control of printed number precision
#include <iomanip>


int main(int argc, char *argv[])
{

  static int const ImageDimension = 3;
  typedef   float PixelType;

  typedef   itk::Image<PixelType,ImageDimension>  ImageType;
  typedef   itk::Image<unsigned int,ImageDimension>  LabelImageType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;
  typedef   itk::ImageFileReader< LabelImageType >    LabelImageReaderType;
  typedef   itk::LabelStatisticsImageFilter< ImageType , LabelImageType > FilterType;

  float maxlabel = 10.0;

  if ( (argc != 4) ) {
    std::cout << 
      "Usage: crlImageStatsLabelled inimage labelimage maxlabel" << std::endl;
    exit(1);
  }

  maxlabel = ::atof(argv[3]);
  std::cout << "maxlabel is " << maxlabel << std::endl;

  ImageReaderType::Pointer r;
  LabelImageReaderType::Pointer labelReader = LabelImageReaderType::New();
  FilterType::Pointer f;

  r = ImageReaderType::New();
  f = FilterType::New();

  r->SetFileName(argv[1]);
  labelReader->SetFileName(argv[2]);

  try {
    r->Update();
    labelReader->Update();
  } catch (itk::ExceptionObject &e) {
      std::cerr << "Caught ITK exception: " << e << std::endl;
      exit(1);
  }

  f->SetInput( r->GetOutput() );
  f->SetLabelInput( labelReader->GetOutput() );

  f->UseHistogramsOff();

  // Now run the filter
  f->Update();

  std::cout << setiosflags(std::ios::fixed) << std::setprecision(15);

  std::cout << "Intensity properties; " << std::endl;
  std::cout << "SPREADSHEET, Label";
  for (unsigned int i = 0; i <= maxlabel; i++) {
    std::cout << ", " << i;
  }
  std::cout << std::endl;
  std::cout << "SPREADSHEET, Mean";
  for (unsigned int i = 0; i <= maxlabel; i++) {
    std::cout << ", " << f->GetMean(i);
  }
  std::cout << std::endl;
  std::cout << "SPREADSHEET, Count";
  for (unsigned int i = 0; i <= maxlabel; i++) {
    std::cout << ", " << f->GetCount(i);
  }
  std::cout << std::endl;
  std::cout << "SPREADSHEET, Minimum";
  for (unsigned int i = 0; i <= maxlabel; i++) {
    std::cout << ", " << f->GetMinimum(i);
  }
  std::cout << std::endl;
  std::cout << "SPREADSHEET, Maximum";
  for (unsigned int i = 0; i <= maxlabel; i++) {
    std::cout << ", " << f->GetMaximum(i);
  }
  std::cout << std::endl;
  std::cout << "SPREADSHEET, Variance";
  for (unsigned int i = 0; i <= maxlabel; i++) {
    std::cout << ", " << f->GetVariance(i);
  }
  std::cout << std::endl;
  std::cout << "SPREADSHEET, Sum";
  for (unsigned int i = 0; i <= maxlabel; i++) {
    std::cout << ", " << f->GetSum(i);
  }
  std::cout << std::endl;

  ImageType::RegionType inputRegion =
                   r->GetOutput()->GetLargestPossibleRegion();
  ImageType::SizeType size = inputRegion.GetSize();
  ImageType::IndexType start = inputRegion.GetIndex();
  ImageType::SpacingType spacing = r->GetOutput()->GetSpacing();

  std::cout << "Geometric properties; " << std::endl;
  std::cout << "Size, Start, Spacing, Origin  " << std::endl;
  std::cout << size << ", " ;
  std::cout << start << ", " ;
  std::cout << r->GetOutput()->GetSpacing() << ", " ;
  std::cout << r->GetOutput()->GetOrigin()  << "  " << std::endl;

  exit(0); // success
}

