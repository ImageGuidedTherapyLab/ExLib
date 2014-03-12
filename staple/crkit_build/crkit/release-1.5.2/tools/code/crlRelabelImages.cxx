/*
 * Copyright 2008 Simon Warfield simon.warfield@childrens.harvard.edu
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://www.crl.med.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#include "crlRelabelImages.h"

int main( int argc, char ** argv )
{
  static int const ImageDimension = 3;
  typedef   float PixelType;

  typedef   itk::Image<PixelType,ImageDimension>  ImageType;
  typedef   itk::ImageFileReader< ImageType >    ImageReaderType;

  if ( (argc != 6) && (argc != 7) )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " inputImageFile1 inputImageFile2 "
              << " inputlabels outputlabels "
              << " outputImageFile"
              << " ?nomatchoutputvalue?"
              << std::endl;
      std::cerr << "change the values of inputImage1 that match the inputlabels to the corresponding labels of outputlabels and use the values of inputImage2 when there are no matches in inputlabels UNLESS a default nomatch value is supplied in which case use that." << std::endl;
    return EXIT_FAILURE;
    }

  crl::RelabelImagesBase::Pointer crlRelabelImagesPtr = 0;
  ImageType::Pointer image1;
  image1 = ImageType::New();

  ImageReaderType::Pointer image1reader;
  image1reader = ImageReaderType::New();
  image1reader->SetFileName( argv[1] );

  try {
    image1reader->GenerateOutputInformation();
  } catch ( itk::ExceptionObject &err )
  {
    std::cout << "ExceptionObject caught!" << std::endl;
    std::cout << err << std::endl;
    return -1;
  }

  itk::ImageIOBase *imageIO = image1reader->GetImageIO();
  if (imageIO->GetPixelType() == itk::ImageIOBase::SCALAR) {
     switch(imageIO->GetComponentType()) {
       case itk::ImageIOBase::UCHAR :
         if (imageIO->GetNumberOfDimensions() == 1) {
crlRelabelImagesPtr = crl::RelabelImages< 1, unsigned char >::New();
         } else if (imageIO->GetNumberOfDimensions() == 2) {
crlRelabelImagesPtr = crl::RelabelImages< 2, unsigned char >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlRelabelImagesPtr = crl::RelabelImages< 3, unsigned char >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlRelabelImagesPtr = crl::RelabelImages< 4, unsigned char >::New();
         } else if (imageIO->GetNumberOfDimensions() == 5) {
crlRelabelImagesPtr = crl::RelabelImages< 5, unsigned char >::New();
         } else if (imageIO->GetNumberOfDimensions() == 6) {
crlRelabelImagesPtr = crl::RelabelImages< 6, unsigned char >::New();
         }
         break;
       case itk::ImageIOBase::CHAR :
         if (imageIO->GetNumberOfDimensions() == 1) {
crlRelabelImagesPtr = crl::RelabelImages< 1, signed char >::New();
         } else if (imageIO->GetNumberOfDimensions() == 2) {
crlRelabelImagesPtr = crl::RelabelImages< 2, signed char >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlRelabelImagesPtr = crl::RelabelImages< 3, signed char >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlRelabelImagesPtr = crl::RelabelImages< 4, signed char >::New();
         } else if (imageIO->GetNumberOfDimensions() == 5) {
crlRelabelImagesPtr = crl::RelabelImages< 5, signed char >::New();
         } else if (imageIO->GetNumberOfDimensions() == 6) {
crlRelabelImagesPtr = crl::RelabelImages< 6, signed char >::New();
         }
         break;
       case itk::ImageIOBase::USHORT :
         if (imageIO->GetNumberOfDimensions() == 1) {
crlRelabelImagesPtr = crl::RelabelImages< 1, unsigned short >::New();
         } else if (imageIO->GetNumberOfDimensions() == 2) {
crlRelabelImagesPtr = crl::RelabelImages< 2, unsigned short >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlRelabelImagesPtr = crl::RelabelImages< 3, unsigned short >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlRelabelImagesPtr = crl::RelabelImages< 4, unsigned short >::New();
         } else if (imageIO->GetNumberOfDimensions() == 5) {
crlRelabelImagesPtr = crl::RelabelImages< 5, unsigned short >::New();
         } else if (imageIO->GetNumberOfDimensions() == 6) {
crlRelabelImagesPtr = crl::RelabelImages< 6, unsigned short >::New();
         }
         break;
       case itk::ImageIOBase::SHORT :
         if (imageIO->GetNumberOfDimensions() == 1) {
crlRelabelImagesPtr = crl::RelabelImages< 1, signed short >::New();
         } else if (imageIO->GetNumberOfDimensions() == 2) {
crlRelabelImagesPtr = crl::RelabelImages< 2, signed short >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlRelabelImagesPtr = crl::RelabelImages< 3, signed short >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlRelabelImagesPtr = crl::RelabelImages< 4, signed short >::New();
         } else if (imageIO->GetNumberOfDimensions() == 5) {
crlRelabelImagesPtr = crl::RelabelImages< 5, signed short >::New();
         } else if (imageIO->GetNumberOfDimensions() == 6) {
crlRelabelImagesPtr = crl::RelabelImages< 6, signed short >::New();
         }
         break;
       case itk::ImageIOBase::UINT :
         if (imageIO->GetNumberOfDimensions() == 1) {
crlRelabelImagesPtr = crl::RelabelImages< 1, unsigned int >::New();
         } else if (imageIO->GetNumberOfDimensions() == 2) {
crlRelabelImagesPtr = crl::RelabelImages< 2, unsigned int >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlRelabelImagesPtr = crl::RelabelImages< 3, unsigned int >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlRelabelImagesPtr = crl::RelabelImages< 4, unsigned int >::New();
         } else if (imageIO->GetNumberOfDimensions() == 5) {
crlRelabelImagesPtr = crl::RelabelImages< 5, unsigned int >::New();
         } else if (imageIO->GetNumberOfDimensions() == 6) {
crlRelabelImagesPtr = crl::RelabelImages< 6, unsigned int >::New();
         }
         break;
       case itk::ImageIOBase::INT :
         if (imageIO->GetNumberOfDimensions() == 1) {
crlRelabelImagesPtr = crl::RelabelImages< 1, signed int >::New();
         } else if (imageIO->GetNumberOfDimensions() == 2) {
crlRelabelImagesPtr = crl::RelabelImages< 2, signed int >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlRelabelImagesPtr = crl::RelabelImages< 3, signed int >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlRelabelImagesPtr = crl::RelabelImages< 4, signed int >::New();
         } else if (imageIO->GetNumberOfDimensions() == 5) {
crlRelabelImagesPtr = crl::RelabelImages< 5, signed int >::New();
         } else if (imageIO->GetNumberOfDimensions() == 6) {
crlRelabelImagesPtr = crl::RelabelImages< 6, signed int >::New();
         }
         break;
       case itk::ImageIOBase::ULONG :
         if (imageIO->GetNumberOfDimensions() == 1) {
crlRelabelImagesPtr = crl::RelabelImages< 1, unsigned long >::New();
         } else if (imageIO->GetNumberOfDimensions() == 2) {
crlRelabelImagesPtr = crl::RelabelImages< 2, unsigned long >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlRelabelImagesPtr = crl::RelabelImages< 3, unsigned long >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlRelabelImagesPtr = crl::RelabelImages< 4, unsigned long >::New();
         } else if (imageIO->GetNumberOfDimensions() == 5) {
crlRelabelImagesPtr = crl::RelabelImages< 5, unsigned long >::New();
         } else if (imageIO->GetNumberOfDimensions() == 6) {
crlRelabelImagesPtr = crl::RelabelImages< 6, unsigned long >::New();
         }
         break;
       case itk::ImageIOBase::LONG :
         if (imageIO->GetNumberOfDimensions() == 1) {
crlRelabelImagesPtr = crl::RelabelImages< 1, signed long >::New();
         } else if (imageIO->GetNumberOfDimensions() == 2) {
crlRelabelImagesPtr = crl::RelabelImages< 2, signed long >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlRelabelImagesPtr = crl::RelabelImages< 3, signed long >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlRelabelImagesPtr = crl::RelabelImages< 4, signed long >::New();
         } else if (imageIO->GetNumberOfDimensions() == 5) {
crlRelabelImagesPtr = crl::RelabelImages< 5, signed long >::New();
         } else if (imageIO->GetNumberOfDimensions() == 6) {
crlRelabelImagesPtr = crl::RelabelImages< 6, signed long >::New();
         }
         break;
       case itk::ImageIOBase::FLOAT :
         if (imageIO->GetNumberOfDimensions() == 1) {
crlRelabelImagesPtr = crl::RelabelImages< 1, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 2) {
crlRelabelImagesPtr = crl::RelabelImages< 2, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlRelabelImagesPtr = crl::RelabelImages< 3, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlRelabelImagesPtr = crl::RelabelImages< 4, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 5) {
crlRelabelImagesPtr = crl::RelabelImages< 5, float >::New();
         } else if (imageIO->GetNumberOfDimensions() == 6) {
crlRelabelImagesPtr = crl::RelabelImages< 6, float >::New();
         }
         break;
       case itk::ImageIOBase::DOUBLE :
         if (imageIO->GetNumberOfDimensions() == 1) {
crlRelabelImagesPtr = crl::RelabelImages< 1, double >::New();
         } else if (imageIO->GetNumberOfDimensions() == 2) {
crlRelabelImagesPtr = crl::RelabelImages< 2, double >::New();
         } else if (imageIO->GetNumberOfDimensions() == 3) {
crlRelabelImagesPtr = crl::RelabelImages< 3, double >::New();
         } else if (imageIO->GetNumberOfDimensions() == 4) {
crlRelabelImagesPtr = crl::RelabelImages< 4, double >::New();
         } else if (imageIO->GetNumberOfDimensions() == 5) {
crlRelabelImagesPtr = crl::RelabelImages< 5, double >::New();
         } else if (imageIO->GetNumberOfDimensions() == 6) {
crlRelabelImagesPtr = crl::RelabelImages< 6, double >::New();
         }
         break;
      default : break;
     }
   }

  // In case the switch statement doesn't initialize the pointer
  if (!crlRelabelImagesPtr) {
    std::cout << "crlRelabelImages not implemented for pixel type " <<
        imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) <<
         " with dimensions " << imageIO->GetNumberOfDimensions() <<
        std::endl;
    return EXIT_FAILURE;
  }

  crlRelabelImagesPtr->SetImage1FileName( argv[1] );
  crlRelabelImagesPtr->SetImage2FileName( argv[2] );
  crlRelabelImagesPtr->SetInputLabels( argv[3] );
  crlRelabelImagesPtr->SetOutputLabels( argv[4] );
  crlRelabelImagesPtr->SetOutputImageFileName( argv[5] );
  if (argc == 7) {
    crlRelabelImagesPtr->SetDefaultOutputValue( std::atof(argv[6]) );
  }

  int ret = crlRelabelImagesPtr->Execute();
  if (ret != 0) {
    std::cerr << "crlRelabelImages returned an error code " << ret <<
    std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}
