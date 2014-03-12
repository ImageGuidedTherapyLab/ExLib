/*
 * Copyright 2008-2009 Children's Hospital Boston
 * Contact: Simon Warfield simon.warfield@childrens.harvard.edu
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#include "crlZeroImage.h"

int main(int argc, char *argv[])
{

  typedef   float PixelType;

  if (argc < 3) {
    std::cout <<
      "Usage: " << argv[0] << " outputFileName dimensionality " <<
      "startx starty startz " <<
      "numx numy numz " << 
      "spacingx spacingy spacingz " <<
      "originx originy originz [value]"
      << std::endl;
    return EXIT_FAILURE;
  }

  crlZeroImageBase::Pointer crlZeroImagePtr = 0;

  signed int newDimensionality = atoi(argv[2]);

  if ((argc != (3 + newDimensionality*4)) && (argc != (4 + newDimensionality*4))) {
    std::cout <<
      "Usage: " << argv[0] << " outputFileName dimensionality " <<
      "startx starty startz " <<
      "numx numy numz " << 
      "spacingx spacingy spacingz " <<
      "originx originy originz [value]"
      << std::endl;
    return EXIT_FAILURE;
  }

  switch( newDimensionality ) {
    case 1 : crlZeroImagePtr = crlZeroImage<1, PixelType>::New();
             break;
    case 2 : crlZeroImagePtr = crlZeroImage<2, PixelType>::New();
             break;
    case 3 : crlZeroImagePtr = crlZeroImage<3, PixelType>::New();
             break;
    case 4 : crlZeroImagePtr = crlZeroImage<4, PixelType>::New();
             break;
    case 5 : crlZeroImagePtr = crlZeroImage<5, PixelType>::New();
             break;
    default : std::cerr << "No template specialization for dimensionality "
                    << argv[2] << std::endl;
  };

  if (!crlZeroImagePtr) {
    std::cerr << "No template specialization for dimensionality "
                    << argv[2] << std::endl;
    return EXIT_FAILURE;
  }

  crlZeroImagePtr->SetOutputFileName( argv[1] );
  if (argc == (4 + newDimensionality*4))
    crlZeroImagePtr->SetValue(atof(argv[3+4*newDimensionality]));

  for (signed int i = 0; i < newDimensionality; i++) {
    crlZeroImagePtr->SetStartIndex( i, atoi(argv[i+3]) );
    crlZeroImagePtr->SetSize( i, atoi(argv[i+3+newDimensionality]) );
    crlZeroImagePtr->SetSpacing( i, atof(argv[i+3+2*newDimensionality]) );
    crlZeroImagePtr->SetOrigin( i, atof(argv[i+3+3*newDimensionality]) );
  }

  int ret = crlZeroImagePtr->Execute();

  if (ret != EXIT_SUCCESS)
  {
    std::cerr << argv[0] << "Returned error code " << ret << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
