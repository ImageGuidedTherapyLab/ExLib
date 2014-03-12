
#include <itkOrientedImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <iostream>
#include <string>
#include <vector>
#include "vtkPolyDataWriter.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkZLibDataCompressor.h"
#include "vtkPolyDataReader.h"
#include "vtkCleanPolyData.h"

#include "crlVtkMeshIO.h"

#include "tclap/CmdLine.h"
#include "configuration.h"
#include <limits.h>

#ifdef WIN32
#define snprintf _snprintf
#define PATH_MAX 512
#endif

int main(int argc, char* argv[]){

  static int const ImageDimension = 3;

  // We will read in a VTK file representing tracts,
  // a set of labels each tract must touch,
  // a set of labels each tract must not touch,
  // and a minimum size each tract must have.

  vtkIdType minimumTractSize = 0;
  std::vector<signed int> touchlabels;
  std::vector<signed int> donttouchlabels;
  std::string inputTractFileName;
  std::string roiImageFileName;
  std::string outputTractFileName;

  // Exceptions are thrown by TCLAP for bad arguments in the CLI.
  try {
    // Define the command line object, and insert a message
    // that describes the program.
    TCLAP::CmdLine cmd("Computational Radiology Laboratory", ' ', 
        CRKIT_VERSION_STRING );

    // Define a value argument for the input image on the command line.
    // This can be specified without any flag. It is simply taken from the 
    // first unflagged position.
    TCLAP::UnlabeledValueArg<std::string> inputImageArg("input", 
            "Input tract file name (.vtk, .vtp or .gii)", true, 
            "", 
            "Input tract file name",
            cmd); 

    TCLAP::UnlabeledValueArg<std::string> roiImageArg("roiFileName", 
            "ROI image file name", true, 
            "", 
            "ROI image file name",
            cmd); 

      TCLAP::ValueArg<std::string> doTouchArg( "t", "touch", 
    "A string of integer labels that each tract must touch", false, "", 
    "String of integer label values",cmd);

      TCLAP::ValueArg<std::string> dontTouchArg( "d", "donttouch", 
    "A string of integer labels that each tract must not touch", false, "", 
    "String of integer label values",cmd);

      TCLAP::ValueArg<unsigned int> lengthArg( "l", "length", 
    "An integer describing the minimum length a tract must have", 
    false, 0, "Integer length",cmd);

    TCLAP::UnlabeledValueArg<std::string> outputImageArg("output", 
            "Output tract file name (.vtk, .vtp or .gii)", true, 
            "", 
            "Output tract file name",
            cmd); 

    // Now parse the command line and assign variable values as required.
    cmd.parse( argc, argv );

  if ( doTouchArg.isSet() ) {
    std::string ss = doTouchArg.getValue();
    const char *s = ss.c_str();
    char *end;
    unsigned int l = 0;
    while ( (static_cast<unsigned char>(*s) != '\0') &&
            isspace(static_cast<unsigned char>(*s)) ) s++;
    do {
      l = static_cast<unsigned int>( strtol(s, &end, 10));
      if (s == end) break; // we are done converting
      touchlabels.push_back( l );
      s = end;  // if not done, then start again where we left off
      while ( (static_cast<unsigned char>(*s) != '\0') &&
            isspace(static_cast<unsigned char>(*s)) ) s++;
    } while (static_cast<unsigned char>(*s) != '\0');
  }
//std::cout << "Touch labels are " << touchlabels << std::endl;

  if ( dontTouchArg.isSet() ) {
    std::string ss = dontTouchArg.getValue();
    const char *s = ss.c_str();
    char *end;
    unsigned int l = 0;
    while ( (static_cast<unsigned char>(*s) != '\0') &&
            isspace(static_cast<unsigned char>(*s)) ) s++;
    do {
      l = static_cast<unsigned int>( strtol(s, &end, 10));
      if (s == end) break; // we are done converting
      donttouchlabels.push_back( l );
      s = end;  // if not done, then start again where we left off
      while ( (static_cast<unsigned char>(*s) != '\0') &&
            isspace(static_cast<unsigned char>(*s)) ) s++;
    } while (static_cast<unsigned char>(*s) != '\0');
  }
//std::cout << "Don't touch labels are " << donttouchlabels << std::endl;

  if ( lengthArg.isSet() ) {
    minimumTractSize = lengthArg.getValue();
  }
  std::cout << "minimumTractSize is " << minimumTractSize << std::endl;

    inputTractFileName = inputImageArg.getValue();
    std::cout << "Input file name is " << inputImageArg.getValue();
std::cout << endl;

    roiImageFileName = roiImageArg.getValue();
    std::cout << "ROI file name is " << roiImageArg.getValue();
std::cout << endl;

    outputTractFileName = outputImageArg.getValue();
    std::cout << "Output file name is " << outputImageArg.getValue();
std::cout << endl;

  } catch (TCLAP::ArgException &e) 
  {
    std::cerr << "error: " << e.error() << " for arg " << 
                  e.argId() << std::endl;
    std::cerr << "This program computes parametric images from tensors." 
        << std::endl;
    exit(1);
  }

  //define the input/output types
  typedef itk::OrientedImage< short int, ImageDimension > ROIImageType;

  //define reader and writer
  typedef itk::ImageFileReader< ROIImageType > ROIImageReaderType;

  //setup the ROI image reader
  ROIImageReaderType::Pointer roireaderPtr = ROIImageReaderType::New();
  roireaderPtr->SetFileName(roiImageFileName.c_str());
  roireaderPtr->Update();

  //Load in the vtk tracts
  vtkPolyData *polyDataTracts;
  if ( (polyDataTracts=crlVtkMeshIO::ReadMesh(inputTractFileName))==NULL ) exit(1);
  vtkCellArray* loadedtracts = polyDataTracts->GetLines();
  std::cout<<"Total Input Tracts: "<<loadedtracts->GetNumberOfCells()<<std::endl;

  //allocate new VTK Polydata to output the filtered tracts
  vtkPolyData* filteredtracts = vtkPolyData::New();
  //vtkPoints* filteredpoints = vtkPoints::New();
  vtkCellArray* filteredtractarray = vtkCellArray::New();

  ROIImageType::IndexType index;  //preallocate for efficiency

  vtkIdType npts;
  vtkIdType* pts;
  vtkPoints* points = polyDataTracts->GetPoints();
  loadedtracts->InitTraversal();

  while( loadedtracts->GetNextCell( npts, pts ) ){
//    std::cout<<std::endl;
//    std::cout << "Number of points here is " << npts << std::endl;

    if (npts < minimumTractSize) {
      // Skip to the next tract.
      // std::cout << "Skipping this tract since " << npts << " < " << minimumTractSize << std::endl;
      continue;
    }

    std::vector<bool> touches(touchlabels.size());
    for (unsigned int i = 0; i < touches.size(); i++) {
      touches[i] = false;
    }
    std::vector<bool> donttouches(donttouchlabels.size());
    for (unsigned int i = 0; i < donttouches.size(); i++) {
      donttouches[i] = false;
    }

    itk::Point<double,3> point;

    for(int currentpointIDindex=0; currentpointIDindex<npts; currentpointIDindex++) {
      double* vertex = points->GetPoint( pts[currentpointIDindex] );
      point[0] = vertex[0]; point[1] = vertex[1]; point[2] = vertex[2];
      roireaderPtr->GetOutput()->TransformPhysicalPointToIndex(point,index);

      //index[0]=static_cast<long int>(vertex[0]);
      //index[1]=static_cast<long int>(vertex[1]);
      //index[2]=static_cast<long int>(vertex[2]);
      ROIImageType::PixelType& roiimagepix = roireaderPtr->GetOutput()->GetPixel( index );

// std::cout << "vertex " << point << " is at index " << index <<
//  " and has label " << roiimagepix << std::endl;

      for (unsigned int i = 0; i < touchlabels.size(); i++) {
        if (roiimagepix == touchlabels[i]) {
          touches[i] = true;
// std::cout << "Setting touches " << touchlabels[i] << " to be true." << std::endl;
          break;
        }
      }
      for (unsigned int i = 0; i < donttouchlabels.size(); i++) {
        if (roiimagepix == donttouchlabels[i]) {
          donttouches[i] = true;
// std::cout << "Setting dont touches " << donttouchlabels[i] << " to be true." << std::endl;
          break;
        }
      }
    } // End loop over this tract.

    // Now decide if we keep the tract or not.
    bool keep = true;
    for (unsigned int i = 0; i < touches.size(); i++) {
      if (touches[i] != true) {
        keep = false;
// std::cout << "Tract does not touch label " << touchlabels[i] << std::endl;
        break;
      }
    }
    for (unsigned int i = 0; i < donttouches.size(); i++) {
      if (donttouches[i] == true) {
        keep = false;
// std::cout << "Tract touches label " << donttouchlabels[i] << " when it shouldnt" << std::endl;
        break;
      }
    }
    if (!keep) {
//      std::cout << "Tract does not meet keeping criteria." << std::endl;
      continue;
    }

    // Now insert the tract into the set of tracts to save.
    // I think this will save the entire tract:
// std::cout << "Keeping the entire tract of " << npts << " points." << std::endl;
    filteredtractarray->InsertNextCell( npts, pts );

  }

  //finish up the vtk polydata
  filteredtracts->SetPoints( polyDataTracts->GetPoints() );
  filteredtracts->SetLines( filteredtractarray );

  //clean up the poly data to remove redundant points
  vtkCleanPolyData* cleaner = vtkCleanPolyData::New();
  cleaner->SetInput( filteredtracts );
  cleaner->SetAbsoluteTolerance( 0.0 );
  cleaner->Update();

  std::cout<<"Total Output Tracts: "<<filteredtracts->GetNumberOfCells()<<std::endl;
  //output the vtk tract container
  //vtkZLibDataCompressor* compressor = vtkZLibDataCompressor::New();
  //vtkXMLPolyDataWriter* tractswriter = vtkXMLPolyDataWriter::New();
  //tractswriter->SetCompressor( compressor );
  //tractswriter->SetDataModeToBinary();

  //vtkPolyDataWriter *tractswriter = vtkPolyDataWriter::New();
 // tractswriter->SetInput( cleaner->GetOutput() );
  //tractswriter->SetFileName( outputTractFileName.c_str() );
  //tractswriter->Write();
  crlVtkMeshIO::WriteMesh( cleaner->GetOutput(), outputTractFileName);

  //cleanup vtk stuff
  //tractsreader->Delete();
  filteredtracts->Delete();
  cleaner->Delete();
  //filteredpoints->Delete();
  filteredtractarray->Delete();
  //tractswriter->Delete();
  //compressor->Delete();

  return EXIT_SUCCESS;
}
