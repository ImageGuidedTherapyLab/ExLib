/*
 * Copyright (c) 2008-2011 Children's Hospital Boston.
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
 * For more information contact: simon.warfield@childrens.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#include "configuration.h"
#include "tclap/CmdLine.h"
#include <iostream>
#include <fstream>

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkVectorImage.h"

#include "crlDWIRaw.h"
#include "crlFileName.h"
#include "crlVectorImageUtils.h"
#include "crlConvertAndWriteScalarImage.h"
#include "crlDocumentedToolOutput.h"
#include "crlDWICommon.h"

using namespace std;
using namespace crl::DWI;

/**********************************************************************************************//**
 * \fn	int main(int argc, char **argv)
 *
 * \brief	Main entry-point for this application. 
 *
 * \author	Benoit Scherrer
 * \date	December 2009
 *
 * \param	argc	Number of command-line arguments. 
 * \param	argv	Array of command-line argument strings. 
 *
 * \return	Exit-code for the process - 0 for success, else an error code. 
*************************************************************************************************/
int main(int argc, char **argv)
{
	string fileInput, fileOutput, fileB0, namePatient;
	bool tensorMode;

	/*-------------------------------------------------
	First parse the command line (with TCLAP)
	-------------------------------------------------*/
	try {
		TCLAP::CmdLine cmd("", ' ', CRKIT_VERSION_STRING );
		crl::DocumentedToolOutput crlToolOutput(
			"CRL Diffusion Weighted Imaging ToolBox.\n"
			"Convert NHDR RAW data to the MedINRIA RAW data.",
			"Benoit Scherrer", 
			"2010" );
		crlToolOutput.SetKeepOriginalDataType(false);
		cmd.setOutput( &crlToolOutput );

		TCLAP::SwitchArg	argTMode	("t","tensor",		"Tensor mode (if not, ODF mode)", cmd, false);
		TCLAP::ValueArg<string> argPatient     ("p","patientname",	"Name of the patient",	false, "Unknown", "name", cmd);
		TCLAP::ValueArg<string> argB0     ("0","B0image",	"If set, use only one file for the B=0 image (BUG MedINRIA 1.9 when more than one B=0 image)",	false, "", "File name", cmd);

		TCLAP::ValueArg<string> argOutput     ("o","output",	"Output DTS file",	true, "", "Input DTS file", cmd);
		TCLAP::ValueArg<string> argInput     ("i","input",	"Input NHDR file",	true, "", "Input NHDR file", cmd);

		// Parse !
		cmd.parse(argc,argv);

		// Get parameters
		fileInput = argInput.getValue();
		fileOutput = argOutput.getValue();
		tensorMode = argTMode.getValue();
		namePatient = argPatient.getValue();
		fileB0 = argB0.getValue();

	}
	catch ( TCLAP::ArgException& e)
	{
		std::cerr << "Error: "<<e.error()<< " for argument " <<e.argId() << std::endl;
		exit(1);
	}

	/*-------------------------------------
	Beginning of the program!
	-------------------------------------*/
	try {
		crl::FileName oFileName(fileOutput);
		crl::FileName iFileName(fileInput);

		if ( iFileName.getExtension()!="nhdr" ) 
			throw itk::ExceptionObject(__FILE__, __LINE__, "The input file must be an NHDR file");

		//-------------------------------------
		//  Read the diffusion data
		//-------------------------------------
		crl::DWI::DWIRaw<unsigned short> rawData;
		rawData.SetVerbose(false);
		cout<<"- Read <"<<fileInput.c_str()<<">..."<<endl;
		rawData.LoadData(fileInput);

		//-------------------------------------
		//  Open output file for writing
		//-------------------------------------
		std::ofstream ofile;
		ofile.open ( fileOutput.c_str(), ofstream::out);
		if  ( !ofile.is_open() ) 
			throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot open the output file for writing", "");


		//-------------------------------------
		// Begin to write
		//-------------------------------------
		ofile<<"<StudyName>"<<endl;
		ofile<<oFileName.getFileName()<<endl;

		ofile<<"<PatientName>"<<endl;
		ofile<<namePatient.c_str()<<endl;


		//-------------------------------------
		// Number of DWI images
		//-------------------------------------
		ofile<<"<NumberOfDTIs>"<<endl;
		if (fileB0=="")
			ofile<< rawData.GetData()->GetNumberOfComponentsPerPixel() <<endl;
		else
			ofile<< rawData.GetNumberOfGradientImages()+1<<endl;


		//-------------------------------------
		// DWI images to include
		//-------------------------------------
		ofile<<"<DTIs>"<<endl;
		std::vector<string> DTItoInclude;

		int idImage = 0;
		char szBuffer[4096];

		crl::FileName outputFileBaseName(fileOutput);
		outputFileBaseName.setExtension("nrrd");	//ensure nrrd extension

		// add the B0 if specified
		if ( fileB0!="" ) 
		{
			// Open the file
			typedef itk::Image<double, 3> B0ImageType;
			itk::ImageFileReader<B0ImageType>::Pointer readerB0 = itk::ImageFileReader<B0ImageType>::New();
			readerB0->SetFileName(fileB0);
			readerB0->Update();

			// Get the filename
			sprintf(szBuffer,"_%04d_B0", idImage);
			crl::FileName fn( outputFileBaseName.getCompleteFilePath_WithSufix(string(szBuffer)));
			DTItoInclude.push_back(fn.getFileNameWithExt());
			idImage++;

			// Write!
			crl::WriteScalarImage<double, 3>( fn.getCompleteFilePath(),  readerB0->GetOutput(), readerB0->GetImageIO()->GetComponentType()  );
		}

		int id=0;
		crl::DWI::GradientDirectionContainerType::Pointer gradientDirection = rawData.GetDiffusionVectors();
		for(crl::DWI::GradientDirectionContainerType::Iterator it = gradientDirection->Begin() ; it != gradientDirection->End(); it++ )
		{
			// Take only zero-norm gradients if fileB0!=""
			if(fileB0=="" || it.Value().one_norm() > 0.0)
			{
				typedef itk::Image<unsigned short, 3> OneImageType;
				OneImageType::Pointer oneImage = crl::ExtractComponentOfVectorImage<unsigned short, 3>( rawData.GetData(),  id );

				// Get the filename
				if ( it.Value().one_norm()==0 ) sprintf(szBuffer,"_%04d_B0", idImage);
				else sprintf(szBuffer,"_%04d", idImage);

				crl::FileName fn( outputFileBaseName.getCompleteFilePath_WithSufix(string(szBuffer)));
				DTItoInclude.push_back(fn.getFileNameWithExt());
				idImage++;

				crl::WriteScalarImage<unsigned short, 3>( fn.getCompleteFilePath(),  oneImage, itk::ImageIOBase::USHORT  );

			}

			id++;
		}	


		for (unsigned int i=0 ; i<DTItoInclude.size() ; i++ )
			ofile<<DTItoInclude[i]<<endl;

		//-------------------------------------
		// Gradient list
		//-------------------------------------
		ofile<<"<GradientList>"<<endl;

		//-------------------------------------
		// Loop on all the gradients and write data
		//-------------------------------------
		GradientDirectionContainerType* gradContainer = rawData.GetDiffusionVectors();


		bool already_extracted=false;	//for single B0 mode
		for(GradientDirectionContainerType::Iterator it = gradContainer->Begin(); it != gradContainer->End(); it++)
		{

			vnl_vector_fixed<double, 3> v = it.Value();	
			double nn = v.two_norm();
			sprintf(szBuffer, "%f %f %f", v[0], v[1], v[2]); 

			if (nn==0) 
			{	
				if (fileB0=="") ofile << szBuffer << endl;
				else if (!already_extracted)
				{
					ofile << szBuffer << endl;
					already_extracted=true;
				}

			}
			else
			{	
				ofile << szBuffer << endl;
			}
		}


		//-------------------------------------
		// Write other information
		//-------------------------------------
		ofile << "<ModelType>" << endl;
		if ( tensorMode )
			ofile << "Tensor" << endl;
		else
			ofile << "fODF" << endl;


		ofile << "<Tensors>" << endl;
		ofile << "Not_Computed" << endl;
		ofile << "<AllFibers>" << endl;
		ofile << "Not_Computed" << endl;
		ofile << "<NumberOfFibers>" << endl;
		ofile << "0" << endl;
		ofile << "<BST>" << endl;
		ofile << "200" << endl;
		ofile << "<Smoothness>" << endl;
		ofile << "0.2" << endl;
		ofile << "<MinLength>" << endl;
		ofile << "10" << endl;
		ofile << "<MaxLength>" << endl;
		ofile << "200" << endl;
		ofile << "<FAThreshold>" << endl;
		ofile << "0.3" << endl;
		ofile << "<TimeStep>" << endl;
		ofile << "2" << endl;
		ofile << "<UseTriLinearInterpolation>" << endl;
		ofile << "1" << endl;
		ofile << "<Sampling>" << endl;
		ofile << "1" << endl;
		ofile << "<AffineTransform>" << endl;
		ofile << "1 0 0 0 1 0 0 0 1 0 0 0 " << endl;
		ofile << "<FlipAxes>" << endl;
		ofile << "0 0 0" << endl;
		ofile.close();


	}
	/*-------------------------------------
	If catched an exception, show a message
	-------------------------------------*/
	catch (itk::ExceptionObject& e)
	{
		std::cout<< "ITK EXCEPTION ERROR CAUGHT"<<std::endl<< e.GetDescription() << std::endl << "Cannot continue." << std::endl ;
		return 1;
	}
	catch ( ... )
	{
		std::cout<<"UNKNOWN EXCEPTION ERROR." << std::endl << "Cannot continue."<< std::endl;
		return 1;
	}
	cout << endl;
	return 0;
}
