/*
 * crlLogEuclideanScaleAffineTransform.cxx
 *
 *  Created on: Jan 23, 2010
 *      Author: weisen
 */

#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>
#include <itkScaleSkewVersor3DTransform.h>
#include <itkSimilarity3DTransform.h>
#include <itkVersorRigid3DTransform.h>
#include <itkTranslationTransform.h>
#include <vnl/algo/vnl_real_eigensystem.h>
#include <vnl/algo/vnl_complex_eigensystem.h>
#include <vnl/algo/vnl_matrix_inverse.h>


const unsigned int Dimension = 3;
typedef itk::AffineTransform< double, Dimension >  AffineTransformType;

typedef itk::TransformFileReader::TransformListType* TransformListTypePointer;


void logEuclideanScaleAffineTransform( AffineTransformType::Pointer& transp, double scale )
{
    typedef vnl_matrix<double >	RealMatrixType;
    typedef vcl_complex<double> ComplexType;
    typedef vnl_matrix<ComplexType > ComplexMatrixType;

    RealMatrixType A(4,4);		// we're going to form a homogeneous matrix
    ComplexMatrixType cA(4,4);

    for ( unsigned int i = 0; i < 3; ++i )
	for ( unsigned int j = 0; j < 3; ++j )
	    cA(i,j) = A(i,j) = transp->GetMatrix()(i,j);

    for ( unsigned int i = 0; i < 3; ++i )
	cA(i,3) = A(i,3) = transp->GetTranslation()[i];

    cA(3,3) = A(3,3) = 1.0;
    A(3,0) = A(3,1) = A(3,2) = 0.;
    cA(3,0) = cA(3,1) = cA(3,2) = 0.;

#if 0
    // debugging NIW
    std::cerr << "original transform: " << std::endl;
    transp->Print(std::cerr);

    std::cerr << "new transform: " << std::endl;
    A.print(std::cerr);

    std::cerr << "new complex transform: " << std::endl;
    cA.print(std::cerr);
#endif

    vnl_real_eigensystem eig(A);

#if 0
    std::cerr << "eigenvectors:" << std::endl;
    eig.V.print(std::cerr);

    std::cerr << "eigenvalues: " << std::endl;
    eig.D.asMatrix().print(std::cerr);
#endif


    ComplexMatrixType Vinv = vnl_matrix_inverse<ComplexType>( eig.V );
    ComplexMatrixType logm = eig.D.asMatrix();
    for ( unsigned int i = 0; i < 4; ++i )
	logm(i,i) = vcl_exp( scale* vcl_log( logm(i,i) ) );

    logm = eig.V * logm * Vinv;

#if 0
    std::cerr << "scaled matrix:" << std::endl;
    logm.print(std::cerr);
#endif

    // repopulate the transform
    AffineTransformType::TranslationType newTrans;
    AffineTransformType::MatrixType newMatrix;

    for ( unsigned int i = 0; i < 3; ++i )
	newTrans[i] = vcl_real<double>( logm(i,3) );

    for ( unsigned int i = 0; i < 3; ++i )
	for ( unsigned int j = 0; j < 3; ++j )
	    newMatrix(i,j) = vcl_real<double>( logm(i,j) );

    transp->SetTranslation(newTrans);
    transp->SetMatrix(newMatrix);


#if 0
    std::cerr << "New transform:" << std::endl;
    transp->Print( std::cerr );
#endif
}




AffineTransformType::Pointer convertTransformListToAffineTransform( TransformListTypePointer transforms )
{

    typedef itk::ScaleSkewVersor3DTransform< double >
                     ScaleSkewVersor3DTransformType;
    typedef itk::Similarity3DTransform< double >
                     Similarity3DTransformType;
    typedef itk::VersorRigid3DTransform< double >
                     VersorRigid3DTransformType;
    typedef itk::TranslationTransform< double, Dimension >
                     TranslationTransformType;

    AffineTransformType::Pointer outputTransform = AffineTransformType::New();
    outputTransform->SetIdentity();

    std::cout << "Number of transforms = " << transforms->size() << std::endl;
    itk::TransformFileReader::TransformListType::const_iterator it =
                 transforms->begin();
    std::ostringstream errs;
    if (transforms->size() <= 0 || transforms->size() > 1)
	{
	errs << "Read " << transforms->size() << " transforms but want 1." << std::endl;
	throw std::runtime_error( errs.str() );
	}
    else if (!strcmp((*it)->GetNameOfClass(), "AffineTransform"))
	{
	AffineTransformType::Pointer affine_read =
		static_cast<AffineTransformType*> ((*it).GetPointer());
	affine_read->Print(std::cout);

	// Setting the matrix and offset is sufficient to duplicate the transform
	// performance operating on points,
	// but does not set the center appropriately for follow on registration.
	outputTransform->SetTranslation(affine_read->GetTranslation());
	outputTransform->SetCenter(affine_read->GetCenter());
	outputTransform->SetMatrix(affine_read->GetMatrix());
	}
    else if (!strcmp((*it)->GetNameOfClass(), "TranslationTransform"))
	{
	TranslationTransformType::Pointer  translation_read =
          static_cast<TranslationTransformType*>((*it).GetPointer());

	outputTransform->SetOffset(translation_read->GetOffset());

	}
    else if (!strcmp((*it)->GetNameOfClass(), "VersorRigid3DTransform"))
	{
	VersorRigid3DTransformType::Pointer  versorrigid_read =
          static_cast<VersorRigid3DTransformType*>((*it).GetPointer());
	outputTransform->SetCenter(versorrigid_read->GetCenter());
	outputTransform->SetTranslation(versorrigid_read->GetTranslation());
	outputTransform->SetMatrix(versorrigid_read->GetMatrix());
	// The above is sufficient to duplicate the transform performance
	// but does not set the center appropriately for follow on registration.
	}
    else if (!strcmp((*it)->GetNameOfClass(), "Similarity3DTransform"))
	{
	Similarity3DTransformType::Pointer  similarity_read =
          static_cast<Similarity3DTransformType*>((*it).GetPointer());
	outputTransform->SetCenter(similarity_read->GetCenter());
	outputTransform->SetTranslation(similarity_read->GetTranslation());
	outputTransform->SetMatrix(similarity_read->GetMatrix());
	}
    else if (!strcmp((*it)->GetNameOfClass(), "ScaleSkewVersor3DTransform"))
	{
	ScaleSkewVersor3DTransformType::Pointer  scaleskewversor_read =
          static_cast<ScaleSkewVersor3DTransformType*>((*it).GetPointer());
	outputTransform->SetCenter(scaleskewversor_read->GetCenter());
	outputTransform->SetTranslation(scaleskewversor_read->GetTranslation());
	outputTransform->SetMatrix(scaleskewversor_read->GetMatrix());
	}
    else
	{
	errs << "Don't know how to convert a " <<
              (*it)->GetNameOfClass() << " transform." << std::endl;
	throw std::runtime_error(errs.str());
	}

    return outputTransform;
}



int main( int argc, char* argv[] )
{
    if ( argc != 4 )
	{
	std::cerr << "usage: " << argv[0] << " input-transform.trsf scale output-transform.trsf" << std::endl;
	std::cerr << "         scale an affine transform in Log-Euclidean space." << std::endl;
	std::cerr << "         0.5 as a scalefactor is the L.E. half-way point" << std::endl;
	std::cerr << "         or the square root A^(0.5) of the matrix." << std::endl;

	return 1;
	}

    const char* inputFile = argv[1];
    float scale;
    std::istringstream(argv[2]) >> scale;
    const char* outputFile = argv[3];

    std::cout << "input file=" << inputFile << std::endl;
    std::cout << "scale=" << scale << std::endl;
    std::cout << "output file=" << outputFile << std::endl;


    //
    // read transform list
    //
    typedef itk::TransformFileReader TransformReader;
    TransformReader::Pointer trsfreader = TransformReader::New();
    trsfreader->SetFileName( inputFile );

    try
	{
	trsfreader->Update();
	}
    catch ( itk::ExceptionObject & excp )
	{
	std::cerr << "Error while reading the transform file " <<
                    argv[1] << std::endl;
	std::cerr << excp << std::endl;
	std::cerr << "[FAILED]" << std::endl;
	return 1;
	}

    //
    // extract affine transform
    //
    AffineTransformType::Pointer outputTransform = convertTransformListToAffineTransform( trsfreader->GetTransformList() );

    //
    // convert
    //
    logEuclideanScaleAffineTransform(outputTransform, scale);


    //
    // write transform output
    //
    typedef itk::TransformFileWriter TransformWriter;
    TransformWriter::Pointer trsfwriter = TransformWriter::New();
    trsfwriter->SetFileName( outputFile );
    trsfwriter->SetInput ( outputTransform );


    try
	{
	trsfwriter->Update();
	}
    catch( itk::ExceptionObject & excep )
	{
	std::cerr << "Exception occurred writing out the transform to file "
	    << argv[2] << std::endl;
	std::cerr << excep << std::endl;
	}


    return 0;
}
