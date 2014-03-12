#include <SceneTensors.h>
#include <SceneImageView.h>

// VTK classes
#include <vtkStructuredPoints.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkSuperquadricSource.h>
#include <vtkSphereSource.h>
#include <vtkLineSource.h>
#include <vtkPolyDataNormals.h>
#include <vtkExtractVOI.h>
#include <vtkImageChangeInformation.h>
#include <vtkUnsignedCharArray.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <itkDiffusionTensor3D.h>
#include <itkImageFileReader.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>
#include <itkOrientImageFilter.h>

#include <crlTensorGlyph.h>

typedef float InternalPixelType;

typedef itk::DiffusionTensor3D<InternalPixelType> PixelType;
typedef itk::Image<PixelType,3> TensorImageType;
typedef itk::ImageFileReader<TensorImageType> TensorReaderType;
typedef itk::ImageRegionConstIteratorWithIndex<TensorImageType> TensorIteratorType;

SceneTensors::SceneTensors()
{
  sp = 0;
  glyphs = 0;
  extractVOI = 0;
  ici = 0;
  ellipNormals = 0;
  tpdf = 0;
  orientation = ImageTypeDefinitions::axial;
  scale = 1000.0;
}

// copy constructor
SceneTensors::SceneTensors(const SceneTensors &o) {
  name = o.GetName();
  sp = o.sp;
  glyphs = o.glyphs;
  extractVOI = o.extractVOI;
  ici = o.ici;
  ellipNormals = o.ellipNormals;
  tpdf = o.tpdf;
  orientation = o.orientation;
  scale = o.scale;
}

// Construct with a specified handle (name)
SceneTensors::SceneTensors(const QString newName)
{
  name = newName;
  sp = 0;
  glyphs = 0;
  extractVOI = 0;
  ici = 0;
  ellipNormals = 0;
  tpdf = 0;
  orientation = ImageTypeDefinitions::axial;
  scale = 1000.0;
}

// Construct with a specified handle (name) and load an image
SceneTensors::SceneTensors(QString newName, QString loadFileName)
{
  name = newName;
  sp = 0;
  glyphs = 0;
  if (!LoadTensors(loadFileName)) {
    fileName = "";
  }
  extractVOI = 0;
  ici = 0;
  ellipNormals = 0;
  tpdf = 0;
  orientation = ImageTypeDefinitions::axial;
  scale = 1000.0;
}

// destructor
SceneTensors::~SceneTensors() {
  if (sp)
    sp->Delete();
  if (glyphs)
    glyphs->Delete();
  if (extractVOI)
    extractVOI->Delete();
  if (ellipNormals)
    ellipNormals->Delete();
  if (tpdf)
    tpdf->Delete();
  if (ici)
    ici->Delete();

  sp = 0;
  glyphs = 0;
  extractVOI = 0;
  ici = 0;
  ellipNormals = 0;
}

void SceneTensors::SetName(const char *n1)
{
  name = QString(n1);
}

void SceneTensors::SetName(QString n1)
{
  name = n1;
}

QString SceneTensors::GetName() const
{
  return name;
}

// This loads all tensor images into memory
bool SceneTensors::LoadTensors(QString loadFileName) {
  fileName = loadFileName;

  if (sp)
    sp->Delete();
  sp = vtkImageData::New();

  TensorReaderType::Pointer reader = TensorReaderType::New();

  reader->SetFileName(fileName.toStdString().c_str());
  try
    {
      reader->Update();
    }
  catch (itk::ExceptionObject &e)
    {
      return false;
    }

  itk::OrientImageFilter<TensorImageType,TensorImageType>::Pointer orienter = itk::OrientImageFilter<TensorImageType,TensorImageType>::New();
  orienter->SetInput(reader->GetOutput());
  orienter->UseImageDirectionOn();
  orienter->SetDesiredCoordinateOrientationToAxial();
  orienter->Update();
  TensorImageType::Pointer tensorImage = orienter->GetOutput();
  itk::ImageRegion<3> region = tensorImage->GetLargestPossibleRegion();
  imageDirection = orienter->GetOutput()->GetDirection();
  TensorImageType::IndexType idx;
  TensorImageType::PointType point;
  idx[0]=0; idx[1]=0; idx[2]=0;
  orienter->GetOutput()->TransformIndexToPhysicalPoint(idx,point);
  translation[0] = point[0]; translation[1] = point[1]; translation[2] = point[2];

  extent[0] = region.GetIndex(0);
  extent[1] = region.GetIndex(0)+region.GetSize(0)-1;
  extent[2] = region.GetIndex(1);
  extent[3] = region.GetIndex(1)+region.GetSize(1)-1;
  extent[4] = region.GetIndex(2);
  extent[5] = region.GetIndex(2)+region.GetSize(2)-1;

  currentSlice = 0;
  sp->SetScalarTypeToUnsignedChar();
  sp->SetExtent(extent);
  sp->SetOrigin(tensorImage->GetOrigin()[0],tensorImage->GetOrigin()[1],tensorImage->GetOrigin()[2]);
  sp->SetSpacing(tensorImage->GetSpacing()[0],tensorImage->GetSpacing()[1],tensorImage->GetSpacing()[2]);
  sp->SetNumberOfScalarComponents(4);

  vtkFloatArray* tensorArray = vtkFloatArray::New();
  tensorArray->SetNumberOfComponents(9);
  /* Don't set the number of tuples, and use the insert we do.
     Either one or the other.
     We may want to insert these into the correct position in some way...
     tensorArray->SetNumberOfTuples(
     region.GetSize(0)*
     region.GetSize(1)*
     region.GetSize(2) );
  */
  /*std::cout << "Size of tensor array is " << (region.GetSize(0)*
    region.GetSize(1)*
    region.GetSize(2) ) << std::endl;*/

  vtkUnsignedCharArray *scalarArray = vtkUnsignedCharArray::New();
  scalarArray->SetNumberOfComponents(4);

  itk::DiffusionTensor3D<InternalPixelType> itkTensor;
  float vtkTensor[9];

  TensorIteratorType it(tensorImage, tensorImage->GetLargestPossibleRegion());

  itk::DiffusionTensor3D<InternalPixelType>::EigenValuesArrayType eigenValues;
  itk::DiffusionTensor3D<InternalPixelType>::EigenVectorsMatrixType eigenVectors;

  double maxfa = 0.0;
  double maxmd = 0.0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      itkTensor = it.Get();

      itkTensor.ComputeEigenAnalysis(eigenValues,eigenVectors);
      // Each ROW of the matrix V represents an eigenvector, which is
      // different to matlab. In matlab each column is an eigenvector.
      // The ordering of the eigenvalues generated is from smallest to largest.
      // Index 0 has the smallest eigenvalue, index 2 has the largest.

      float fa = itkTensor.GetFractionalAnisotropy();
      if (fa < 0) fa = 0.0;
      if (fa > 1) fa = 1.0;

      // Compute the mean diffusivity
      double md = (eigenValues[0] + eigenValues[1] + eigenValues[2])/3.0;

      if (fa > maxfa) maxfa = fa;
      if (md > maxmd) maxmd = md;
    }

  TensorImageType::IndexType index;

  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      itkTensor = it.Get();
      index = it.GetIndex();
      /*if ( (index[0] == 51) && (index[1] == 31) && (index[2] == 38) ) {
        std::cout << "tensor at 51,31,38 is " << std::endl;
        std::cout << itkTensor << std::endl;
	}*/

      itkTensor.ComputeEigenAnalysis(eigenValues,eigenVectors);
      // Each ROW of the matrix V represents an eigenvector, which is
      // different to matlab. In matlab each column is an eigenvector.
      // The ordering of the eigenvalues generated is from smallest to largest.
      // Index 0 has the smallest eigenvalue, index 2 has the largest.

      float fa = itkTensor.GetFractionalAnisotropy();
      if (fa < 0) fa = 0.0;
      if (fa > 1) fa = 1.0;

      // Normalize the fa to the max fa of the image, in order to reduce
      // the color intensity in low fa regions.
      fa = fa / maxfa;

      // Compute the mean diffusivity
      //double md = (eigenValues[0] + eigenValues[1] + eigenValues[2])/3.0;

      /*if ( (index[0] == 51) && (index[1] == 31) && (index[2] == 38) ) {
        std::cout << "fa is " << fa << std::endl;
        std::cout << "md is " << md << std::endl;
        std::cout << " e1 is " << eigenValues[2] <<
	" e2 is " << eigenValues[1] <<
	" e3 is " << eigenValues[0] <<
	std::endl;
        std::cout << "eigenvector[2] is " <<
	eigenVectors[2][0] << " " <<
	eigenVectors[2][1] << " " <<
	eigenVectors[2][2] << " " <<
        std::endl;
        std::cout << "eigenvector[1] is " <<
	eigenVectors[1][0] << " " <<
	eigenVectors[1][1] << " " <<
	eigenVectors[1][2] << " " <<
        std::endl;
        std::cout << "eigenvector[0] is " <<
	eigenVectors[0][0] << " " <<
	eigenVectors[0][1] << " " <<
	eigenVectors[0][2] << " " <<
        std::endl;
	}*/

      // Let's ensure the eigenvectors are normalized to have unit magnitude
      // so that dot products can be easily computed.
      for (unsigned int n = 0; n < 3; n++)
	{
          double norm = eigenVectors[n][0] * eigenVectors[n][0] +
	    eigenVectors[n][1] * eigenVectors[n][1] +
	    eigenVectors[n][2] * eigenVectors[n][2] ;

          norm = sqrt(norm);
	  eigenVectors[n][0] /= norm;
	  eigenVectors[n][1] /= norm;
	  eigenVectors[n][2] /= norm;
	}

      float RGB[4];
      // Set the red channel to the projection of the primary eigenvector
      // onto the x axis, y axis and z axis, with appropriate scaling.
      //   Whether the eigenvector is pointing to +x or -x, we want the
      // projection to result in a positive color.
      RGB[0] = fabs(floor(eigenVectors[2][0] * fa * 255));
      RGB[1] = fabs(floor(eigenVectors[2][1] * fa * 255));
      RGB[2] = fabs(floor(eigenVectors[2][2] * fa * 255));
      RGB[3] = 255.0;
      /*if ( (index[0] == 51) && (index[1] == 31) && (index[2] == 38) ) {
        std::cout << " RGB[0] is " << RGB[0]
	<< " RGB[1] is " << RGB[1]
	<< " RGB[2] is " << RGB[2]
	<< std::endl;
	}*/

      // This highlights CSF regions RGB[3] = 255.0 * md / maxmd;
      scalarArray->InsertNextTuple(RGB);

      // Insert the tensors for glyph rendering, using the VTK mode of
      // ExtractEigenvaluesOff. In this mode, we set the eigenvectors and
      // eigenvalues, rather than the tensors.
      //   The column of the tensors is the eigenvector, and its norm,
      // always positive, is the eigenvalue.
      vtkTensor[0] = eigenVectors[2][0]*eigenValues[2];
      vtkTensor[1] = eigenVectors[2][1]*eigenValues[2];
      vtkTensor[2] = eigenVectors[2][2]*eigenValues[2];
      vtkTensor[3] = eigenVectors[1][0]*eigenValues[1];
      vtkTensor[4] = eigenVectors[1][1]*eigenValues[1];
      vtkTensor[5] = eigenVectors[1][2]*eigenValues[1];
      vtkTensor[6] = eigenVectors[0][0]*eigenValues[0];
      vtkTensor[7] = eigenVectors[0][1]*eigenValues[0];
      vtkTensor[8] = eigenVectors[0][2]*eigenValues[0];
      tensorArray->InsertNextTuple(vtkTensor);

    }

  sp->GetPointData()->SetTensors(tensorArray);
  sp->GetPointData()->SetScalars(scalarArray);
  tensorArray->Delete();
  scalarArray->Delete();

  return true; // success
}

vtkPolyData *SceneTensors::GetVTKTensorGlyphs()
{
  if (!sp)
    return NULL;

  if (!extractVOI)
    {
      extractVOI = vtkExtractVOI::New();
      extractVOI->SetInput(sp);
      switch(orientation)
	{
	case ImageTypeDefinitions::axial:
	  extractVOI->SetVOI(extent[0],extent[1],extent[2],extent[3],
			     currentSlice, currentSlice);
	  break;
	case ImageTypeDefinitions::sagittal:
	  extractVOI->SetVOI(currentSlice,currentSlice,extent[2],extent[3],extent[4],extent[5]);
	  break;
	case ImageTypeDefinitions::coronal:
	  extractVOI->SetVOI(extent[0],extent[1],currentSlice,currentSlice,extent[4],extent[5]);
	};
      extractVOI->Update();
    }

  if (!tpdf)
    {
      vtkLineSource *sqs = vtkLineSource::New();
      glyphs = vtkPolyData::New();

      ellipsoids = crlTensorGlyph::New();
      //  extractVOI has been modified to have tensor eigenvectors and eigenvalues
      // rather than directly tensors. This requires ExtractEigenvaluesOff.
      ellipsoids->ExtractEigenvaluesOff();
      ellipsoids->SetInput(extractVOI->GetOutput());
      ellipsoids->SetSource(sqs->GetOutput());
      ellipsoids->SetScaleFactor(scale);
      ellipsoids->ThreeGlyphsOff();
      ellipsoids->SymmetricOff();
      ellipsoids->ColorGlyphsOn();
      ellipsoids->SetColorModeToScalars();
      ellipsoids->Update();

      ellipNormals = vtkPolyDataNormals::New();
      ellipNormals->SetInput(ellipsoids->GetOutput());
      ellipNormals->Update();

      ellipsoids->Delete();
      sqs->Delete();

      // take care of transforms of the input image that were ignored previously
      vtkMatrix4x4 *rotation = vtkMatrix4x4::New();
      rotation->SetElement(0,0,imageDirection[0][0]);
      rotation->SetElement(0,1,imageDirection[0][1]);
      rotation->SetElement(0,2,imageDirection[0][2]);
      rotation->SetElement(1,0,imageDirection[1][0]);
      rotation->SetElement(1,1,imageDirection[1][1]);
      rotation->SetElement(1,2,imageDirection[1][2]);
      rotation->SetElement(2,0,imageDirection[2][0]);
      rotation->SetElement(2,1,imageDirection[2][1]);
      rotation->SetElement(2,2,imageDirection[2][2]);

      tpdf = vtkTransformPolyDataFilter::New();
      vtkTransform *transform = vtkTransform::New();
      tpdf->SetInput(ellipNormals->GetOutput());
      transform->PostMultiply();
      transform->Translate(-translation[0],-translation[1],-translation[2]);
      transform->Concatenate(rotation);
      transform->Translate(translation);
      tpdf->SetTransform(transform);
      tpdf->Update();
      transform->Delete();
    }

  return tpdf->GetOutput();
}

vtkPolyData *SceneTensors::GetVTKTensorGlyphs(SceneImageView *siv)
{
  vtkPolyData *pd = this->GetVTKTensorGlyphs();

  vtkTransform *transform = vtkTransform::New();
  vtkTransformPolyDataFilter *tpdf2 = vtkTransformPolyDataFilter::New();
  vtkMatrix4x4 *rotation = vtkMatrix4x4::New();
  tpdf2->SetInput(pd);
  transform->Translate(-(siv->GetTranslation()[0]),-(siv->GetTranslation()[1]),-(siv->GetTranslation()[2]));
  transform->Concatenate(siv->GetRotation());
  transform->Translate(siv->GetTranslation());
  tpdf2->SetTransform(transform);

  transform->Delete();
  rotation->Delete();
  return tpdf2->GetOutput();
}


vtkImageData *SceneTensors::GetRGBImage()
{
  return sp;
}

int SceneTensors::GetMinSlice()
{
  switch (orientation)
    {
    case ImageTypeDefinitions::sagittal:
      return extent[0];

    case ImageTypeDefinitions::coronal:
      return extent[2];

    default:
      return extent[4];
    }
}

int SceneTensors::GetMaxSlice()
{
  switch (orientation)
    {
    case ImageTypeDefinitions::sagittal:
      return extent[1];

    case ImageTypeDefinitions::coronal:
      return extent[3];

    default:
      return extent[5];
    }
}

void SceneTensors::SetSlice(int slice, ImageTypeDefinitions::OrientationType newOrientation)
{
  if (!extractVOI)
    return;

  orientation = newOrientation;

  switch (orientation)
    {
    case ImageTypeDefinitions::axial:
      if ((slice < extent[4]) || (slice > extent[5]))
	return;
      extractVOI->SetVOI(extent[0],extent[1],extent[2],extent[3],slice,slice);
      break;
    case ImageTypeDefinitions::sagittal:
      if ((slice < extent[0]) || (slice > extent[1]))
	return;
      extractVOI->SetVOI(slice,slice,extent[2],extent[3],extent[4],extent[5]);
      break;
    case ImageTypeDefinitions::coronal:
      if ((slice < extent[2]) || (slice > extent[3]))
	return;
      extractVOI->SetVOI(extent[0],extent[1],slice,slice,extent[4],extent[5]);
    }

  currentSlice = slice;
  tpdf->Update();
}


ImageTypeDefinitions::ImageType::DirectionType SceneTensors::GetImageDirection()
{
  return imageDirection;
}


double SceneTensors::GetScale()
{
  return scale;
}


void SceneTensors::SetScale(double newScale)
{
  scale = newScale;
  if (tpdf)
    {
      ellipsoids->SetScaleFactor(scale);
      ellipsoids->Update();
      tpdf->Update();
    }
}
