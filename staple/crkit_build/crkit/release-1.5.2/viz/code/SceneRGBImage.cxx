#include <SceneRGBImage.h>

#include <itkOrientImageFilter.h>
#include <itkImageRegionConstIterator.h>

#include <vtkImageData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkImageFlip.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageCast.h>

#include <SceneTensors.h>
#include <SceneImageView.h>

SceneRGBImage::SceneRGBImage()
{
  tensors = 0;
  siv = 0;
  cachedImage = 0;
}


SceneRGBImage::SceneRGBImage(const SceneRGBImage &o)
{
  name = o.GetName();
  cachedImage = 0;
  siv = 0;
  tensors = 0;
}


SceneRGBImage::SceneRGBImage(const QString newName)
{
  name = newName;
  cachedImage = 0;
  siv = 0;
  tensors = 0;
}


SceneRGBImage::SceneRGBImage(QString newName, QString loadFileName)
{
  name = newName;
  if (!LoadImage(loadFileName))
    {
      fileName = "";
    }
  cachedImage = 0;
  siv = 0;
  tensors = 0;
}


SceneRGBImage::~SceneRGBImage()
{
  if (cachedImage)
    cachedImage->Delete();
  if (tensors)
    delete tensors;
  if (siv)
    delete siv;
}


void SceneRGBImage::SetName(const char *n1)
{
  name = QString(n1);
}


void SceneRGBImage::SetName(QString n1)
{
  name = n1;
}


QString SceneRGBImage::GetName() const
{
  return name;
}


bool SceneRGBImage::LoadImage(QString loadFileName)
{
  fileName = loadFileName;

  RGBReaderType::Pointer reader = RGBReaderType::New();

  reader->SetFileName(fileName.toStdString().c_str());
  try
    {
      reader->Update();
    }
  catch (itk::ExceptionObject &e)
    {
      return false;
    }

  if (reader->GetOutput()->GetNumberOfComponentsPerPixel() < 4)
    {
      itk::OrientImageFilter<RGBImageType,RGBImageType>::Pointer orienter = itk::OrientImageFilter<RGBImageType,RGBImageType>::New();

      orienter->UseImageDirectionOn();
      orienter->SetDesiredCoordinateOrientationToAxial();
      orienter->SetInput(reader->GetOutput());
      image = orienter->GetOutput();
      image->Update();

      image->DisconnectPipeline();
      imageDirection = image->GetDirection();
    }
  else
    {
      tensors = new SceneTensors();
      tensors->LoadTensors(loadFileName);
      imageDirection = tensors->GetImageDirection();
      image = NULL;
    }

  return true;
}


vtkImageData *SceneRGBImage::GetVTKImageData()
{
  if (cachedImage)
    return cachedImage;

  if (tensors)
    {
      cachedImage = vtkImageData::New();
      cachedImage->DeepCopy(tensors->GetRGBImage());
    }
  else
    {
      cachedImage = vtkImageData::New();
      cachedImage->SetScalarTypeToUnsignedChar();
      itk::ImageRegion<3> region = image->GetLargestPossibleRegion();

      cachedImage->SetExtent(region.GetIndex(0),region.GetIndex(0)+region.GetSize(0)-1,region.GetIndex(1),region.GetIndex(1)+region.GetSize(1)-1,region.GetIndex(2),region.GetIndex(2)+region.GetSize(2)-1);
      //cachedImage->SetOrigin(0,0,0);
      cachedImage->SetOrigin(image->GetOrigin()[0],image->GetOrigin()[1],image->GetOrigin()[2]);
      cachedImage->SetSpacing(image->GetSpacing()[0],image->GetSpacing()[1],image->GetSpacing()[2]);
      cachedImage->SetNumberOfScalarComponents(4);

      vtkUnsignedCharArray *scalarArray = vtkUnsignedCharArray::New();
      scalarArray->SetNumberOfComponents(4);

      itk::ImageRegionConstIterator<RGBImageType> it(image,image->GetLargestPossibleRegion());

      double maxMagnitude[3] = {0.0,0.0,0.0};
      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
	  if (it.Get()[0] > maxMagnitude[0]) maxMagnitude[0] = it.Get()[0];
	  if (it.Get()[1] > maxMagnitude[1]) maxMagnitude[1] = it.Get()[1];
	  if (it.Get()[2] > maxMagnitude[2]) maxMagnitude[2] = it.Get()[2];
	}

      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
	  float RGB[4];
	  RGB[0] = fabs(floor(it.Get()[0] / maxMagnitude[0] * 255));
	  RGB[1] = fabs(floor(it.Get()[1] / maxMagnitude[1] * 255));
	  RGB[2] = fabs(floor(it.Get()[2] / maxMagnitude[2] * 255));
	  RGB[3] = 255.0;
	  scalarArray->InsertNextTuple(RGB);
	}

      cachedImage->GetPointData()->SetScalars(scalarArray);
      scalarArray->Delete();
    }
  return cachedImage;
}


SceneImageView *SceneRGBImage::GetSceneImageView()
{
  if (!siv)
    siv = new SceneImageView(this);

  return siv;
}


RGBImageType::Pointer SceneRGBImage::GetImage()
{
  return image;
}


RGBImageType::DirectionType SceneRGBImage::GetImageDirection()
{
  return imageDirection;
}
