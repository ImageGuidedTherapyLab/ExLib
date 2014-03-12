
/* (c) Simon Warfield simon.warfield@childrens.harvard.edu 2006 */

/*
 * This class generates a VTK renderable representation of a SceneImage
 * suitable for injecting into the Render2D viewing class.
 */

#include <ImageTypeDefinitions.h>
#include <ImageOperations.h>
#include <SceneImageView.h>
#include <SceneImage.h>
#include <SceneRGBImage.h>

#include <vtkImageData.h>
#include <vtkLookupTable.h>
#include <vtkImageMapToWindowLevelColors.h>
#include <vtkImageActor.h>
#include <vtkMatrix4x4.h>
#include <vtkImageResample.h>
#include <vtkImageChangeInformation.h>
#include <itkRegionOfInterestImageFilter.h>


SceneImageView::SceneImageView(SceneImage *ssi)
{
  si = ssi;
  name = ssi->GetName();
  cached = false;
  vtkimage = 0;
  ici = 0;
  dataType = ImageTypeDefinitions::greyscale;
  lut = 0;
}

SceneImageView::SceneImageView(SceneRGBImage *sri)
{
  si = sri;
  name = sri->GetName();
  cached = false;
  vtkimage = 0;
  ici = 0;
  dataType = ImageTypeDefinitions::rgb;
  lut = 0;
}

SceneImageView::~SceneImageView()
{
  if (cached && vtkimage)
    vtkimage->Delete();
  if (ici)
    ici->Delete();
  if (lut)
    lut->Delete();
}

void SceneImageView::GetInitialWindowLevel(double &window, double &level)
{
  if ( (!cached) || (!GetVTKImageData()) ) return;
  level = mean;
  window = 6.0*stddev;
}

double SceneImageView::GetInitialWindow()
{
  if ( (!cached) || (!GetVTKImageData()) ) return 0.0;
  return 6.0*stddev;
}

double SceneImageView::GetInitialLevel()
{
  if ( (!cached) || (!GetVTKImageData()) ) return 0.0;
  return mean;
}

vtkImageData *SceneImageView::GetVTKImageData()
{
  if (dataType == ImageTypeDefinitions::rgb)
    {
      if (!cached)
	{
	  vtkimage = vtkImageData::New();
	  //vtkimage->DeepCopy(dynamic_cast<SceneRGBImage*>(si)->GetVTKImageData());
	  ici = vtkImageChangeInformation::New();
	  ici->SetInput(dynamic_cast<SceneRGBImage*>(si)->GetVTKImageData());
	  ici->CenterImageOff();
	  ici->Update();
	  vtkimage->DeepCopy(ici->GetOutput());
	  cached = true;
	  
	  ImageTypeDefinitions::ImageType::IndexType idx;
	  ImageTypeDefinitions::ImageType::PointType point;
	  idx[0]=0; idx[1]=0; idx[2]=0;
	  if (dynamic_cast<SceneRGBImage*>(si)->GetImage())
	    {
	      dynamic_cast<SceneRGBImage*>(si)->GetImage()->TransformIndexToPhysicalPoint(idx,point);
	      translation[0] = point[0];
	      translation[1] = point[1];
	      translation[2] = point[2];
	    }
	  else
	    {
	      translation[0] = vtkimage->GetOrigin()[0];
	      translation[1] = vtkimage->GetOrigin()[1];
	      translation[2] = vtkimage->GetOrigin()[2];
	    }
	  
	  ImageTypeDefinitions::ImageType::DirectionType direction = dynamic_cast<SceneRGBImage*>(si)->GetImageDirection();
	  rotation = vtkMatrix4x4::New();
	  rotation->SetElement(0,0,direction[0][0]);
	  rotation->SetElement(0,1,direction[0][1]);
	  rotation->SetElement(0,2,direction[0][2]);
	  rotation->SetElement(1,0,direction[1][0]);
	  rotation->SetElement(1,1,direction[1][1]);
	  rotation->SetElement(1,2,direction[1][2]);
	  rotation->SetElement(2,0,direction[2][0]);
	  rotation->SetElement(2,1,direction[2][1]);
	  rotation->SetElement(2,2,direction[2][2]);
	  
	  rotation->Invert();
	}
      return vtkimage;
    }
  
  if (!si) {
    std::cerr << "SceneImage si is NULL" << std::endl;
    return NULL;
  }
  if (!dynamic_cast<SceneImage*>(si)->inputImage) {
    std::cerr << "inputImage SceneImage si->inputImage is NULL" << std::endl;
    return NULL;
  }
  if (!cached) {
    ImageOperations::calculateImageStatistics( dynamic_cast<SceneImage*>(si)->inputImage,
					       min, max, mean, stddev);
    
    // Make a new ITK2VTK connector for this input
    ImageTypeDefinitions::ITK2VTKConnectorFilterType::Pointer
      ITK2VTKConnector =
      ImageTypeDefinitions::ITK2VTKConnectorFilterType::New();
    ITK2VTKConnector->SetInput( dynamic_cast<SceneImage*>(si)->inputImage );
    
    // Now cache the data from the pipeline
    vtkimage = vtkImageData::New();
    //vtkimage->DeepCopy(ITK2VTKConnector->GetOutput());
    
    // ITK2VTKConnector - smart pointers can be allowed to pass out of scope
    // where they delete themselves.
    
    // The correct data is now cached. Unless something happens to invalidate
    // the data, we will keep using it.
    cached = true;
    
    
    // ici is only here for dummy purposes. For some reason the VTK rendering pipeline hangs when using the output
    // of ITK2VTKConnector directly.
    ici = vtkImageChangeInformation::New();
    ici->SetInput(ITK2VTKConnector->GetOutput());
    ici->CenterImageOff();
    ici->Update();
    vtkimage->DeepCopy(ici->GetOutput());
    
    //ici->GetOutput()->GetOrigin(o);
    ImageTypeDefinitions::ImageType::IndexType idx;
    ImageTypeDefinitions::ImageType::PointType point;
    idx[0]=0; idx[1]=0; idx[2]=0;
    dynamic_cast<SceneImage*>(si)->inputImage->TransformIndexToPhysicalPoint(idx,point);
    translation[0] = point[0];
    translation[1] = point[1];
    translation[2] = point[2];
    
    ImageTypeDefinitions::ImageType::DirectionType direction = dynamic_cast<SceneImage*>(si)->GetImageDirection();
    rotation = vtkMatrix4x4::New();
    rotation->SetElement(0,0,direction[0][0]);
    rotation->SetElement(0,1,direction[0][1]);
    rotation->SetElement(0,2,direction[0][2]);
    rotation->SetElement(1,0,direction[1][0]);
    rotation->SetElement(1,1,direction[1][1]);
    rotation->SetElement(1,2,direction[1][2]);
    rotation->SetElement(2,0,direction[2][0]);
    rotation->SetElement(2,1,direction[2][1]);
    rotation->SetElement(2,2,direction[2][2]);
    
    rotation->Invert();
  }
  
  return vtkimage;
}

void SceneImageView::SetDataType(ImageTypeDefinitions::DataType type)
{
  this->dataType = type;
  if (lut)
    {
      lut->Delete();
      lut = NULL;
    }
}

ImageTypeDefinitions::DataType SceneImageView::GetDataType()
{
  return this->dataType;
}

ImageTypeDefinitions::ImageType::Pointer SceneImageView::ExtractROI(int *bounds)
{
  if (dataType == ImageTypeDefinitions::rgb)
    return NULL;
  
  typedef itk::RegionOfInterestImageFilter<ImageTypeDefinitions::ImageType,ImageTypeDefinitions::ImageType> FilterType;
  
  FilterType::Pointer filter = FilterType::New();
  
  ImageTypeDefinitions::ImageType::RegionType regionROI,regionImage;
  regionImage = dynamic_cast<SceneImage*>(si)->inputImage->GetLargestPossibleRegion();
  
  // let's check if the requested ROI is possible
  if (bounds[0] < regionImage.GetIndex()[0])
    bounds[0] = regionImage.GetIndex()[0];
  if (bounds[2] < regionImage.GetIndex()[1])
    bounds[2] = regionImage.GetIndex()[1];
  if (bounds[4] < regionImage.GetIndex()[2])
    bounds[4] = regionImage.GetIndex()[2];
  if (bounds[1] > regionImage.GetIndex()[0]+(int)regionImage.GetSize()[0])
    bounds[1] = regionImage.GetIndex()[0]+regionImage.GetSize()[0];
  if (bounds[3] > regionImage.GetIndex()[1]+(int)regionImage.GetSize()[1])
    bounds[3] = regionImage.GetIndex()[1]+regionImage.GetSize()[1];
  if (bounds[5] > regionImage.GetIndex()[2]+(int)regionImage.GetSize()[2])
    bounds[5] = regionImage.GetIndex()[2]+regionImage.GetSize()[2];
  
  
  regionROI.SetIndex(0,bounds[0]);
  regionROI.SetIndex(1,bounds[2]);
  regionROI.SetIndex(2,bounds[4]);
  regionROI.SetSize(0,bounds[1]-bounds[0]);
  regionROI.SetSize(1,bounds[3]-bounds[2]);
  regionROI.SetSize(2,bounds[5]-bounds[4]);
  
  
  filter->SetRegionOfInterest(regionROI);
  filter->SetInput(dynamic_cast<SceneImage*>(si)->inputImage);
  try
    {
      filter->Update();
    }
  catch (itk::ExceptionObject &e)
    {
      std::cout << "ROI could not be extracted." << std::endl;
      std::cout << e << std::endl;
      return NULL;
    }
  
  return filter->GetOutput();
}

double *SceneImageView::GetTranslation()
{
  return translation;
}

vtkMatrix4x4 *SceneImageView::GetRotation()
{
  return rotation;
}

SceneData *SceneImageView::GetData()
{
  return si;
}

vtkLookupTable* SceneImageView::GetLookupTable(double alpha)
{
  if (!lut)
    {
      lutAlpha = alpha;
      if (dataType == ImageTypeDefinitions::greyscale)
	{
	  return(NULL);
	}
      
      if (dataType == ImageTypeDefinitions::probmap)
	{
	  // create a color lookup table for probability maps
	  lut = vtkLookupTable::New();
	  //vtklut->SetTableRange(siv->GetVTKImageData()->GetScalarRange()[0],siv->GetVTKImageData()->GetScalarRange()[1]);
	  
	  // we can be certain that the data is of type SceneImage* for probability maps
	  SceneImage *image = dynamic_cast<SceneImage*>(si);
	  
	  lut->SetTableRange(image->GetLowerValue(),image->GetUpperValue());
	  lut->SetHueRange(image->GetColor(SceneImage::lower).hueF(),image->GetColor(SceneImage::upper).hueF());
	  lut->SetSaturationRange(image->GetColor(SceneImage::lower).saturationF(),image->GetColor(SceneImage::upper).saturationF());
	  lut->SetValueRange(image->GetColor(SceneImage::lower).valueF(),image->GetColor(SceneImage::upper).valueF());
	  if (image->GetLogLookupTable())
	    lut->SetScaleToLog10();
	  lut->Build();
	  lut->SetTableValue(lut->GetIndex(0),0,0,0,alpha);
	  return(lut);
	}
      
      if (dataType == ImageTypeDefinitions::fmri)
	{
	  // create a color lookup table for fmri
	  lut = vtkLookupTable::New();
	  //vtklut->SetTableRange(0.001,imageData->GetScalarRange()[1]);
	  //vtklut->SetTableRange(0.001,20);
	  
	  // the data has to be a SceneImage* for fMRI images
	  SceneImage *image = dynamic_cast<SceneImage*>(si);
	  lut->SetTableRange(image->GetLowerValue(),image->GetUpperValue());
	  lut->SetHueRange(0.0,0.18);
	  lut->SetSaturationRange(1,1);
	  lut->SetValueRange(1,1);
	  lut->Build();
	  lut->SetTableValue(lut->GetIndex(0),0,0,0,alpha);
	  return(lut);
	}
      
      if (dataType == ImageTypeDefinitions::segmentation)
	{
	  // create a lookup table that matches the standard CRL colors for baby segmentations
	  lut = vtkLookupTable::New();
	  lut->SetNumberOfColors(256);
	  lut->SetTableRange(0,255);
	  lut->Build();
	  lut->SetTableValue(0,0.0,0.0,0.0,alpha);
	  for (unsigned int n=0; n < 255; n+=17)
	    {
	      lut->SetTableValue(n+1,0.7,0.0,1.0);
	      lut->SetTableValue(n+2,0.0,1.0,1.0);
	      lut->SetTableValue(n+3,1.0,0.55,0.55);
	      lut->SetTableValue(n+4,0.59,0.59,0.59);
	      lut->SetTableValue(n+5,0.0,0.0,1.0);
	      lut->SetTableValue(n+6,1.0,0.55,0.0);
	      lut->SetTableValue(n+7,1.0,0.0,0.0);
	      lut->SetTableValue(n+8,0.97,0.97,1.0);
	      lut->SetTableValue(n+9,0.64,0.32,0.0);
	      lut->SetTableValue(n+10,1.0,0.0,1.0);
	      lut->SetTableValue(n+11,0.5,0.5,0.5);
	      lut->SetTableValue(n+12,1.0,0.66,0.14);
	      lut->SetTableValue(n+13,0.5,1.0,0.0);
	      lut->SetTableValue(n+14,0.5,1.0,0.83);
	      lut->SetTableValue(n+15,0.39,0.58,0.93);
	      lut->SetTableValue(n+16,0.54,0.17,0.89);
	      lut->SetTableValue(n+17,0.99,0.99,0.99);
	    }
	  return(lut);
	}
      
      return(NULL);
    }
  else
    return lut;
}

void SceneImageView::UpdateLookupTable()
{
  if (!lut)
    return;
  if ((dataType == ImageTypeDefinitions::greyscale) || (dataType == ImageTypeDefinitions::segmentation))
    return;
  
  if (dataType == ImageTypeDefinitions::fmri)
    {
      SceneImage *image = dynamic_cast<SceneImage*>(si);
      lut->SetTableRange(image->GetLowerValue(),image->GetUpperValue());
      lut->ForceBuild();
      lut->SetTableValue(lut->GetIndex(0),0,0,0,lutAlpha);
    }
  if (dataType == ImageTypeDefinitions::probmap)
    {
      SceneImage *image = dynamic_cast<SceneImage*>(si);
      lut->SetTableRange(image->GetLowerValue(),image->GetUpperValue());
      lut->SetHueRange(image->GetColor(SceneImage::lower).hueF(),image->GetColor(SceneImage::upper).hueF());
      lut->SetSaturationRange(image->GetColor(SceneImage::lower).saturationF(),image->GetColor(SceneImage::upper).saturationF());
      lut->SetValueRange(image->GetColor(SceneImage::lower).valueF(),image->GetColor(SceneImage::upper).valueF());
      if (image->GetLogLookupTable())
	lut->SetScaleToLog10();
      lut->ForceBuild();
      lut->SetTableValue(lut->GetIndex(0),0,0,0,lutAlpha);
    }
}
