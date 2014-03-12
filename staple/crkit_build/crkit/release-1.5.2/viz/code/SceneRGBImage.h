#ifndef _SCENERGBIMAGE_H_
#define _SCENERGBIMAGE_H_

#include <SceneData.h>

#include <QString>

#include <itkOrientedImage.h>
#include <itkImageFileReader.h>

class vtkImageData;
class SceneTensors;
class SceneImageView;


typedef itk::VectorImage<float,3> RGBImageType;
typedef itk::ImageFileReader<RGBImageType> RGBReaderType;


class SceneRGBImage : public SceneData
{
public:
  SceneRGBImage();
  SceneRGBImage(const SceneRGBImage &o);
  SceneRGBImage(const QString newName);
  SceneRGBImage(QString newName, QString loadFileName);

  virtual ~SceneRGBImage();

  void SetName(const char *n1);
  void SetName(QString n1);
  QString GetName() const;

  bool LoadImage(QString loadFileName);

  vtkImageData *GetVTKImageData();
  RGBImageType::DirectionType GetImageDirection();

  SceneImageView *GetSceneImageView();
  RGBImageType::Pointer GetImage();

protected:
  QString name;
  QString fileName;

  SceneTensors *tensors;
  RGBImageType::Pointer image;
  RGBImageType::DirectionType imageDirection;
  vtkImageData *cachedImage;
  SceneImageView *siv;
};

#endif
