
/* (c) Simon Warfield simon.warfield@childrens.harvard.edu 2006 */

/*
 * This class generates a VTK renderable representation of a SceneImage
 * suitable for injecting into the Render2D viewing class.
 */

#ifndef _SCENEIMAGEVIEW_INCLUDED
#define _SCENEIMAGEVIEW_INCLUDED 1

// Forward declaration of SceneImage
class SceneImage;
class SceneRGBImage;
class SceneData;

#include "ImageTypeDefinitions.h"

#include <itkOrientImageFilter.h>

// Qt classes
#include <QString>

// VTK classes
class vtkImageData;
class vtkImageChangeInformation;
class vtkObject;
class vtkMatrix4x4;
class vtkLookupTable;

class SceneImageView
{

protected:
  QString name;   // identifier in the scene

  ImageTypeDefinitions::ImageType::PixelType min;
  ImageTypeDefinitions::ImageType::PixelType max;
  ImageTypeDefinitions::ImageType::PixelType mean;
  ImageTypeDefinitions::ImageType::PixelType stddev;
  ImageTypeDefinitions::DataType dataType;

  bool cached;

  vtkImageData *vtkimage;
  vtkImageChangeInformation *ici;

private:
  SceneData *si;
  double translation[3];
  vtkMatrix4x4 *rotation;
  vtkLookupTable *lut;
  double lutAlpha;

public:

  SceneImageView(SceneImage *ssi);
  SceneImageView(SceneRGBImage *sri);

  // destructor
  virtual ~SceneImageView();

  void SetName(const char *n1) { name = QString(n1); }
  void SetName(QString n1) { name = n1; }
  QString GetName() const { return name; }

  void GetInitialWindowLevel(double &window, double &level);
  double GetInitialWindow();
  double GetInitialLevel();

  vtkImageData *GetVTKImageData();

  void SetDataType(ImageTypeDefinitions::DataType type);
  ImageTypeDefinitions::DataType GetDataType();

  ImageTypeDefinitions::ImageType::Pointer ExtractROI(int *bounds);

  double *GetTranslation();
  vtkMatrix4x4 *GetRotation();

  SceneData *GetData();

  vtkLookupTable *GetLookupTable(double alpha = 0.0);
  void UpdateLookupTable();
};

#endif
