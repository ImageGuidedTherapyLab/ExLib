#ifndef _SCENETENSORS_INCLUDED
#define _SCENETENSORS_INCLUDED 1

#include "ImageTypeDefinitions.h"
#include <SceneData.h>

// Qt classes
#include <QString>

class crlTensorGlyph;
class vtkImageData;
class vtkPolyData;
class vtkExtractVOI;
class vtkPolyDataNormals;
class vtkImageChangeInformation;
class SceneImageView;
class vtkTransformPolyDataFilter;

class SceneTensors : public SceneData
{

  protected:
  QString name;   // identifier in the scene
  QString fileName; // Name of the file that represents the tensor volume on disk.

  public:
  SceneTensors();

  // copy constructor
  SceneTensors(const SceneTensors &o);

  // Construct with a specified handle (name)
  SceneTensors(const QString newName);

  // Construct with a specified handle (name) and load a tensor volume
  SceneTensors(QString newName, QString loadFileName);

  // destructor
  virtual ~SceneTensors();

  void SetName(const char *n1);
  void SetName(QString n1);
  QString GetName() const;

  void SetScale(double newScale);
  double GetScale();

  // This loads the actual tensors
  bool LoadTensors(QString loadFileName);

  int GetMinSlice();
  int GetMaxSlice();

  void SetSlice(int slice, ImageTypeDefinitions::OrientationType newOrientation);

  vtkPolyData *GetVTKTensorGlyphs();
  vtkPolyData *GetVTKTensorGlyphs(SceneImageView *siv);
  vtkImageData *GetRGBImage();
  ImageTypeDefinitions::ImageType::DirectionType GetImageDirection();

private:
  vtkImageData *sp;
  vtkPolyData *glyphs;
  vtkExtractVOI *extractVOI;
  vtkImageChangeInformation *ici;
  crlTensorGlyph *ellipsoids;
  vtkPolyDataNormals *ellipNormals;
  vtkTransformPolyDataFilter *tpdf;

  int extent[6];
  int currentSlice;

  ImageTypeDefinitions::OrientationType orientation;
  ImageTypeDefinitions::ImageType::DirectionType imageDirection;
  double translation[3];
  double scale;
};

#endif
