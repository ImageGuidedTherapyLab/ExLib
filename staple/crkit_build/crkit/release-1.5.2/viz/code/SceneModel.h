
#ifndef _SCENEMODEL_INCLUDED
#define _SCENEMODEL_INCLUDED 1

#include "ImageTypeDefinitions.h"
#include <SceneData.h>

#include <vector>

// Qt classes
#include <QString>
#include <QColor>

class vtkPolyData;
class vtkActor;
class vtkImageData;
class SceneImageView;

class SceneModel : public SceneData
{

  protected:
  QString name;   // identifier in the scene
  QString fileName; // Name of the file that represents the model on disk.

  public:
  SceneModel();
  enum ColorType {noScalars, lower, upper};

  // copy constructor
  SceneModel(const SceneModel &o);

  // Construct with a specified handle (name)
  SceneModel(const QString newName);

  // Construct with a specified handle (name) and load a model
  SceneModel(QString newName, QString loadFileName);

  // destructor
  virtual ~SceneModel();

  void SetName(const char *n1);
  void SetName(QString n1);
  QString GetName() const;

  // This loads the actual model
  bool LoadModel(QString loadFileName, bool hasScalars = false);

  vtkPolyData *GetPolyData();
  void SetPolyData(vtkPolyData *data);

  void SetColor(QColor newColor, ColorType whichOne);
  QColor GetColor(ColorType whichOne);
  void SetLowerValue(double value);
  void SetUpperValue(double value);
  double GetLowerValue();
  double GetUpperValue();
  void SetOpacity(double value);
  double GetOpacity();

  vtkActor *CreateNewActor();
  vtkActor *CreateNewActor(SceneImageView *siv);
  void AddTextureMap(vtkImageData* texture);

  bool isFlippedX, isFlippedY, isFlippedZ;
  bool ignoreScalars;

  void ToggleLookupTable();
  bool GetLogLookupTable();

private:
  vtkPolyData *polyData;
  QColor color,lowerColor,upperColor;
  double lowerValue, upperValue;
  double opacity;
  std::vector<vtkActor*> modelActors;
  bool logLookupTable;
};

#endif
