
/* (c) Simon Warfield simon.warfield@childrens.harvard.edu 2006 */

/*
 * This class represents images loaded into the scene.
 * It provides basic functionality for identifying the image and allowing
 * access to the data of the image.
 *
 * Aspects related to display of the image in the scene must be handled
 * separately.
 */

#ifndef _SCENEIMAGE_INCLUDED
#define _SCENEIMAGE_INCLUDED 1

#include "ImageTypeDefinitions.h"
#include <SceneData.h>

// Qt classes
#include <QString>
#include <QColor>

class SceneImageView;

class SceneImage : public SceneData
{

  protected:
  QString name;   // identifier in the scene
  QString fileName; // Name of the file that represents the image on disk.

  public:
  // Pointer to internal image represented.
  ImageTypeDefinitions::ImageType::Pointer inputImage;
  ImageTypeDefinitions::ImageType::DirectionType inputImageDirection;

  // Details of the view generate are handled by the SceneImageView class.
  SceneImageView *sceneImageView;

  SceneImage();
  enum ColorType {noScalars, lower, upper};

  // copy constructor
  SceneImage(const SceneImage &o);

  // Construct with a specified handle (name)
  SceneImage(const QString newName);

  // Construct with a specified handle (name) and load an image
  SceneImage(QString newName, QString loadFileName);

  // destructor
  virtual ~SceneImage();

  void SetName(const char *n1);
  void SetName(QString n1);
  QString GetName() const;

  // This loads and orients all images in memory into a standard orientation.
  bool LoadImage(QString loadFileName);

  ImageTypeDefinitions::ImageType::Pointer GetImage();
  void SetImage(ImageTypeDefinitions::ImageType::Pointer image);
  ImageTypeDefinitions::ImageType::DirectionType GetImageDirection();

  SceneImageView *GetSceneImageView();

  void SetColor(QColor newColor, ColorType whichOne);
  QColor GetColor(ColorType whichOne);
  void ToggleLookupTable();
  bool GetLogLookupTable();
  void SetLowerValue(double value);
  void SetUpperValue(double value);
  double GetLowerValue();
  double GetUpperValue();

private:
  QColor lowerColor,upperColor;
  bool logLookupTable;
  double lowerValue, upperValue;
};

#endif
