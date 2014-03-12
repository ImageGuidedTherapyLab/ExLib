#ifndef _VIZVIEWERBASE_H_
#define _VIZVIEWERBASE_H_

#include <QMainWindow>

class VizMainWindow;
class SceneData;

// abstract base class for viz viewers
class VizViewerBase : public QMainWindow
{
  Q_OBJECT

public:
  virtual ~VizViewerBase() {};
  virtual void SetMarker(unsigned int _group, double* worldPos) = 0;
  virtual void SetMarkerInViewers(unsigned int _group, double *worldPos) = 0;
  virtual void CursorPositionHasChanged(double pos[3]) = 0;
  virtual unsigned int GetGroup() = 0;
  virtual unsigned int GetNumberOfRenderersForObject(SceneData *sd) = 0;
  virtual void SaveModel(SceneData *model,std::string fileName) = 0;
  virtual void RotateModel(SceneData *model) = 0;
  virtual void Render() = 0;

  VizMainWindow *mainWin;
};

#endif
