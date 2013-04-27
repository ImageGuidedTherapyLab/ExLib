#ifndef _SCENEDATA_H_
#define _SCENEDATA_H_

#include <QString>

class SceneData
{
public:
  SceneData() {};
  virtual ~SceneData() {};

  virtual QString GetName() const = 0;
};

#endif
