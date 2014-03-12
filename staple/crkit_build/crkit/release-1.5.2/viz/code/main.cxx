
#include <QtGui/QApplication>

#include "VizMainWindow.h"
#include "MultiVizViewerWindow.h"
#include "version.h"

/* Apple appears to define a check macro in OS X 10.5 in /usr/include/AssertMacros.h that conflicts with tclap. */
#undef check
#include <tclap/CmdLine.h>

int main( int argc, char *argv[] )
{
  std::vector<std::string> tensorFiles, rgbFiles, segFiles, probFiles, modelFiles, volFiles, fmriFiles, scalarModelFiles, sceneFiles;
  bool renderStereo = false;
  bool multiViz = false;
  bool ignoreLastDirectory = false;

  try
    {
      TCLAP::CmdLine cmd("Computational Radiology Laboratory Viewer",' ',Version::GetVersion());
      TCLAP::MultiArg<std::string> tensorFileArg("t","tensorVolume","a file containing a tensor volume",false,"tensor file name",cmd);
      TCLAP::MultiArg<std::string> rgbVolumeFileArg("r","RGBVolume","a file containing a RGB volume",false,"RGB file name",cmd);
      TCLAP::MultiArg<std::string> segmentationFileArg("s","segmentedVolume","a file containing a segmented volume",false,"segmented volume file name",cmd);
      TCLAP::MultiArg<std::string> probabilityFileArg("p","probabilityMap","a file containing a probability map",false,"probability map file name",cmd);
      TCLAP::MultiArg<std::string> modelFileArg("m","model","a file containing a polydata model or surface without scalars",false,"model file name",cmd);
      TCLAP::MultiArg<std::string> scalarModelFileArg("M","scalarModel","a file containing a polydata model or surface with scalars",false,"scalar model file name",cmd);
      TCLAP::MultiArg<std::string> fmriFileArg("f","fmriMap","a file containing an fMRI activation map",false,"fmri map file name",cmd);
      TCLAP::MultiArg<std::string> sceneFileArg("x","xmlScene","a file containing an XML scene description",false,"XML scene map file name",cmd);
      TCLAP::SwitchArg stereoArg("3","stereo","turn on stereoscopic rendering",cmd,false);
      TCLAP::SwitchArg multiVizArg("c","multiviz","use MultiViz mode",cmd,false);
      TCLAP::SwitchArg ignoreLastDirectoryArg("i","ignore","use current directory in open/save dialogs",cmd,false);
      TCLAP::UnlabeledMultiArg<std::string> volumeFilesArg("volumes","file(s) containing greyscale volume(s)",false,"volume file name(s)",cmd);

      cmd.parse(argc,argv);

      if (tensorFileArg.isSet()) tensorFiles = tensorFileArg.getValue();
      if (rgbVolumeFileArg.isSet()) rgbFiles = rgbVolumeFileArg.getValue();
      if (segmentationFileArg.isSet()) segFiles = segmentationFileArg.getValue();
      if (probabilityFileArg.isSet()) probFiles = probabilityFileArg.getValue();
      if (modelFileArg.isSet()) modelFiles = modelFileArg.getValue();
      if (scalarModelFileArg.isSet()) scalarModelFiles = scalarModelFileArg.getValue();
      if (fmriFileArg.isSet()) fmriFiles = fmriFileArg.getValue();
      if (sceneFileArg.isSet()) sceneFiles = sceneFileArg.getValue();
      renderStereo = stereoArg.getValue();
      multiViz = multiVizArg.isSet();
      ignoreLastDirectory = ignoreLastDirectoryArg.isSet();
      if (volumeFilesArg.isSet()) volFiles = volumeFilesArg.getValue();
    }
  catch (TCLAP::ArgException &e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(1);
    }

  QApplication app( argc, argv );

  if (!multiViz)
    {
      VizMainWindow w;

      w.stereo = renderStereo;
      w.ignoreLastDirectory = ignoreLastDirectory;
      w.show();

      std::vector<std::string>::iterator it;

      // open all files
      for (it = volFiles.begin(); it != volFiles.end(); it++)
	w.FileOpenWithName(QString(it->c_str()),ImageTypeDefinitions::greyscale);
      for (it = tensorFiles.begin(); it != tensorFiles.end(); it++)
	w.FileOpenWithName(QString(it->c_str()),ImageTypeDefinitions::tensor);
      for (it = rgbFiles.begin(); it != rgbFiles.end(); it ++)
	w.FileOpenWithName(QString(it->c_str()),ImageTypeDefinitions::rgb);
      for (it = segFiles.begin(); it != segFiles.end(); it++)
	w.FileOpenWithName(QString(it->c_str()),ImageTypeDefinitions::segmentation);
      for (it = probFiles.begin(); it != probFiles.end(); it++)
	w.FileOpenWithName(QString(it->c_str()),ImageTypeDefinitions::probmap);
      for (it = modelFiles.begin(); it != modelFiles.end(); it++)
	w.FileOpenWithName(QString(it->c_str()),ImageTypeDefinitions::polydata);
      for (it = scalarModelFiles.begin(); it != scalarModelFiles.end(); it++)
	w.FileOpenWithName(QString(it->c_str()),ImageTypeDefinitions::scalarpolydata);
      for (it = fmriFiles.begin(); it != fmriFiles.end(); it++)
	w.FileOpenWithName(QString(it->c_str()),ImageTypeDefinitions::fmri);
      for (it = sceneFiles.begin(); it != sceneFiles.end(); it++)
	w.OpenScene(QString(it->c_str()));

      // Initiate the event loop.
      return app.exec();

    }
  else
    {
      MultiVizViewerWindow w;
      w.show();
      for (std::vector<std::string>::iterator it = volFiles.begin(); it != volFiles.end(); it++)
	w.LoadImage(*it,ImageTypeDefinitions::greyscale);

      // Initiate the event loop.
      return app.exec();

    }
}
