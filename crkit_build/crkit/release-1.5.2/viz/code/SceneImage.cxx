#include <SceneImage.h>
#include <SceneImageView.h>

// ITK classes
#include <itkOrientImageFilter.h>

SceneImage::SceneImage()
{
	inputImage = 0;
	sceneImageView = 0;
	logLookupTable = false;
	lowerColor.setHsv(180,255,255,255);
	upperColor.setHsvF(0.0,1.0,1.0,1.0);
	lowerValue = 0;
	upperValue = 0;
}

// copy constructor
SceneImage::SceneImage(const SceneImage &o) {
	name = o.GetName();
	inputImage = o.inputImage;
	sceneImageView = o.sceneImageView;
	logLookupTable = o.logLookupTable;
	lowerColor = o.lowerColor;
	upperColor = o.upperColor;
	lowerValue = o.lowerValue;
	upperValue = o.upperValue;
}

// Construct with a specified handle (name)
SceneImage::SceneImage(const QString newName)
{
	name = newName;
	inputImage = 0;
	sceneImageView = 0;
	logLookupTable = false;
	lowerColor.setHsv(180,255,255,255);
	upperColor.setHsvF(0.0,1.0,1.0,1.0);
	lowerValue = 0;
	upperValue = 0;
}

// Construct with a specified handle (name) and load an image
SceneImage::SceneImage(QString newName, QString loadFileName)
{
	name = newName;
	inputImage = 0;
	sceneImageView = 0;
	if (!LoadImage(loadFileName)) {
		fileName = "";
	}
	logLookupTable = false;
	lowerColor.setHsv(180,255,255,255);
	upperColor.setHsvF(0.0,1.0,1.0,1.0);
	lowerValue = 0;
	upperValue = 0;
}

// destructor
SceneImage::~SceneImage() {
	inputImage = 0;
	if (sceneImageView)
		delete sceneImageView;
	sceneImageView = 0;
}

void SceneImage::SetName(const char *n1)
{
	name = QString(n1);
}

void SceneImage::SetName(QString n1)
{
	name = n1;
}

QString SceneImage::GetName() const
{
	return name;
}

// This loads and orients all images in memory into a standard orientation.
bool SceneImage::LoadImage(QString loadFileName) {
	fileName = loadFileName;

	// Declare and initialize an image orienter filter
	itk::OrientImageFilter< ImageTypeDefinitions::ImageType,
		ImageTypeDefinitions::ImageType >::Pointer imageOrienter;
	imageOrienter = itk::OrientImageFilter< ImageTypeDefinitions::ImageType,
		ImageTypeDefinitions::ImageType >::New();

	// Now load the image off disk and initialize the image data.
	ImageTypeDefinitions::ReaderType::Pointer reader =
		ImageTypeDefinitions::ReaderType::New();

	reader->SetFileName( loadFileName.toStdString().c_str() );
	try {
		reader->Update();
	} catch (itk::ExceptionObject & err) {
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return false;
	}
	// store original direction of the input image
	//inputImageDirection = reader->GetOutput()->GetDirection();

	imageOrienter->UseImageDirectionOn();
	imageOrienter->SetDesiredCoordinateOrientationToAxial();
	imageOrienter->SetInput( reader->GetOutput() );
	inputImage = imageOrienter->GetOutput();
	try {
		inputImage->Update();
	} catch (itk::ExceptionObject & err) {
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return false;
	}

	// Disconnect the DataObject inputImage from the pipeline.
	inputImage->DisconnectPipeline();

	// store direction of the reoriented image
	inputImageDirection = inputImage->GetDirection();

	// The smart pointers initialized using the ::New() object factory
	// will automatically UnRegister() and Delete() themselves when they go
	// out of scope.

	return true; // success
}

ImageTypeDefinitions::ImageType::Pointer SceneImage::GetImage() {
	return inputImage;
}

void SceneImage::SetImage(ImageTypeDefinitions::ImageType::Pointer image)
{
  inputImage = image;
}

ImageTypeDefinitions::ImageType::DirectionType SceneImage::GetImageDirection()
{
  return inputImageDirection;
}

SceneImageView *SceneImage::GetSceneImageView() {
	// The scene image provides a class that manages view properties
	// for the renderer.

	if (!sceneImageView) {
		sceneImageView = new SceneImageView(this);
		assert(sceneImageView != 0);
	}
	lowerValue = sceneImageView->GetVTKImageData()->GetScalarRange()[0];
	upperValue = sceneImageView->GetVTKImageData()->GetScalarRange()[1];
	return sceneImageView;
}

void SceneImage::SetColor(QColor newColor, ColorType whichOne)
{
  if (whichOne == lower)
    lowerColor = newColor;
  else
    upperColor = newColor;

  if (sceneImageView)
    sceneImageView->UpdateLookupTable();
}

QColor SceneImage::GetColor(ColorType whichOne)
{
  if (whichOne == lower)
    return lowerColor;
  else
    return upperColor;
}

void SceneImage::ToggleLookupTable()
{
  logLookupTable = !logLookupTable;
}

bool SceneImage::GetLogLookupTable()
{
  return logLookupTable;
}

double SceneImage::GetLowerValue()
{
  return lowerValue;
}

double SceneImage::GetUpperValue()
{
  return upperValue;
}

void SceneImage::SetLowerValue(double value)
{
  lowerValue = value;
  if (sceneImageView)
    sceneImageView->UpdateLookupTable();
}

void SceneImage::SetUpperValue(double value)
{
  upperValue = value;
  if (sceneImageView)
    sceneImageView->UpdateLookupTable();
}
