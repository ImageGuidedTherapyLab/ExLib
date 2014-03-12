#include "VizWizardSTAPLE.h"

#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QListWidget>
#include <QMessageBox>
#include <QLineEdit>
#include <QPushButton>
#include <QFileDialog>
#include <QCheckBox>
#include <QFileInfo>

#include <VizMainWindow.h>
#include <ImageTypeDefinitions.h>

#include <MSTAPLE.h>
#include <crlIndexOfMaxComponent.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageIOBase.h>


// These macros require certain variables to have been initialized
// by calls to imageIOBase
#define INSTANTIATECLASS(dim, scalar) do {     \
    if (imageIODimensionality == dim) {		       \
      crlIndexOfMaxComponentPtr =			\
	crl::IndexOfMaxComponent< dim, scalar >::New(); \
    }							\
  } while (0)

#define INSTANTIATECLASSOVERDIM(scalar) do {	\
    INSTANTIATECLASS(2 , scalar);		\
    INSTANTIATECLASS(3 , scalar);		\
    INSTANTIATECLASS(4 , scalar);		\
    INSTANTIATECLASS(5 , scalar);		\
    INSTANTIATECLASS(6 , scalar);		\
  } while (0)

#define INSTANTIATECLASSOVERDIMANDTYPE(componentType) do {	\
    switch(componentType) {					\
    case itk::ImageIOBase::UCHAR :				\
      INSTANTIATECLASSOVERDIM(unsigned char); break;		\
    case itk::ImageIOBase::CHAR :				\
      INSTANTIATECLASSOVERDIM(char); break;			\
    case itk::ImageIOBase::USHORT :				\
      INSTANTIATECLASSOVERDIM(unsigned short); break;		\
    case itk::ImageIOBase::SHORT :				\
      INSTANTIATECLASSOVERDIM(signed short); break;		\
    case itk::ImageIOBase::UINT :				\
      INSTANTIATECLASSOVERDIM(unsigned int); break;		\
    case itk::ImageIOBase::INT :				\
      INSTANTIATECLASSOVERDIM(signed int); break;		\
    case itk::ImageIOBase::ULONG :				\
      INSTANTIATECLASSOVERDIM(unsigned long); break;		\
    case itk::ImageIOBase::LONG :				\
      INSTANTIATECLASSOVERDIM(signed long); break;		\
    case itk::ImageIOBase::FLOAT :				\
      INSTANTIATECLASSOVERDIM(float); break;			\
    case itk::ImageIOBase::DOUBLE :				\
      INSTANTIATECLASSOVERDIM(double); break;			\
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE :		\
    default : break;						\
    }								\
  } while (0)


VizWizardSTAPLE::VizWizardSTAPLE(VizMainWindow *parent) : QWizard(parent)
{
  addPage(new VizWizardSTAPLEIntroPage);
  addPage(new VizWizardSTAPLEInputPage(this,this));
  addPage(new VizWizardSTAPLEOutputPage(this,this));
  addPage(new VizWizardSTAPLEOptionsPage);

  setWindowTitle("CRL STAPLE Wizard");
  mainWin = parent;
}

void VizWizardSTAPLE::accept()
{
  // initialize STAPLE
  MSTAPLEBase * staple;
  typedef itk::Image<unsigned char,3> ImageType;

  itk::ImageFileReader<ImageType>::Pointer reader = itk::ImageFileReader<ImageType>::New();
  reader->SetFileName(inputFiles[0]);
  reader->GenerateOutputInformation();

  itk::ImageIOBase *iob = reader->GetImageIO();
  int dimensionality = iob->GetNumberOfDimensions();

  if (iob->GetPixelType() != itk::ImageIOBase::SCALAR)
    {
      // TODO: this shouldn't happen, so add some error handling here...
    }

  switch (dimensionality)
    {
    case 1:
      staple = new MSTAPLE<1>;
      break;
    case 2:
      staple = new MSTAPLE<2>;
      break;
    case 3:
      staple = new MSTAPLE<3>;
      break;
    case 4:
      staple = new MSTAPLE<4>;
      break;
    default:
      // TODO: this also shouldn't happen...
      staple = new MSTAPLE<3>;
    }

  // file names
  std::vector<std::string>::iterator it;
  for (it = inputFiles.begin(); it != inputFiles.end(); it++)
    staple->AddFileName(it->c_str());


  staple->SetOutputFileName(outputFileWeights.c_str());

  staple->Execute();

  // now run maximum component computation
  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;
  itk::ImageIOBase::Pointer imageIO;
  unsigned int imageIODimensionality;

  crl::IndexOfMaxComponentBase::Pointer crlIndexOfMaxComponentPtr = 0;
  GetImageType(outputFileWeights.c_str(),pixelType,componentType,imageIO,imageIODimensionality);
  if (pixelType == itk::ImageIOBase::SCALAR)
    {
      INSTANTIATECLASSOVERDIMANDTYPE(componentType);
    }
  crlIndexOfMaxComponentPtr->SetInputFileName(outputFileWeights.c_str());
  crlIndexOfMaxComponentPtr->SetOutputFileName(outputFileMaxComponent.c_str());
  std::cout << "Calculating index of maximum component...";
  crlIndexOfMaxComponentPtr->Execute();

  // add output image to file list
  mainWin->FileOpenWithName(QString(outputFileMaxComponent.c_str()),ImageTypeDefinitions::segmentation);
  std::cout << "done." << std::endl;
  QDialog::accept();
}


VizWizardSTAPLEIntroPage::VizWizardSTAPLEIntroPage(QWidget *parent) : QWizardPage(parent)
{
  setTitle("CRL STAPLE");

  QLabel *label = new QLabel("This wizard will help you run STAPLE on images loaded in crlViz.");
  label->setWordWrap(true);

  QVBoxLayout *layout = new QVBoxLayout();
  layout->addWidget(label);
  setLayout(layout);

}

VizWizardSTAPLEInputPage::VizWizardSTAPLEInputPage(QWidget *parent, VizWizardSTAPLE *vws) : QWizardPage(parent)
{
  setTitle("Step 1: Input files");

  QLabel *label = new QLabel("Please drag your input files from the crLViz main window and drop them here.");
  label->setWordWrap(true);

  files = new QListWidget();

  QVBoxLayout *layout = new QVBoxLayout();
  layout->addWidget(label);
  layout->addWidget(files);

  setLayout(layout);
  setAcceptDrops(true);
  parentWizard = vws;
}

void VizWizardSTAPLEInputPage::dragEnterEvent(QDragEnterEvent *e)
{
  if (e->mimeData()->hasFormat("application/x-qabstractitemmodeldatalist")) {
    e->acceptProposedAction();
    return;
  }
  e->ignore();
}

void VizWizardSTAPLEInputPage::dropEvent(QDropEvent *e)
{
  if (!e->mimeData()->hasFormat("application/x-qabstractitemmodeldatalist")) {
    e->ignore();
    return;
  }

  QByteArray data = e->mimeData()->data("application/x-qabstractitemmodeldatalist");
  QDataStream stream(&data, QIODevice::ReadOnly);
  QMap<int, QVariant> map;

  if (!stream.atEnd())
    {
      int r,c;
      stream >> r >> c >> map;
      QString handle = map.begin().value().toString();
      ImageInformation img = parentWizard->mainWin->GetImageInformationForHandle(handle);

      bool found = false;
      for (int i = 0; i < files->count(); i++)
	{
	  if (files->item(i)->text() == handle)
	    found = true;
	}
      if (!found)
	files->addItem(handle);
      else
	QMessageBox::information(this,"Image exists","This image has already been added to the file list.");
    }

  e->acceptProposedAction();
}

bool VizWizardSTAPLEInputPage::validatePage()
{
  if (files->count() == 0)
    {
      QMessageBox::information(this,"Select files","Please select one or more input files.");
      return false;
    }
  for (int i = 0; i < files->count(); i++)
    {
      parentWizard->inputFiles.push_back(parentWizard->mainWin->GetImageInformationForHandle(files->item(i)->text()).fileName);
    }
  return true;
}


VizWizardSTAPLEOutputPage::VizWizardSTAPLEOutputPage(QWidget *parent, VizWizardSTAPLE *vws) : QWizardPage(parent)
{
  setTitle("Step 2: Output file");

  fileNameWeights = new QLineEdit();
  fileNameMaxComponent = new QLineEdit();

  QPushButton *browseButton1 = new QPushButton("Browse...");
  QPushButton *browseButton2 = new QPushButton("Browse...");

  QHBoxLayout *layout1 = new QHBoxLayout();
  QHBoxLayout *layout2 = new QHBoxLayout();
  layout1->addWidget(fileNameWeights);
  layout1->addWidget(browseButton1);
  layout2->addWidget(fileNameMaxComponent);
  layout2->addWidget(browseButton2);

  QLabel *label1 = new QLabel("Weights file:");
  QLabel *label2 = new QLabel("Maximum component file:");
  QVBoxLayout *layout = new QVBoxLayout();
  layout->addWidget(label1);
  layout->addLayout(layout1);
  layout->addWidget(label2);
  layout->addLayout(layout2);

  setLayout(layout);
  parentWizard = vws;
  connect(browseButton1,SIGNAL(clicked()),this,SLOT(on_browseButton1_clicked()));
  connect(browseButton2,SIGNAL(clicked()),this,SLOT(on_browseButton2_clicked()));
}

bool VizWizardSTAPLEOutputPage::validatePage()
{
  if (fileNameWeights->text().isEmpty() || fileNameMaxComponent->text().isEmpty())
    {
      QMessageBox::information(this,"Select file","Please select output files.");
      return false;
    }
  parentWizard->outputFileWeights = fileNameWeights->text().toStdString();
  parentWizard->outputFileMaxComponent = fileNameMaxComponent->text().toStdString();
  return true;
}

void VizWizardSTAPLEOutputPage::on_browseButton1_clicked()
{
  QString s = QFileDialog::getSaveFileName(this,"Weights output file",".","NRRD files (*.nrrd)");
  if (!s.isEmpty())
    {
      if (QFileInfo(s).suffix() != "nrrd")
	s = s+".nrrd";
      fileNameWeights->setText(s);
    }
}

void VizWizardSTAPLEOutputPage::on_browseButton2_clicked()
{
  QString s = QFileDialog::getSaveFileName(this,"Maximum component output file",".","NRRD files (*.nrrd)");
  if (!s.isEmpty())
    {
      if (QFileInfo(s).suffix() != "nrrd")
	s = s+".nrrd";
      fileNameMaxComponent->setText(s);
    }
}

VizWizardSTAPLEOptionsPage::VizWizardSTAPLEOptionsPage(QWidget *parent)
{
  setTitle("Step 3: Options");

  QCheckBox *compression = new QCheckBox("Use compression");

  QVBoxLayout *layout = new QVBoxLayout();
  layout->addWidget(compression);

  setLayout(layout);
}
