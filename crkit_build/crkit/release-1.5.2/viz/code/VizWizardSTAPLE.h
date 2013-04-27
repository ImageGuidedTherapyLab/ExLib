#ifndef _VIZWIZARDSTAPLE_H_
#define _VIZWIZARDSTAPLE_H_

#include <QWizard>
#include <QWizardPage>

class VizMainWindow;
class QListWidget;
class QLineEdit;

class VizWizardSTAPLE : public QWizard
{
  Q_OBJECT

public:
  VizWizardSTAPLE(VizMainWindow *parent);

  void accept();

  std::vector<std::string> inputFiles;
  std::string outputFileWeights;
  std::string outputFileMaxComponent;

  VizMainWindow *mainWin;
};


class VizWizardSTAPLEIntroPage : public QWizardPage
{
  Q_OBJECT

public:
  VizWizardSTAPLEIntroPage(QWidget *parent = 0);
};


class VizWizardSTAPLEInputPage : public QWizardPage
{
  Q_OBJECT

public:
  VizWizardSTAPLEInputPage(QWidget *parent, VizWizardSTAPLE *vws);
  virtual bool validatePage();

private:
  void dragEnterEvent(QDragEnterEvent *e);
  void dropEvent(QDropEvent *e);

  QListWidget *files;
  VizWizardSTAPLE *parentWizard;
};


class VizWizardSTAPLEOutputPage : public QWizardPage
{
  Q_OBJECT

public:
  VizWizardSTAPLEOutputPage(QWidget *parent, VizWizardSTAPLE *vws);
  virtual bool validatePage();

private slots:
  void on_browseButton1_clicked();
  void on_browseButton2_clicked();

private:
  QLineEdit *fileNameWeights;
  QLineEdit *fileNameMaxComponent;
  VizWizardSTAPLE *parentWizard;
};


class VizWizardSTAPLEOptionsPage : public QWizardPage
{
  Q_OBJECT

public:
  VizWizardSTAPLEOptionsPage(QWidget *parent = 0);
};

#endif
