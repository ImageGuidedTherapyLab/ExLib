#ifndef _CRLDICOMEXTENDDICTIONARY_H_
#define _CRLDICOMEXTENDDICTIONARY_H_

#include <gdcmGlobal.h>
#include <gdcmDictSet.h>


void crlDICOMExtendDictionary()
{
  // Add DICOM dictionary entries for private Siemens tags
  gdcm::Dict *pubDict = gdcm::Global::GetDicts()->GetDict(PUB_DICT_NAME);
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x0019),uint16_t(0x0010),"LO","1","PrivateCreator"));
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x0019),uint16_t(0x100a),"US","1","NumberOfImagesInMosaic"));
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x0019),uint16_t(0x100b),"DS","1","SliceMeasurementDuration"));
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x0019),uint16_t(0x100c),"IS","1","B_value"));
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x0019),uint16_t(0x100d),"CS","1","DiffusionDirectionality"));
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x0019),uint16_t(0x100e),"FD","3","DiffusionGradientDirection"));
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x0019),uint16_t(0x100f),"SH","1","GradientMode"));
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x0019),uint16_t(0x1027),"FD","6","B_matrix"));
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x0019),uint16_t(0x1028),"FD","1","BandwidthPerPixelPhaseEncode"));

  // Add DICOM dictionary entries for private GE tags
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x0043),uint16_t(0x102f),"SS","1","Image Type"));

  // Add DICOM dictionary entries for private Philips tags
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x2005),uint16_t(0x10b0),"DS","1","DiffusionGradient1"));
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x2005),uint16_t(0x10b1),"DS","1","DiffusionGradient2"));
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x2005),uint16_t(0x10b2),"DS","1","DiffusionGradient3"));
  pubDict->AddEntry(gdcm::DictEntry(uint16_t(0x2001),uint16_t(0x1003),"IS","1","B_valueP"));
}

#endif
