
#ifndef _version_h
#define _version_h

// Pull in compile time configuration options
#include "configuration.h"

// These macro antics pull in the version defined at compile time,
// and place the value of the constant into a quoted string literal.
#define CRKIT_VERSION_TO_STRING(x) CRKIT_VERSION_TO_STRING0(x)
#define CRKIT_VERSION_TO_STRING0(x) #x
#define CRKIT_VERSION CRKIT_VERSION_TO_STRING(CRKIT_VERSION_MAJOR) "." \
                    CRKIT_VERSION_TO_STRING(CRKIT_VERSION_MINOR) "." \
                    CRKIT_VERSION_TO_STRING(CRKIT_VERSION_PATCH)

#define CRKIT_SOURCE_VERSION "CRkit version " CRKIT_VERSION ", CRkit source $Revision: 832 $, $Date: 2007-07-17 10:55:45 -0400 (Tue, 17 Jul 2007) $ (GMT)"

class Version
{
  public: 
  static const char *GetVersion() { return CRKIT_VERSION_STRING; };
  static int GetMajorVersion() { return CRKIT_VERSION_MAJOR; };
  static int GetMinorVersion() { return CRKIT_VERSION_MINOR; };
  static int GetBuildVersion() { return CRKIT_VERSION_PATCH; };
  static const char *GetSourceVersion() { return CRKIT_SOURCE_VERSION; };

  protected:
  Version();

  ~Version();

  private:
  Version(const Version &); // purposely not implemented
  void operator=(const Version &); // purposely not implemented

};


#endif

