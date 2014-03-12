#ifndef __crlTensorGlyph_h
#define __crlTensorGlyph_h

#include <vtkTensorGlyph.h>

class crlTensorGlyph : public vtkTensorGlyph
{
public:
  vtkTypeRevisionMacro(crlTensorGlyph,vtkTensorGlyph);

  // Description
  // Construct object with scaling on and scale factor 1.0. Eigenvalues are
  // extracted, glyphs are colored with input scalar data, and logarithmic
  // scaling is turned off.
  static crlTensorGlyph *New();

protected:
  crlTensorGlyph();
  ~crlTensorGlyph();

  virtual int RequestData(vtkInformation*,vtkInformationVector**,vtkInformationVector*);

private:
  crlTensorGlyph(const crlTensorGlyph&);  // Not implemented.
  void operator=(const crlTensorGlyph&);  // Not implemented.
};

#endif
