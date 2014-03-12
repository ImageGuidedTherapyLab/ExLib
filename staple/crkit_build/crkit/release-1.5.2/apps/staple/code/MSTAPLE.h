/*
 * Copyright (c) 2008-2011 Children's Hospital Boston.
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://www.crl.med.harvard.edu
 * For more information contact: simon.warfield@childrens.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#ifndef _MSTAPLE_h_
#define _MSTAPLE_h_

#include <itkImage.h>
#include "crlMSTAPLEImageFilter.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

class MSTAPLEBase
{
public:
  MSTAPLEBase()  { }
  virtual ~MSTAPLEBase() {}

  void SetOutputFileName( const char *s )
  { m_OutputFile = s; }

  const std::string &GetOutputFileName() const
  { return m_OutputFile; }

  void AddFileName( const char *s )
  {
    std::string tmp(s);
    m_Files.push_back(tmp);
  }

  const std::string &GetFileName( unsigned int i ) const
  { return m_Files[i]; }

  void ClearFileNames()
  { m_Files.clear(); }

  unsigned int GetNumberOfFiles() const
  { return m_Files.size(); }

  virtual void SetUseWriteCompression(bool useCompression) = 0;
  virtual bool GetUseWriteCompression( ) = 0;

  virtual std::vector<double> GetPerformance( unsigned int ) = 0;

  virtual int Execute() = 0;
  virtual unsigned int GetElapsedIterations() = 0;
  virtual void SetMaximumIterations( unsigned int ) = 0;

  virtual double GetRelativeConvergenceThreshold( ) = 0;
  virtual void SetRelativeConvergenceThreshold( double ) = 0;
  virtual void SetUnderflowProtection( int ) = 0;
  virtual int GetUnderflowProtection( ) = 0;
  virtual void SetStationaryPrior( std::vector<double> priors ) = 0;
  virtual void SetAssignConsensusVoxels(bool on) = 0;
  virtual void SetMAPStaple( bool ) = 0;
  virtual bool GetMAPStaple() = 0;
  virtual void SetMAPAlpha( double) = 0;
  virtual double GetMAPAlpha() = 0;
  virtual void SetMAPBeta( double ) = 0;
  virtual double GetMAPBeta() = 0;
  virtual void SetInitialExpertPerformanceParameters( std::vector<double> perf) = 0;
  virtual void SetStationaryPriorWeight( double ) = 0;
  virtual void SetStartAtMStep( bool ) = 0;
  virtual void SetStartAtEStep( bool ) = 0;

protected:
  std::vector< std::string > m_Files;
  std::string m_OutputFile;
};

template< unsigned int VDimension >
class MSTAPLE : public MSTAPLEBase
{

private:
  const static unsigned int ODimension = VDimension + 1;

public:

  /* The output image also stores calculations while iterating. */
  // PixelType of double is more accurate than using 'float'.
  typedef itk::Image< double , ODimension > OutputImageType;
  typedef itk::Image< unsigned short, VDimension > InputImageType;
  typedef crl::MSTAPLEImageFilter<InputImageType, OutputImageType> MSTAPLEFilterType;

  MSTAPLE()
  {
    m_MSTAPLE = MSTAPLEFilterType::New();
  }
  virtual ~MSTAPLE() {}

  /** This replaces the default calculation with specified values.
   * In most circumstances, the default calculation is appropriate.
   */
  virtual void SetStationaryPrior( std::vector<double> priors )
  {
    m_MSTAPLE->SetStationaryPrior(priors);
  }

  virtual void SetStationaryPriorWeight( double w )
  {
    m_MSTAPLE->SetStationaryPriorWeight(w);
  }

  virtual void SetInitialExpertPerformanceParameters( std::vector<double> perf)
  {
    m_MSTAPLE->SetInitialExpertPerformanceParameters( perf );
  }

  virtual std::vector<double> GetPerformance( unsigned int i )
  { return m_MSTAPLE->GetPerformance(i); }

  virtual unsigned int GetElapsedIterations()
  { return m_MSTAPLE->GetElapsedIterations(); }

  virtual void SetMaximumIterations(unsigned int maxiterations)
  { m_MSTAPLE->SetMaximumIterations(maxiterations); }

  virtual void SetRelativeConvergenceThreshold(double relconvergence)
  { m_MSTAPLE->SetRelativeConvergenceThreshold(relconvergence); }
  virtual double GetRelativeConvergenceThreshold()
  { return m_MSTAPLE->GetRelativeConvergenceThreshold(); }

  virtual void SetUnderflowProtection( int u ) {
    m_MSTAPLE->SetUnderflowProtection(u);
  };
  virtual int GetUnderflowProtection( ) {
    return m_MSTAPLE->GetUnderflowProtection();
  };

  virtual void SetUseWriteCompression(bool useCompression)
  { m_MSTAPLE->SetUseWriteCompression( useCompression ); };

  virtual bool GetUseWriteCompression( )
  { return m_MSTAPLE->GetUseWriteCompression(); }

  virtual void SetAssignConsensusVoxels(bool on)
  { m_MSTAPLE->SetAssignConsensusVoxels(on); };

  virtual void SetMAPStaple(bool map)
  { m_MSTAPLE->SetMAPStaple(map); };

  virtual bool GetMAPStaple()
  { return m_MSTAPLE->GetMAPStaple(); };

  virtual void SetMAPAlpha(double alpha)
  { m_MSTAPLE->SetMAPAlpha(alpha); };

  virtual double GetMAPAlpha()
  { return m_MSTAPLE->GetMAPAlpha(); };

  virtual void SetMAPBeta(double beta)
  { m_MSTAPLE->SetMAPBeta(beta); };

  virtual double GetMAPBeta()
  { return m_MSTAPLE->GetMAPBeta(); };

  virtual void SetStartAtMStep( bool go ) 
    { return m_MSTAPLE->SetStartAtMStep(go); };

  virtual void SetStartAtEStep( bool go ) 
    { return m_MSTAPLE->SetStartAtEStep(go); };

  virtual int Execute();

private:
  typename MSTAPLEFilterType::Pointer m_MSTAPLE;
};

#ifndef ITK_MANUAL_INSTATIATION
#include "MSTAPLE.txx"
#endif

#endif
