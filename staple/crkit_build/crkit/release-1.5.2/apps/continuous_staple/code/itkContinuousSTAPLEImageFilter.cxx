#include "itkContinuousSTAPLEImageFilter.h"

#include "itkTimeProbe.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"

// This brings in types needed for modifying the pipeline.
#include "itkProcessObject.h"

namespace itk
{
  void ContinuousSTAPLEImageFilter::InitializeExpertBias(double biasValue)
  {
    m_Bias.clear();
    m_OldBias.clear();

    BiasType tmpVal(m_NDim,biasValue);

    for (unsigned int i = 0; i < this->GetNumberOfInputs();++i)
    {
      m_Bias.push_back(tmpVal);
      m_OldBias.push_back(tmpVal);
    }
  }

  void ContinuousSTAPLEImageFilter::InitializeExpertCovariance(double covValue)
  {
    m_Covariance.clear();
    m_InvCovariances.clear();
    m_OldCovariance.clear();

    CovMatrixType tmpVal(m_NDim,m_NDim);
    for (unsigned int i = 0; i < m_NDim; ++i)
      for (unsigned int j = 0; j < m_NDim; ++j)
      {
        if (i == j)
	  tmpVal(i,j) = covValue;
	else
	  tmpVal(i,j) = 0;
      }

    for (unsigned int i = 0; i < this->GetNumberOfInputs();++i)
    {
      m_Covariance.push_back(tmpVal);
      m_InvCovariances.push_back(tmpVal);
      m_OldCovariance.push_back(tmpVal);
    }
  }
  
  void ContinuousSTAPLEImageFilter::CreateComputationMaskFromInputs()
  {
    if (m_ComputationMask)
      return;

    if (m_Verbose)
      std::cout << "No computation mask specified... Creating a default one..." << std::endl;

    m_ComputationMask = MaskImageType::New();
    m_ComputationMask->Initialize();
    InputImageRegionType region = this->GetInput(0)->GetLargestPossibleRegion();
    m_ComputationMask->SetRegions(region);
    m_ComputationMask->SetSpacing (this->GetInput(0)->GetSpacing());
    m_ComputationMask->SetOrigin (this->GetInput(0)->GetOrigin());
    m_ComputationMask->SetDirection (this->GetInput(0)->GetDirection());
    m_ComputationMask->Allocate();
    
    m_ComputationMask->FillBuffer(1);
    
    typedef ImageRegionIteratorWithIndex< MaskImageType > MaskRegionIteratorType;
    typedef ImageRegionConstIteratorWithIndex< TInputImage > InIteratorType;

    MaskRegionIteratorType maskItr(m_ComputationMask,m_ComputationMask->GetLargestPossibleRegion());
    InputPixelType tmpVec;
    
    for (unsigned int i = 0;i < this->GetNumberOfInputs();++i)
    {
      maskItr.GoToBegin();

      InIteratorType inputItr(this->GetInput(i), region);
      inputItr.GoToBegin();

      while (!maskItr.IsAtEnd())
      {
        if (maskItr.Get() == 0)
        {
          ++maskItr;
          ++inputItr;
          continue;
        }

        tmpVec = inputItr.Get();

        bool isNull = true;

        for (unsigned int j = 0;j < m_NDim;++j)
        {
          if (tmpVec[j] != 0)
          {
            isNull = false;
            break;
          }
        }

        if (isNull)
          maskItr.Set(0);

        ++maskItr;
        ++inputItr;
      }
    }
  }

  void ContinuousSTAPLEImageFilter::InitializeComputationRegionFromMask()
  {
    if (!m_ComputationMask)
      this->CreateComputationMaskFromInputs();

    typedef ImageRegionIteratorWithIndex< MaskImageType > MaskRegionIteratorType;

    MaskRegionIteratorType maskItr(m_ComputationMask,m_ComputationMask->GetLargestPossibleRegion());
    maskItr.GoToBegin();

    MaskImageType::IndexType minPos, maxPos;

    for (unsigned int i = 0;i < m_ComputationMask->GetImageDimension();++i)
    {
      minPos[i] = m_ComputationMask->GetLargestPossibleRegion().GetIndex()[i] + m_ComputationMask->GetLargestPossibleRegion().GetSize()[i];
      maxPos[i] = 0;
    }

    while (!maskItr.IsAtEnd())
    {
      if (maskItr.Get() != 0)
      {
        MaskImageType::IndexType tmpInd = maskItr.GetIndex();

        for (unsigned int i = 0;i < m_ComputationMask->GetImageDimension();++i)
        {
          if (minPos[i] > tmpInd[i])
            minPos[i] = tmpInd[i];

          if (maxPos[i] < tmpInd[i])
            maxPos[i] = tmpInd[i];
        }
      }

      ++maskItr;
    }

    m_ComputationRegion.SetIndex(minPos);

    MaskImageType::SizeType tmpSize;
    for (unsigned int i = 0;i < m_ComputationMask->GetImageDimension();++i)
      tmpSize[i] = maxPos[i] - minPos[i] + 1;

    m_ComputationRegion.SetSize(tmpSize);

    this->GetOutput()->SetRequestedRegion(m_ComputationRegion);
  }
  
  void ContinuousSTAPLEImageFilter::InitializeSplitRegionsFromMask()
  {
    if (!m_ComputationMask)
      this->CreateComputationMaskFromInputs();
    
    typedef ImageRegionIteratorWithIndex< MaskImageType > MaskRegionIteratorType;
    
    MaskRegionIteratorType maskItr(m_ComputationMask,this->GetOutput()->GetRequestedRegion());
    maskItr.GoToBegin();
    
    bool isImage3D = (this->GetOutput()->GetRequestedRegion().GetSize()[2] > 1);
    
    if (this->GetNumberOfThreads() == 1)
    {
      m_lowerLimits.clear();
      m_upperLimits.clear();
      
      m_lowerLimits.push_back(this->GetOutput()->GetRequestedRegion().GetIndex()[1 + isImage3D]);
      m_upperLimits.push_back(this->GetOutput()->GetRequestedRegion().GetSize()[1 + isImage3D] + this->GetOutput()->GetRequestedRegion().GetIndex()[1 + isImage3D] - 1);
      
      return;
    }
    
    if ((unsigned int)this->GetNumberOfThreads() > this->GetOutput()->GetRequestedRegion().GetSize()[1 + isImage3D])
    {
      m_lowerLimits.clear();
      m_upperLimits.clear();
      
      for (unsigned int i = 0;i < this->GetOutput()->GetRequestedRegion().GetSize()[1 + isImage3D];++i)
      {
        m_lowerLimits.push_back(i);
        m_upperLimits.push_back(i);
      }
      
      return;
    }
    
    unsigned int nbPts = 0;
    while (!maskItr.IsAtEnd())
    {
      if (maskItr.Get() != 0)
        nbPts++;
      ++maskItr;
    }
    
    unsigned int approxNumPtsPerThread = (unsigned int)floor(nbPts/((double)this->GetNumberOfThreads()));
    
    //std::cout << this->GetOutput()->GetRequestedRegion() << std::endl;
    //std::cout << approxNumPtsPerThread << std::endl;
    
    maskItr.GoToBegin();
    
    m_lowerLimits.clear();
    m_upperLimits.clear();
    
    m_lowerLimits.push_back(this->GetOutput()->GetRequestedRegion().GetIndex()[1 + isImage3D]);
    unsigned int fromPrevious = 0;
    
    for (unsigned int j = 0;j < (unsigned int)this->GetNumberOfThreads();++j)
    {
      unsigned int nbPtsThread = fromPrevious;
      unsigned int borneSup = m_lowerLimits[j];
      unsigned int lastSliceNbPts = 0;
      if (fromPrevious != 0)
        borneSup++;
      
      while ((nbPtsThread < approxNumPtsPerThread)&&(!maskItr.IsAtEnd()))
      {
        lastSliceNbPts = 0;
				unsigned int tmpVal = maskItr.GetIndex()[1 + isImage3D];
        // Add one slice
        while ((tmpVal == borneSup)&&(!maskItr.IsAtEnd()))
        {
          if (maskItr.Get() != 0)
            lastSliceNbPts++;

          ++maskItr;
					tmpVal = maskItr.GetIndex()[1 + isImage3D];
        }
        nbPtsThread += lastSliceNbPts;
        borneSup++;
      }
      
      if (j != (unsigned int)this->GetNumberOfThreads() - 1)
      {
        unsigned int nbPtsSliceBefore = nbPtsThread - lastSliceNbPts;
      
        if (approxNumPtsPerThread - nbPtsSliceBefore < nbPtsThread - approxNumPtsPerThread)
        {
          borneSup--;
          fromPrevious = lastSliceNbPts;
          nbPtsThread = nbPtsThread - lastSliceNbPts;
        }
        else
          fromPrevious = 0;
      }

      m_upperLimits.push_back(borneSup - 1);
      m_lowerLimits.push_back(borneSup);
      //std::cout << m_lowerLimits[j] << " " << m_upperLimits[j] << " " << nbPtsThread << std::endl;
    }
    m_lowerLimits.pop_back();
  }
  
  int ContinuousSTAPLEImageFilter::SplitRequestedRegion(int i, int num, OutputImageRegionType& splitRegion)
  {
    splitRegion = this->GetOutput()->GetRequestedRegion();
    bool isImage3D = (this->GetOutput()->GetRequestedRegion().GetSize()[2] > 1);
    
    splitRegion.SetIndex(1 + isImage3D,m_lowerLimits[i]);
    splitRegion.SetSize(1 + isImage3D,m_upperLimits[i] - m_lowerLimits[i] + 1);
    
    return m_upperLimits.size();
  }

  void ContinuousSTAPLEImageFilter::GenerateData()
  {
    /* This function is responsible for allocating the output.
     * ThreadedGenerateData is not.
     */
    if (!m_ComputationMask)
      this->CreateComputationMaskFromInputs();
    
    if (!m_OverrideComputationRegion)
      this->InitializeComputationRegionFromMask();
    else
    {
      m_ComputationRegion = this->GetInput(0)->GetLargestPossibleRegion();
      this->GetOutput()->SetRequestedRegion(m_ComputationRegion);      
    }
    
    this->InitializeSplitRegionsFromMask();
    
    try 
    {
      this->GetOutput()->SetNumberOfComponentsPerPixel (m_NDim);
      this->AllocateOutputs();
    }
    catch (std::bad_alloc &exp)
    {
      std::cerr << "std::bad_alloc exception caught !" << std::endl;
      std::cerr << "too much memory needed..." << std::endl;
      exit(1);
    }
        
    if (m_Bias.size() != this->GetNumberOfInputs())
      this->InitializeExpertBias(0.0);

    if (m_Covariance.size() != this->GetNumberOfInputs())
      this->InitializeExpertCovariance(1.0);

    m_ReferenceCovariance.set_size(m_NDim,m_NDim);
    
    unsigned int itncount = 0;
    bool continueLoop = true;
    while ((itncount < m_MaximumIterations)&&(continueLoop))
    {
      if (m_Verbose)
        std::cout << "Iteration " << itncount + 1 << "..." << std::endl;
      EstimateReferenceStandard();
      EstimatePerformanceParameters();
      
      ++itncount;

      if (itncount != 1)
        continueLoop = !endConditionReached();

      if (continueLoop)
      {
        for (unsigned int i = 0;i < this->GetNumberOfInputs();++i)
	{
	  for (unsigned int j = 0;j < m_NDim;++j)
	  {
	    m_OldBias[i][j] = m_Bias[i][j];
	    for (unsigned int k = 0;k < m_NDim;++k)
	      m_OldCovariance[i](j,k) = m_Covariance[i](j,k);
	  }
	}
      }
    }

    m_ElapsedIterations = itncount;
  }  

  void ContinuousSTAPLEImageFilter::EstimateReferenceStandard()
  {
    itk::TimeProbe tmpTime;
    tmpTime.Start();
    // Update common part
    for (unsigned int i = 0;i < m_NDim;++i)
      for (unsigned int j = 0;j < m_NDim;++j)
        m_ReferenceCovariance(i,j) = 0;

    unsigned int nbInputs = this->GetNumberOfInputs();

    // Computes inverse covariances
    for (unsigned int k = 0;k < nbInputs;++k)
    {
      vnl_matrix_inverse<double> tmpInv(m_Covariance[k]);

      m_InvCovariances[k] = tmpInv.inverse();

      for (unsigned int i = 0;i < m_NDim;++i)
        for (unsigned int j = i;j < m_NDim;++j)
          m_ReferenceCovariance(i,j) += m_InvCovariances[k](i,j);
    }

    for (unsigned int i = 0;i < m_NDim;++i)
      for (unsigned int j = i+1;j < m_NDim;++j)
        m_ReferenceCovariance(j,i) = m_ReferenceCovariance(i,j);    
    
    // Computes reference standard covariance
    vnl_matrix_inverse<double> tmpInv(m_ReferenceCovariance);
    CovMatrixType tmpMat(m_NDim,m_NDim);
    tmpMat = tmpInv.inverse();
    for (unsigned int i = 0;i < m_NDim;++i)
      for (unsigned int j = i;j < m_NDim;++j)
    {
      m_ReferenceCovariance(i,j) = tmpMat(i,j);
      if (i != j)
        m_ReferenceCovariance(j,i) = m_ReferenceCovariance(i,j);
    }
    
    // And now do the multithreaded step
    itk::MultiThreader::Pointer threaderEstep = itk::MultiThreader::New();
    
    EMStepThreadStruct *tmpStr = new EMStepThreadStruct;
    tmpStr->Filter = this;

    threaderEstep->SetNumberOfThreads(this->GetNumberOfThreads());
    threaderEstep->SetSingleMethod(this->ThreadEstimateRefStd,tmpStr);
    threaderEstep->SingleMethodExecute();

    delete tmpStr;
    
    tmpTime.Stop();
    if (m_Verbose)
      std::cout << "Reference standard estimated... " << tmpTime.GetTotal() << std::endl;
  }

  ITK_THREAD_RETURN_TYPE ContinuousSTAPLEImageFilter::ThreadEstimateRefStd(void *arg)
  {
    MultiThreader::ThreadInfoStruct *threadArgs = (MultiThreader::ThreadInfoStruct *)arg;

    unsigned int nbThread = threadArgs->ThreadID;
    unsigned int nbProcs = threadArgs->NumberOfThreads;

    EMStepThreadStruct *tmpStr = (EMStepThreadStruct *)threadArgs->UserData;
    InputImageRegionType threadRegion;
    int total = tmpStr->Filter->SplitRequestedRegion(nbThread,nbProcs,threadRegion);

    if ((int)nbThread < total)
      tmpStr->Filter->EstimateReferenceStandard(threadRegion);

    return NULL;
  }


  void ContinuousSTAPLEImageFilter::EstimateReferenceStandard(TOutputImage::RegionType &region)
  {
    typedef ImageRegionConstIteratorWithIndex< TInputImage > InIteratorType;
    typedef ImageRegionIteratorWithIndex< TOutputImage > OutRegionIteratorType;

    typedef ImageRegionIteratorWithIndex< MaskImageType > MaskRegionIteratorType;
    typedef MaskImageType::RegionType MaskRegionType;

    unsigned int nbInputs = this->GetNumberOfInputs();
    
    OutRegionIteratorType outItr(this->GetOutput(), region);
    MaskRegionIteratorType maskItr(m_ComputationMask, region);

    outItr.GoToBegin();
    maskItr.GoToBegin();
    std::vector < InIteratorType > inputIterators;
    for (unsigned int i = 0; i < nbInputs;++i)
      inputIterators.push_back(InIteratorType(this->GetInput(i), region));

    std::vector <double> tmpVal(m_NDim,0), tmpComp(m_NDim,0);
    OutputPixelType resVec;
    resVec.SetSize(m_NDim);

    while (!outItr.IsAtEnd())
    {
      if (maskItr.Get() == 0)
      {
        ++outItr;
	for (unsigned int i = 0; i < nbInputs;++i)
	  ++inputIterators[i];
	++maskItr;
	continue;
      }

      for (unsigned int j = 0;j < m_NDim;++j)
        tmpVal[j] = 0;

      for (unsigned int j = 0;j < nbInputs;++j)
      {
        for (unsigned int k = 0;k < m_NDim;++k)
	  tmpComp[k] = inputIterators[j].Get()[k] - m_Bias[j][k];

	for (unsigned int k = 0;k < m_NDim;++k)
	{
	  for (unsigned int l = 0;l < m_NDim;++l)
            tmpVal[k] += m_InvCovariances[j](k,l)*tmpComp[l];
	}
      }

      for (unsigned int k = 0;k < m_NDim;++k)
      {
	resVec[k] = 0;
	for (unsigned int l = 0;l < m_NDim;++l)
  	  resVec[k] +=  m_ReferenceCovariance(k,l)*tmpVal[l];
      }

      outItr.Set(resVec);

      ++outItr;
      ++maskItr;
      for (unsigned int i = 0; i < nbInputs;++i)
	++inputIterators[i];
    }
  }  

  void ContinuousSTAPLEImageFilter::EstimatePerformanceParameters()
  {
    itk::TimeProbe tmpTime;
    tmpTime.Start();
    itk::MultiThreader::Pointer threaderMstep = itk::MultiThreader::New();
    
    EMStepThreadStruct *tmpStr = new EMStepThreadStruct;
    tmpStr->Filter = this;

    unsigned int actualNumberOfThreads = MIN((unsigned int)this->GetNumberOfThreads(),this->GetNumberOfInputs());

    threaderMstep->SetNumberOfThreads(actualNumberOfThreads);
    threaderMstep->SetSingleMethod(this->ThreadEstimatePerfParams,tmpStr);
    threaderMstep->SingleMethodExecute();

    delete tmpStr;

    tmpTime.Stop();
    if (m_Verbose)
      std::cout << "Performance parameters estimated... " << tmpTime.GetTotal() << std::endl;
  }

  ITK_THREAD_RETURN_TYPE ContinuousSTAPLEImageFilter::ThreadEstimatePerfParams(void *arg)
  {
    MultiThreader::ThreadInfoStruct *threadArgs = (MultiThreader::ThreadInfoStruct *)arg;

    unsigned int nbThread = threadArgs->ThreadID;
    unsigned int nbProcs = threadArgs->NumberOfThreads;

    EMStepThreadStruct *tmpStr = (EMStepThreadStruct *)threadArgs->UserData;
    unsigned int nbExperts = tmpStr->Filter->GetNumberOfInputs();

    unsigned int minExp = (unsigned int)floor((double)nbThread*nbExperts/nbProcs);
    unsigned int maxExp = (unsigned int)floor((double)(nbThread + 1.0)*nbExperts/nbProcs);

    maxExp = MIN(nbExperts,maxExp);

    tmpStr->Filter->EstimatePerformanceParameters(minExp,maxExp);

    return NULL;
  }

  void ContinuousSTAPLEImageFilter::EstimatePerformanceParameters(unsigned int minExp, unsigned int maxExp)
  {
    typedef ImageRegionConstIteratorWithIndex< TInputImage > InIteratorType;
    typedef ImageRegionIteratorWithIndex< TOutputImage > OutRegionIteratorType;

    typedef ImageRegionIteratorWithIndex< MaskImageType > MaskRegionIteratorType;
    typedef MaskImageType::RegionType MaskRegionType;

    std::vector < std::vector <double> > tmpComp;
    
    InputPixelType tmpData;
    OutputPixelType outputPoint;

    OutputImagePointer W = this->GetOutput();
    OutRegionIteratorType outItr(W, m_ComputationRegion);

    MaskRegionIteratorType maskItr(m_ComputationMask, m_ComputationRegion);
    
    std::vector < InIteratorType > inputIterators;
    for (unsigned int i = 0; i < maxExp;++i)
      inputIterators.push_back(InIteratorType(this->GetInput(i), m_ComputationRegion));
    
    unsigned int nbPts = 0;
    std::vector <double> tmpVec(m_NDim,0);
    
    tmpComp.clear();
    for (unsigned int i = 0;i < maxExp;++i)
      tmpComp.push_back(tmpVec);

    // Compute bias parameters
    while (!outItr.IsAtEnd())
    {
      if (maskItr.Get() == 0)
      {
        ++outItr;
        for (unsigned int i = minExp;i < maxExp;++i)
          ++inputIterators[i];
        ++maskItr;
        continue;
      }
      
      nbPts++;
      outputPoint = outItr.Get();
      for (unsigned int i = minExp;i < maxExp;++i)
      {
        tmpData = inputIterators[i].Get();

        for (unsigned int k = 0;k < m_NDim;++k)
          tmpComp[i][k] += tmpData[k] - outputPoint[k];
      }
      
      ++outItr;
      ++maskItr;
      for (unsigned int i = minExp;i < maxExp;++i)
        ++inputIterators[i];
    }
    
    for (unsigned int i = minExp;i < maxExp;++i)
      for (unsigned int k = 0;k < m_NDim;++k)
        m_Bias[i][k] = tmpComp[i][k]/nbPts;

    // Computes covariance parameters
    std::vector < CovMatrixType > tmpVal;
    for (unsigned int i = 0;i < maxExp;++i)
    {
      CovMatrixType tmpMat(m_NDim,m_NDim,0);
      tmpVal.push_back(tmpMat);
    }

    outItr.GoToBegin();
    maskItr.GoToBegin();
    for (unsigned int i = minExp;i < maxExp;++i)
      inputIterators[i].GoToBegin();
    
    while (!outItr.IsAtEnd())
    {
      if (maskItr.Get() == 0)
      {
        ++outItr;
        for (unsigned int i = minExp;i < maxExp;++i)
          ++inputIterators[i];
        ++maskItr;
        continue;
      }
      
      outputPoint = outItr.Get();
      for (unsigned int i = minExp;i < maxExp;++i)
      {
        tmpData = inputIterators[i].Get();
        
        for (unsigned int k = 0;k < m_NDim;++k)
        {
          for (unsigned int l = k;l < m_NDim;++l)
          {
            tmpVal[i](k,l) += (tmpData[k] - m_Bias[i][k] - outputPoint[k])*(tmpData[l] - m_Bias[i][l] - outputPoint[l]);
          }
        }
      }
      
      ++outItr;
      ++maskItr;
      for (unsigned int i = minExp;i < maxExp;++i)
        ++inputIterators[i];
    }
    
    for (unsigned int i = minExp;i < maxExp;++i)
      for (unsigned int k = 0;k < m_NDim;++k)
        for (unsigned int l = k;l < m_NDim;++l)
        {
          m_Covariance[i](k,l) = tmpVal[i](k,l)/nbPts + m_ReferenceCovariance(k,l);
          m_Covariance[i](l,k) = m_Covariance[i](k,l);
        }
  }

  void ContinuousSTAPLEImageFilter::PrintSelf(std::ostream& os, itk::Indent indent) const
  {
    os << indent << "covParams = zeros(" << this->GetNumberOfInputs() << "," << m_NDim << "," << m_NDim << ");" << std::endl;
    for (unsigned int i = 0;i < this->GetNumberOfInputs();++i)
    {
      os << indent << "biasParams(" << i+1 << ",:) = [";
      for (unsigned int j = 0;j < m_NDim;++j)
	os << indent << m_Bias[i][j] << " ";

      os << indent << "];" << std::endl;
      
      os << indent << "covParams(" << i+1 << ",:,:) = [";
      for (unsigned int j = 0;j < m_NDim;++j)
      {
	for (unsigned int k = 0;k < m_NDim;++k)
	  os << indent << m_Covariance[i](j,k) << " ";

	if (j == m_NDim - 1)
	  os << indent << "];" << std::endl;
        else
          os << indent << std::endl;
      }
      os << indent << std::endl;
    }

    os << indent << "classifSigma = [";

    for (unsigned int j = 0;j < m_NDim;++j)
    {
      for (unsigned int k = 0;k < m_NDim;++k)
	os << indent << m_ReferenceCovariance(j,k) << " ";

	if (j == m_NDim - 1)
	  os << indent << "];" << std::endl;
        else
          os << indent << std::endl;
    }
  }

  void ContinuousSTAPLEImageFilter::PrintPerformanceParameters(std::string &matlabFileName)
  {
    itk::Indent tmp;
    if (strcmp(matlabFileName.c_str(),"") == 0)
    {
      this->PrintSelf(std::cout,tmp);
    }
    else
    {
      std::ofstream fileOut(matlabFileName.c_str());
      this->PrintSelf(fileOut,tmp);
      fileOut.close();
    }
  }

  bool ContinuousSTAPLEImageFilter::endConditionReached()
  {
    double absDiff = 0;

    for (unsigned int i = 0;i < this->GetNumberOfInputs();++i)
    {
      for (unsigned int k = 0;k < m_NDim;++k)
      {
        if (fabs(m_OldBias[i][k] - m_Bias[i][k]) > absDiff)
	  absDiff = fabs(m_OldBias[i][k] - m_Bias[i][k]);
	  
	for (unsigned int l = k;l < m_NDim;++l)
	{
	  if (fabs(m_OldCovariance[i](k,l) - m_Covariance[i](k,l)) > absDiff)
	    absDiff = fabs(m_OldCovariance[i](k,l) - m_Covariance[i](k,l));
	}
      }
    }

    if (absDiff > m_RelativeConvergenceThreshold)
      return false;
    else
      return true;
  }

} // end namespace itk
