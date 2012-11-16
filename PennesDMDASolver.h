#ifndef __realtimeimagingKalman_h
#define __realtimeimagingKalman_h
/* ------------------------------------------------------------------------
   Main data structure for Kalman filter based off of the following matlab code
 
   disp('Prediction for state vector and covariance');
   s.x = s.A*s.x + s.B*s.u;
   s.P = s.A * s.P * s.A' + s.Q;

   disp('Compute factorization needed for Kalman gain factor');
   %K = s.P*s.H'*inv(s.H*s.P*s.H'+s.R);
   Xtmp = s.P* ((s.P+s.R)\[(s.z-s.x) s.P]);

   disp('Correction based on observation');
   %s.x = s.x + K*(s.z-s.H*s.x);
   %s.P = s.P - K*s.H*s.P;
   s.x = s.x + Xtmp(:,1);
   s.P = s.P - Xtmp( : , 2:size(Xtmp,2) );
  ------------------------------------------------------------------------- */
class KalmanFilter
{

public:
  KalmanFilter(); // constructor

  void DebugOn(); // enable debugging

  // echo data
  virtual void printSelf(std::ostream& os=std::cout) ; 


  //setup structured grid infrastructure
  PetscErrorCode SetupDA( InputImageType::PointType   &,
                          InputImageType::SpacingType &,
                          InputImageType::RegionType::SizeType &);

  //extract image data into petsc Vec
  PetscErrorCode ExtractImageData(InputImageType::Pointer, 
                                  InputImageType::RegionType &,Vec);
  //extract image data into petsc Vec
  PetscErrorCode ExtractCovarianceData(const int, 
                                       InputImageType::Pointer, 
                                       InputImageType::Pointer, 
                                       InputImageType::RegionType &);

  /* Evaluate initial conditions for state and covariance */
  virtual PetscErrorCode InitializeState(Vec) = 0;

  // prediction of state 
  PetscErrorCode StatePredict(Vec,int);

  // apply dirichlet data 
  PetscErrorCode ApplyDirichletData(Vec,PetscScalar);

  // prediction of covariance
  virtual PetscErrorCode CovariancePredict();

  // update of state prediction from measurement
  PetscErrorCode StateUpdate(Vec,Vec);

  // update of covariance from measurement
  virtual PetscErrorCode CovarianceUpdate() = 0;

  //close structured grid infrastructure
  PetscErrorCode FinalizeDA();  

  // write a paraview file from the DA data structures
  PetscErrorCode WriteDAImage(const int, const char *,Vec);

  // processor bounding box
  PetscInt ProcStart[3], ProcStartGhost[3];
  PetscInt ProcWidth[3], ProcWidthGhost[3];

  // structured grid infrastructure
  DA        dac;
  Vec globalVec; 

  // control over direct solver
  MatFactorType factortype;
  char solvertype[PETSC_MAX_PATH_LEN],
       ordertype[PETSC_MAX_PATH_LEN];

protected:
  Mat           CovP,       //covariance of the state vector estimate
                CovR,       //measurement noise covariance 
                CovQ,       //model error covariance 
                SumCov,     //data structure used to hold P + R and factor
                tmpMatDen,  //dense matrix used as temporary storage
                StateX,     // matrix to hold solution for covariance update
                StateXtrans;//    StateA StateX = StateB 
  PetscTruth  noDenseUpdate;
  PetscScalar modelcov,  // initial model covariance
              statecov;  // initial state covariance
  PetscScalar zerotol; // pointer for raw memory allocation

  PetscScalar bodyTemp; // IC for raw memory allocation

  // member functions to initialize covariance
  PetscErrorCode InitializeSparseUncorrCov(Mat&,PetscScalar);
  PetscErrorCode InitializeDenseUncorrCov( Mat&,PetscScalar);

  // member functions to convert whatever is store in the dense temporary
  // storage to the sparse format
  PetscErrorCode GetSparseMatrixFromTmp(Mat &);

  MatFactorInfo factorInfo;

private:

  // data structures for plotting
  Vec          imageVec, loadVec,tmpVec,covMaxVec;
  VecScatter     gather;

  // parallel solver
  KSP           kspPre;  /* linear solver for state prediction*/
  Mat           StateA,     // matrix operators for state prediction 
                StateB;     // A u_n+1 =  f - B  u_n 
  PetscInt maxiter; // max ksp iterations
 
  PetscScalar      mu_tr,mu_a,mu_s,anfact,mu_eff,gflux,
                   u_a,w_0,c_blood,
                   c_p,rho,
                   k_0,
                   deltat,volumeFraction;

  std::vector< PetscScalar > X_0,Y_0,Z_0,Power;
  std::vector< PetscInt >    dirichletNodes;

  InputImageType::PointType   ROIorgn;
  InputImageType::RegionType::SizeType ROIsize;

};
/*-------------------------------------------------------------------*/
class SparseKalmanFilter : public KalmanFilter
{

public:
  SparseKalmanFilter(); // constructor

  /* Evaluate initial conditions for state and covariance */
  virtual PetscErrorCode InitializeState(Vec);

  // update of covariance from measurement
  virtual PetscErrorCode CovarianceUpdate(); 

private:
  // none
};
/*-------------------------------------------------------------------*/
class UncorrelatedKalmanFilter : public SparseKalmanFilter
{

public:
  UncorrelatedKalmanFilter();// constructor

  // update of covariance from measurement
  virtual PetscErrorCode CovariancePredict(); 

  // update of covariance from measurement
  virtual PetscErrorCode CovarianceUpdate(); 

private:
  // none
};
/*-------------------------------------------------------------------*/
class DenseKalmanFilter : public KalmanFilter
{

public:
  DenseKalmanFilter(); // constructor

  /* Evaluate initial conditions for state and covariance */
  virtual PetscErrorCode InitializeState(Vec);

  // update of covariance from measurement
  virtual PetscErrorCode CovarianceUpdate(); 

private:
  // none
};
#endif
