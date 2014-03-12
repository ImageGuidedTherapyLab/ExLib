// system includes
#include <math.h>
#include <string>
#include <vector>

// libmesh includes
#include "libmesh.h"
#include "getpot.h"

// local includes
#include "baseInfo.h"
#include "Setup.h"
#include "Imaging.h"
#include "KalmanFilter.h"
#include "itkVTKImageVariableNameIO.h"

// global initial tolerance
const PetscScalar rtol0 = 1.e-9;

/* ------------------------------------------------------------------------

    FD pennes model
    the partial differential equation
  
            \rho dudt - k Laplacian u = q_laser
  
    with dirichlet boundary conditions
   
             u = 0  for  x = x_min, x = x_max, 
                         y = y_min, y = y_max, 
                         z = z_min, z = z_max
  
    A finite difference approximation with the usual 7-point stencil
    is used to discretize the boundary value problem to obtain a linear 
    system of equations.


  ------------------------------------------------------------------------- */

/* constructor */
KalmanFilter::KalmanFilter( )
{
  rho      = 1000.0;   // [kg/m^3]
  c_p      = 3600.0;   // [j / kg / k]
  k_0      = 0.57;     // [j /s / m /k] 
  w_0      = 6.00;     // [kg /s / m^3] 
  u_a      = 37.0;     // [K]
  c_blood  = 4180.0;   // [j / kg / k]
  mu_a     = 5.0e+2;   // [m^-1]
  mu_s     = 140.0e+2; // [m^-1]
  anfact   = 0.862;   // dimensionless
  deltat   = 5.0;     // seconds
  /*
    Gaussian of the form
     f(x) = 1/\sqrt{2 \pi P} \exp(-1/2/P (x -\mu)^2 ) 

    sigma^2 = P

    (\mu- \sigma,\mu+ \sigma)  ---> 0.683 confidence
    (\mu-2\sigma,\mu+2\sigma)  ---> 0.997 confidence
 
  */
  modelcov = 2.5*2.5;     // [deg C]
  statecov = 1.5    ;     // [deg C]
  gflux    = 0.0;     // [watt / m^2]
  SumCov =NULL; tmpMatDen=NULL;StateX=NULL; 
  zerotol=1.e-6;StateXtrans=NULL; 

  PetscTruth  flg;
  noDenseUpdate = PETSC_FALSE;
  PetscErrorCode ierr = PetscOptionsGetTruth(PETSC_NULL,"-noupdate",&noDenseUpdate,&flg);
}

void KalmanFilter::printSelf(std::ostream& os)
{
 PetscErrorCode ierr;
 ierr=PetscPrintf(PETSC_COMM_WORLD,
      "rho=%f c_p=%f k_0=%f \n",
                     rho,c_p,k_0); CHKERRV(ierr);

 ierr=PetscPrintf(PETSC_COMM_WORLD,
      "c_blood=%f u_a=%f w_0=%f mu_a=%f mu_s=%f\n",
                     c_blood,u_a,w_0,mu_a,mu_s); CHKERRV(ierr);

 ierr=PetscPrintf(PETSC_COMM_WORLD,
    "anfact=%f deltat=%f volumeFrac=%f\n", 
     anfact,deltat,volumeFraction);CHKERRV(ierr);
 printStdVector(      os, "          X_0[",X_0);
 printStdVector(      os, "          Y_0[",Y_0);
 printStdVector(      os, "          Z_0[",Z_0);

 return;
}
SparseKalmanFilter::SparseKalmanFilter( ):KalmanFilter( )
{
}

DenseKalmanFilter::DenseKalmanFilter( ):KalmanFilter( )
{
}

UncorrelatedKalmanFilter::UncorrelatedKalmanFilter():SparseKalmanFilter() // constructor
{
  noDenseUpdate = PETSC_TRUE;
}

//setup Petsc DA data structures
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::SetupDA"
PetscErrorCode 
KalmanFilter::SetupDA( InputImageType::PointType  &orgn,
                       InputImageType::SpacingType  &sp,
                       InputImageType::RegionType::SizeType &size)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // command line parameters
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-rho"    ,&rho    ,PETSC_NULL);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-c_p"    ,&c_p    ,PETSC_NULL);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-k_0"    ,&k_0    ,PETSC_NULL);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-w_0"    ,&w_0    ,PETSC_NULL);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-u_a"    ,&u_a    ,PETSC_NULL);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-c_blood",&c_blood,PETSC_NULL);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-mu_a"   ,&mu_a   ,PETSC_NULL);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-mu_s"   ,&mu_s   ,PETSC_NULL);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-anfact" ,&anfact ,PETSC_NULL);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-statecov",&statecov,PETSC_NULL);
  ierr=PetscOptionsGetScalar(PETSC_NULL,"-modelcov",&modelcov,PETSC_NULL);
  deltat  = static_cast<PetscScalar>(Setup::acquisitionTime);

  // open using ITK library of file types
  InputReaderType::Pointer  reader = InputReaderType::New();
  reader->SetFileName( "files/probeData.mha" );
  try
   {
     reader->Update();

     // setup iterator
     InputIteratorType   tempIt( reader->GetOutput()   ,
                                 reader->GetOutput()->GetRequestedRegion() );
     tempIt.SetFirstDirection(  0 );  tempIt.SetSecondDirection( 1 );
     tempIt.GoToBegin();

     std::cout << "   Setting up WFS Model" << std::endl << std::endl ;

     // reset data structures
     X_0.clear(); Y_0.clear(); Z_0.clear(); 
     PetscScalar pixelTotal = 0.0;

     /* loop through parallel data structures */
     for (PetscInt k=0; k<Setup::size[2]; k++) 
     {
       for (PetscInt j=0; j<Setup::size[1]; j++) 
       {
         for (PetscInt i=0; i<Setup::size[0]; i++) 
         {
           if(tempIt.Get() >= 1.0)
             {
             if(tempIt.Get() >= 2.0) // diffusing fiber 
               {
                X_0.push_back(Setup::orgn[0] + sp[0] * (i + 0.5) );
                Y_0.push_back(Setup::orgn[1] + sp[1] * (j + 0.5) );
                Z_0.push_back(Setup::orgn[2] + sp[2] *    k      );
                // update pixel count
                pixelTotal = pixelTotal + 1.0;
               }

              // store dirichlet data
              PetscInt roiIndex[3] = { i - Setup::index[0],
                                       j - Setup::index[1],
                                       k - Setup::index[2]};
              if( roiIndex[0] &&  roiIndex[0] < size[0] &&
                  roiIndex[1] &&  roiIndex[1] < size[1] &&
                  roiIndex[2] &&  roiIndex[2] < size[2]  )
                {
                  dirichletNodes.push_back( roiIndex[2] * size[1] * size[0]
                                           +roiIndex[1] * size[0]
                                           +roiIndex[0] );  
                }
             }
           ++tempIt; // update iterators
         }
         // get next line
         tempIt.NextLine(); 
       }
       // get next slice
       tempIt.NextSlice(); 
     }
     CHKMEMQ; // check for memory corruption use -malloc_debug to enable

     // compute the volume fraction
     volumeFraction = 1.0/pixelTotal;
   }
  catch (itk::ExceptionObject &excp)
    { // default is No "probeDomain" about the middle
      PetscInt nelemtip = 1;
      ierr=PetscOptionsGetInt(PETSC_NULL , "-nelemtip" , &nelemtip, PETSC_NULL);
      X_0.resize(nelemtip, orgn[0]+sp[0]*size[0]/2 + 0.5 * sp[0] );     // [m]
      Y_0.resize(nelemtip, orgn[1]+sp[1]*size[1]/2 + 0.5 * sp[1] );     // [m]
      Z_0.resize(nelemtip, orgn[2]+sp[2]*size[2]/2               );     // [m]
      ierr=PetscOptionsGetScalar(PETSC_NULL , "-X_0" , &X_0[0] , PETSC_NULL);
      ierr=PetscOptionsGetScalar(PETSC_NULL , "-Y_0" , &Y_0[0] , PETSC_NULL);
      ierr=PetscOptionsGetScalar(PETSC_NULL , "-Z_0" , &Z_0[0] , PETSC_NULL);
      if( nelemtip == 1 )
       {
        std::cout << "   Setting up SDA Model" << std::endl << std::endl ;
        volumeFraction = 1.0;
       }
      else if( nelemtip > 1 )
       {
        std::cout << "   Setting up line source WFS " << std::endl << std::endl;
        PetscScalar X_1 = X_0[0]    , // default to y-direction
                    Y_1 = Y_0[0]+1.0, // default to y-direction
                    Z_1 = Z_0[0]    , // default to y-direction
                    diffusinglength = 0.01,   // default to 1.0cm
                    diffusingradius = 0.00075; // default to diameter 1.5mm
        ierr=PetscOptionsGetScalar(PETSC_NULL , "-X_1" , &X_1 , PETSC_NULL);
        ierr=PetscOptionsGetScalar(PETSC_NULL , "-Y_1" , &Y_1 , PETSC_NULL);
        ierr=PetscOptionsGetScalar(PETSC_NULL , "-Z_1" , &Z_1 , PETSC_NULL);
        ierr=PetscOptionsGetScalar(PETSC_NULL , "-diffusinglength" , 
                                                 &diffusinglength , PETSC_NULL);
        ierr=PetscOptionsGetScalar(PETSC_NULL , "-diffusingradius" , 
                                                 &diffusingradius , PETSC_NULL);
        PetscScalar elemLength = diffusinglength / nelemtip ;
        // compute direction vector
        // (X_1,Y_1,Z_1) - (X_0,Y_0,Z_0) 
        //  should point towards laser tip as the reference
        PetscScalar differenceNorm = std::sqrt(
                                           std::pow(X_1 - X_0[0],2)   +
                                           std::pow(Y_1 - Y_0[0],2)   +
                                           std::pow(Z_1 - Z_0[0],2)   
                                              );
        PetscScalar unitVec[3] = 
                             { (X_1 - X_0[0]) / differenceNorm,
                               (Y_1 - Y_0[0]) / differenceNorm,
                               (Z_1 - Z_0[0]) / differenceNorm};
        // compute laser starting point 
        // (xhat_dir,yhat_dir,zhat_dir) should point away from laser tip
        PetscScalar laserTip[3] = 
                           { X_0[0] - unitVec[0]*0.5*diffusinglength , 
                             Y_0[0] - unitVec[1]*0.5*diffusinglength , 
                             Z_0[0] - unitVec[2]*0.5*diffusinglength}; 

        volumeFraction = elemLength / diffusinglength ;
        /**
         *   put points at centroids of elements along the tip
         *    |--x--|--x--|--x--|--x--|--x--|
         */
        for ( PetscInt Ii = 0 ; Ii < nelemtip ; Ii++  )
          {
            X_0[Ii] = laserTip[0] + unitVec[0]* (Ii+0.5) * elemLength;
            Y_0[Ii] = laserTip[1] + unitVec[1]* (Ii+0.5) * elemLength;
            Z_0[Ii] = laserTip[2] + unitVec[2]* (Ii+0.5) * elemLength;
          }
       for (PetscInt k=0; k<Setup::size[2]; k++) 
        for (PetscInt j=0; j<Setup::size[1]; j++) 
         for (PetscInt i=0; i<Setup::size[0]; i++) 
          {
           PetscScalar locpos[3]=
                             {Setup::orgn[0] + (i+0.5)*Setup::sp[0], 
                              Setup::orgn[1] + (j+0.5)*Setup::sp[1], 
                              Setup::orgn[2] + (k+0.5)*Setup::sp[2]};
           PetscScalar locposRefTip[3]=
                             {locpos[0] - laserTip[0],
                              locpos[1] - laserTip[1],
                              locpos[2] - laserTip[2]};
           PetscScalar axialComponent =
                            locposRefTip[0]* unitVec[0] +
                            locposRefTip[1]* unitVec[1] +
                            locposRefTip[2]* unitVec[2] ;
           if( axialComponent > 0.0)
             { // point  is along the fiber tract
               PetscScalar radialVec[3]=
                   {locposRefTip[0] - axialComponent * unitVec[0],
                    locposRefTip[1] - axialComponent * unitVec[1],
                    locposRefTip[2] - axialComponent * unitVec[2]};
               // check if satisfy the radius criteria
               PetscScalar locRadius = std::sqrt(
                                       std::pow(radialVec[0],2)   +
                                       std::pow(radialVec[1],2)   +
                                       std::pow(radialVec[2],2)   
                                                );
               if( locRadius < diffusingradius )
                  dirichletNodes.push_back( k * size[1] * size[0]
                                           +j * size[0]
                                           +i );  
             }
          }
       }
      else 
       {
        std::cout << "   error input laser model..." << std::endl << std::endl ;
        abort();
       }

    } 

  ierr=PetscOptionsGetScalar(PETSC_NULL,"-zerotol" ,&zerotol,PETSC_NULL);
  // get the power data
  std::string powerFile("files/power.ini");
  GetPowerData(powerFile, Power);

  // default solver to superlu_dist
  sprintf(solvertype,"%s",MAT_SOLVER_PLAPACK);
  sprintf(solvertype,"%s",MAT_SOLVER_MUMPS);
  sprintf(solvertype,"%s",MAT_SOLVER_SUPERLU_DIST);
  PetscTruth  flg;
  ierr = PetscOptionsGetString(PETSC_NULL,"-solver",solvertype,
                                           PETSC_MAX_PATH_LEN,&flg);

  // default factorization to LU
  factortype = MAT_FACTOR_CHOLESKY;
  factortype = MAT_FACTOR_LU;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-factor",(PetscInt *)&factortype,&flg);

  // default ordering to nested dissection
  sprintf(ordertype,"%s",MATORDERING_RCM);
  sprintf(ordertype,"%s",MATORDERING_ND);
  ierr = PetscOptionsGetString(PETSC_NULL,"-ordering",ordertype,
                                           PETSC_MAX_PATH_LEN,&flg);

  // params for matrix factorization
  ierr = MatFactorInfoInitialize(&factorInfo); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DA) to manage parallel grid and vectors
       use -da_view to print out info about the DA 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DAPeriodicType ptype = DA_NONPERIODIC;
  DAStencilType  stype = DA_STENCIL_STAR;
  ROIsize = size ; 

  // extend out of plane dimension 
  // TODO fix bugs caused by this
  PetscInt extendPlane = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-extendplane",&extendPlane,PETSC_NULL);
  if(extendPlane) ROIsize[2] = ROIsize[2] + extendPlane; 

  ierr = DACreate3d(PETSC_COMM_WORLD,ptype,stype, ROIsize[0], ROIsize[1], ROIsize[2],
                    PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,nvarplot,1,
                    PETSC_NULL,PETSC_NULL,PETSC_NULL,&dac);CHKERRQ(ierr);
  ierr = DAGetCorners(dac,&ProcStart[0],&ProcStart[1],&ProcStart[2], 
                          &ProcWidth[0],&ProcWidth[1],&ProcWidth[2]); 
  CHKERRQ(ierr);

  // get dimensions w/ ghost cell info
  ierr = DAGetGhostCorners(dac, 
          &ProcStartGhost[0],&ProcStartGhost[1],&ProcStartGhost[2], 
          &ProcWidthGhost[0],&ProcWidthGhost[1],&ProcWidthGhost[2] ); 
  CHKERRQ(ierr);
  
  ROIorgn = orgn ; 
  ierr = DASetUniformCoordinates(dac, ROIorgn[0], ROIorgn[0]+sp[0]*ROIsize[0] ,
                                      ROIorgn[1], ROIorgn[1]+sp[1]*ROIsize[1] ,
                                      ROIorgn[2], ROIorgn[2]+sp[2]*ROIsize[2]  );
  CHKERRQ(ierr);

  // Extract global vectors from DA; 
  DACreateGlobalVector(dac,&globalVec);
  VecDuplicate(globalVec,    &loadVec); 
  VecDuplicate(globalVec,     &tmpVec); 
  VecDuplicate(globalVec,  &covMaxVec); 
  VecScatterCreateToZero(globalVec,&gather,&imageVec);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create matrix data structure; 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DAGetMatrix(dac,MATAIJ,   &StateA  );CHKERRQ(ierr);
  ierr = DAGetMatrix(dac,MATAIJ,   &StateB  );CHKERRQ(ierr);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Build State Transition Matrix
     Compute entries for the locally owned part of the Jacobian.
      - Currently, all PETSc parallel matrix formats are partitioned by
        contiguous chunks of rows across the processors. 
      - Each processor needs to insert only elements that it owns
        locally (but any non-local elements will be sent to the
        appropriate processor during matrix assembly). 
      - Here, we set all entries for a particular row at once.
      - We can set matrix entries either using either
        MatSetValuesLocal() or MatSetValues()
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscScalar    v[7],two = 2.0,hx,hy,hz;
  MatStencil     col[7],row;

  // spacing parameters
  hx     = Setup::sp[0]; // 1.0/(PetscReal)(size[0]-1);
  hy     = Setup::sp[1]; // 1.0/(PetscReal)(size[1]-1);
  hz     = Setup::sp[2]; // 1.0/(PetscReal)(size[2]-1);

  // fd parameters
  PetscScalar sc      = hx*hz*hy;
  PetscScalar hxhzdhy = hx*hz/hy;
  PetscScalar hyhzdhx = hy*hz/hx;
  PetscScalar hxhydhz = hx*hy/hz;
  ierr=PetscPrintf(PETSC_COMM_WORLD,"sc=%12.5e hxhzdhy=%12.5e hyhzdhx=%12.5e hxhydhz=%12.5e \n",
                                     sc,hxhzdhy,hyhzdhx,hxhydhz);CHKERRQ(ierr);

  // build matrix entries
  //  take FD equations and multiply by hx * hy * hz before assembling
  if(ROIsize[2]==1){ // 2D
     for (int j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) {
       std::cout << " " << j << std::flush;
       for (int i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) {
         row.k = 0; row.j = j; row.i = i;
         /* boundary points or 2D */
         if (i==0 || j==0 || i==ROIsize[0]-1 || j==ROIsize[1]-1 ) {
           v[0] = 1.0;
           ierr = MatSetValuesStencil(StateA,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
           ierr = MatSetValuesStencil(StateB,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
         } else {
         /* interior grid points */
           v[0] = -0.5 * k_0 * hxhzdhy; col[0].k=0;  col[0].j=j-1;col[0].i = i;
           v[1] = -0.5 * k_0 * hyhzdhx; col[1].k=0;  col[1].j=j;  col[1].i = i-1;
           v[2] =  sc*( rho*c_p/deltat + 0.5 * w_0 * c_blood) 
                      + 1.0 * k_0 * (hyhzdhx+hxhzdhy);
                            col[2].k=row.k;col[2].j=row.j;col[2].i = row.i;
           v[3] = -0.5 * k_0 * hyhzdhx; col[3].k=0;  col[3].j=j;  col[3].i = i+1;
           v[4] = -0.5 * k_0 * hxhzdhy; col[4].k=0;  col[4].j=j+1;col[4].i = i;
           ierr = MatSetValuesStencil(StateA,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
           v[2] =  sc*(-rho*c_p/deltat + 0.5 * w_0 * c_blood)
                      + 1.0 * k_0 * (hyhzdhx+hxhzdhy);
           ierr = MatSetValuesStencil(StateB,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
         }
       } // end for (int i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) 
     } // end for (int j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) 
   }else{ // 3D
     for (int k=ProcStart[2]; k<ProcStart[2]+ProcWidth[2]; k++) {
       for (int j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) {
         for (int i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) {
           row.k = k; row.j = j; row.i = i;
           /* boundary points or 2D */
           if (i==0 || j==0 || i==ROIsize[0]-1 || j==ROIsize[1]-1 ) {
             v[0] = 1.0;
             ierr = MatSetValuesStencil(StateA,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
             ierr = MatSetValuesStencil(StateB,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
           } else if (k ==       0  ) {// out of plane zero flux bc
             v[0] = -0.5 * k_0 / hz; col[0].k=k+1;col[0].j=j;  col[0].i = i;
             v[1] =  0.5 * k_0 / hz; col[1].k=k  ;col[1].j=j;  col[1].i = i;
             ierr = MatSetValuesStencil(StateA,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
             ierr = MatSetValuesStencil(StateB,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
           } else if (k == ROIsize[2]-1) {// out of plane zero flux bc
             v[0] = -0.5 * k_0 / hz; col[0].k=k  ;col[0].j=j;  col[0].i = i;
             v[1] =  0.5 * k_0 / hz; col[1].k=k-1;col[1].j=j;  col[1].i = i;
             ierr = MatSetValuesStencil(StateA,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
             ierr = MatSetValuesStencil(StateB,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
           } else {
           /* interior grid points */
             v[0] = -0.5 * k_0 * hxhydhz; col[0].k=k-1;col[0].j=j;  col[0].i = i;
             v[1] = -0.5 * k_0 * hxhzdhy; col[1].k=k;  col[1].j=j-1;col[1].i = i;
             v[2] = -0.5 * k_0 * hyhzdhx; col[2].k=k;  col[2].j=j;  col[2].i = i-1;
             v[3] =  sc*( rho*c_p/deltat + 0.5 * w_0 * c_blood) 
                       + 1.0 * k_0 * (hyhzdhx+hxhzdhy+hxhydhz);
                              col[3].k=row.k;col[3].j=row.j;col[3].i = row.i;
             v[4] = -0.5 * k_0 * hyhzdhx; col[4].k=k;  col[4].j=j;  col[4].i = i+1;
             v[5] = -0.5 * k_0 * hxhzdhy; col[5].k=k;  col[5].j=j+1;col[5].i = i;
             v[6] = -0.5 * k_0 * hxhydhz; col[6].k=k+1;col[6].j=j;  col[6].i = i;
             ierr = MatSetValuesStencil(StateA,1,&row,7,col,v,INSERT_VALUES);CHKERRQ(ierr);
             v[3] =  sc*(-rho*c_p/deltat + 0.5 * w_0 * c_blood) 
                       + 1.0 * k_0 * (hyhzdhx+hxhzdhy+hxhydhz);
             ierr = MatSetValuesStencil(StateB,1,&row,7,col,v,INSERT_VALUES);CHKERRQ(ierr);
           }
         } // end for (int i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) 
       } // end for (int j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) 
     } // end for (int k=ProcStart[2]; k<ProcStart[2]+ProcWidth[2]; k++) 
   }

  /* Assemble matrix, using the 2-step process:
       MatAssemblyBegin(), MatAssemblyEnd().  */
  ierr = MatAssemblyBegin(StateA,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  StateA,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(StateB,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  StateB,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatScale(StateB,-1.0);CHKERRQ(ierr);

  if(dirichletNodes.size())
   {
    ierr = MatZeroRows(StateA,dirichletNodes.size(),&dirichletNodes[0],1.0);
    ierr = MatZeroRows(StateB,dirichletNodes.size(),&dirichletNodes[0],1.0);
   }

  ierr = MatDataInfo(StateA,"(SetupDA) StateA");CHKERRQ(ierr);
  ierr = MatDataInfo(StateB,"(SetupDA) StateB");CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create linear solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = KSPCreate(PETSC_COMM_WORLD,&kspPre);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(kspPre,PETSC_TRUE);

  ierr = KSPSetOperators(kspPre,StateA,StateA,DIFFERENT_NONZERO_PATTERN);
  CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize linear solver; set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PC pcPre;
  ierr = KSPGetPC(kspPre,&pcPre); CHKERRQ(ierr);
  ierr = PCSetType(pcPre,PCBJACOBI); CHKERRQ(ierr);
  maxiter= ROIsize[0] * ROIsize[1] * ROIsize[2] ;
  ierr = KSPSetTolerances(kspPre,rtol0,PETSC_DEFAULT,PETSC_DEFAULT, maxiter); 
  CHKERRQ(ierr);
  ierr = KSPSetFromOptions(kspPre);CHKERRQ(ierr);

  // no covariance prediction
  if(noDenseUpdate) 
   { 
    PetscScalar dummy = 1.0; // dummy initialize for static gaussian model
    ierr = InitializeSparseUncorrCov(SumCov, dummy);CHKERRQ(ierr);
    PetscFunctionReturn(0);
   }

  // scratch matrix to hold solution for covariance update
  PetscInt M,N,m,n;
  ierr = MatGetSize(StateA,&M,&N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(StateA,&m,&n);CHKERRQ(ierr);
  ierr = MatCreateMPIDense(PETSC_COMM_WORLD,m,n,M,N,
                            PETSC_NULL,&tmpMatDen);CHKERRQ(ierr);
  // factor the matrix StateAFact
  Mat StateAFact;
  ierr = PetscPrintf(PETSC_COMM_WORLD,
               "Factor Matrix SOLVER %s FACTOR %d ORDER %s size (%dx%d)...\n",
                             solvertype,factortype,ordertype,M,N);CHKERRQ(ierr);

  IS isrow,iscol;
  if(!strcasecmp(solvertype,MAT_SOLVER_PLAPACK)) { // plapack uses a dense format
      Mat StateAtmp;
      ierr = MatConvert(StateA,MATMPIDENSE,MAT_INITIAL_MATRIX,&StateAtmp);
      ierr = MatGetFactor(StateAtmp,solvertype,
                             factortype,&StateAFact); CHKERRQ(ierr);
      ierr = MatGetOrdering(StateAtmp,ordertype,&isrow,&iscol); CHKERRQ(ierr);
      ierr = MatLUFactorSymbolic(StateAFact,StateAtmp,isrow,iscol,&factorInfo);CHKERRQ(ierr);
      ierr = MatLUFactorNumeric(StateAFact,StateAtmp,&factorInfo);CHKERRQ(ierr);
      ierr = MatDestroy(StateAtmp);CHKERRQ(ierr);
  } else { // default is mpiiaj 
      ierr = MatGetFactor(StateA,solvertype,
                             factortype,&StateAFact); CHKERRQ(ierr);
      ierr = MatGetOrdering(StateA,ordertype,&isrow,&iscol); CHKERRQ(ierr);
      ierr = MatLUFactorSymbolic(StateAFact,StateA,isrow,iscol,&factorInfo);CHKERRQ(ierr);
      ierr = MatLUFactorNumeric(StateAFact,StateA,&factorInfo);CHKERRQ(ierr);
  }

  Mat StateBtmp;  // covariance of the state vector estimate
  ierr = MatConvert(StateB,MATMPIDENSE,MAT_INITIAL_MATRIX,&StateBtmp);
  // store tmpMatDen = A^-1 B
  ierr = MatMatSolve(StateAFact,StateBtmp,tmpMatDen);CHKERRQ(ierr);

  // assemble if necessary 
  PetscTruth  assembled ; 
  ierr = MatAssembled(tmpMatDen,&assembled); CHKERRQ(ierr);
  if(!assembled)
    { 
     ierr = MatAssemblyBegin(tmpMatDen,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
     ierr = MatAssemblyEnd(tmpMatDen,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }
 
  // free memory
  ierr = MatDestroy(StateBtmp);CHKERRQ(ierr);
  ierr = MatDestroy(StateAFact);CHKERRQ(ierr);

  // setup data structure for covariance update
  ierr=PetscPrintf(PETSC_COMM_WORLD,
            "Setup Work Data Struc for Covariance Update...\n");CHKERRQ(ierr);

  if( !strcasecmp(solvertype,MAT_SOLVER_PLAPACK)) { //plapack is dense
    ierr = MatCreateMPIDense(PETSC_COMM_WORLD,m,n,M,N,
                              PETSC_NULL,&SumCov);CHKERRQ(ierr);
  } else { // default is mpiiaj 
    ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,m,n,M,N,n,PETSC_NULL,
                                        N-n,PETSC_NULL,&SumCov);CHKERRQ(ierr);
    PetscInt Istart,Iend;
    ierr = MatGetOwnershipRange(SumCov,&Istart,&Iend);CHKERRQ(ierr);
    /* Insert Ones and assemble to ensure full matrix is allocated */
    for (PetscInt Ii=Istart; Ii<Iend; Ii++)
       for (PetscInt Jj=0; Jj<N; Jj++)
          MatSetValue(SumCov,Ii,Jj,1.0,INSERT_VALUES);
    ierr = MatAssemblyBegin(SumCov,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(SumCov,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  }

  // Software Guide : EndCodeSnippet
  PetscFunctionReturn(0);
}

// Free Petsc Data structures
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::FinalizeDA"
PetscErrorCode KalmanFilter::FinalizeDA()
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  // destroy scatter context, vector, when no longer needed
  if(StateX      !=NULL) {ierr = MatDestroy(StateX);CHKERRQ(ierr); }
  if(StateXtrans !=NULL) {ierr = MatDestroy(StateXtrans);CHKERRQ(ierr); }
  if(tmpMatDen   !=NULL) {ierr = MatDestroy(tmpMatDen);CHKERRQ(ierr); }
  if(SumCov      !=NULL) {ierr = MatDestroy(SumCov    );CHKERRQ(ierr); }
  ierr = VecScatterDestroy(gather);CHKERRQ(ierr);
  ierr = VecDestroy(loadVec);CHKERRQ(ierr);
  ierr = VecDestroy(tmpVec);CHKERRQ(ierr);
  ierr = VecDestroy(covMaxVec);CHKERRQ(ierr);
  ierr = VecDestroy(imageVec);CHKERRQ(ierr);
  ierr = KSPDestroy(kspPre);CHKERRQ(ierr);
  //ierr = VecDestroy(globalVec);CHKERRQ(ierr); DA will destroy this
  ierr = DADestroy(dac);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// Prediction of state 
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::StatePredict"
PetscErrorCode KalmanFilter::StatePredict(Vec State, int istep)
{

  PetscErrorCode ierr;
  PetscFunctionBegin;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    s.x = s.A*s.x + s.B*s.u;
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  // spacing parameters
  PetscScalar    hx,hy,hz,sc;
  hx     = Setup::sp[0]; // 1.0/(PetscReal)(size[0]-1);
  hy     = Setup::sp[1]; // 1.0/(PetscReal)(size[1]-1);
  hz     = Setup::sp[2]; // 1.0/(PetscReal)(size[2]-1);
  sc      = hx*hz*hy;
  //SOURCE TERM
  PetscScalar dist=0.0;
  mu_tr=mu_a+mu_s*(1.0-anfact);
  mu_eff=sqrt(3.0*mu_a*mu_tr);

  /* Get pointers to vector data */
  PetscScalar    ***source;
  ierr = DAVecGetArray(dac,globalVec,&source);CHKERRQ(ierr);

  /* Compute function over the locally owned part of the grid */
  for (int k=ProcStart[2]; k<ProcStart[2]+ProcWidth[2]; k++) {
    for (int j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) {
      for (int i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) {
        source[k][j][i] = 0.0;
        // out of plane flux bc
        if (k ==  0 || k == ROIsize[2]-1 ) 
           source[k][j][i] = source[k][j][i] + gflux; 
        // laser model
        for (unsigned int Ii = 0; Ii < X_0.size(); Ii++ )
         {
          dist=std::sqrt( std::pow(ROIorgn[0] + (i+0.5) * Setup::sp[0]  - X_0[Ii],2) 
                        + std::pow(ROIorgn[1] + (j+0.5) * Setup::sp[1]  - Y_0[Ii],2) 
                        + std::pow(ROIorgn[2] +    k    * Setup::sp[2]  - Z_0[Ii],2) );
          /* libmesh_assert (dist != 0.0);  NOT needed. check for nan in norm */
          source[k][j][i] = source[k][j][i] + 
                             sc * 0.75 * volumeFraction * Power[istep] 
                                * mu_a * mu_tr
                                * std::exp(-mu_eff*dist) / libMesh::pi / dist;
         }

        //if(source[k][j][i] != 0.0) 
        //       std::cout<< source[k][j][i] << " i= " << i 
        //                                   << " j= " << j 
        //                                   << " k= " << k
        //                << std::endl << std::flush;

        // perfusion term
        source[k][j][i] = source[k][j][i] + sc * w_0 * c_blood * u_a;

      }
    }
  }

  /* Restore vectors */
  ierr = DAVecRestoreArray(dac,globalVec,&source);CHKERRQ(ierr);

  // compute loadVec = globalVec + B * State
  ierr = MatMultAdd(StateB, State, globalVec, loadVec);CHKERRQ(ierr);

  /* apply dirichlet data */
  PetscScalar    ***load;
  ierr = DAVecGetArray(dac,loadVec,&load);CHKERRQ(ierr);

  for (int k=ProcStart[2]; k<ProcStart[2]+ProcWidth[2]; k++) {
    for (int j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) {
      for (int i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) {
        if (i == 0 || j == 0 || i == ROIsize[0]-1 || j == ROIsize[1]-1 ) 
                                         load[k][j][i] = Setup::bodyTemp; 
      }
    }
  }

  /* Restore vectors */
  ierr = DAVecRestoreArray(dac,loadVec,&load);CHKERRQ(ierr);

  // dirichlet data
  ierr = this->ApplyDirichletData(loadVec,Setup::probeTemp);CHKERRQ(ierr);

  // print info
  ierr = VecDataInfo(globalVec,"(Prediction) globalVec");CHKERRQ(ierr);
  ierr = VecDataInfo(State    ,"(Prediction) State_n-1");CHKERRQ(ierr);
  ierr = VecDataInfo(loadVec  ,"(Prediction) loadVec  ");CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve linear system overwrite current state
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr=PetscPrintf(PETSC_COMM_WORLD,
                 "Prediction for state vector...\n");CHKERRQ(ierr);
  
  // nan check in load vector
  PetscScalar   loadNorm;
  ierr = VecNorm(loadVec, NORM_2, &loadNorm); CHKERRQ(ierr);
  if (loadNorm != loadNorm)
   {
      std::cout << "(Prediction) nan detected in load vector... "
                << "cannot Solve! at time step " << istep
                << std::endl << std::flush ; 
      abort();
   }
  ierr = KSPSolve(kspPre,loadVec,State);CHKERRQ(ierr);
  ierr = VecDataInfo(State, "(Prediction) State");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                  s.P = s.A * s.P * s.A' + s.Q
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::CovariancePredict"
PetscErrorCode KalmanFilter::CovariancePredict()
{

  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr=PetscPrintf(PETSC_COMM_WORLD,
                   "Prediction for covariance...\n");CHKERRQ(ierr);

  if(noDenseUpdate) PetscFunctionReturn(0);

  Mat locMat;

  ierr=PetscPrintf(PETSC_COMM_WORLD,"MatMatMult... \n");CHKERRQ(ierr);
  ierr = MatDataInfo(CovP,"(Prediction) CovP");CHKERRQ(ierr);
  ierr = MatMatMult(StateX,CovP,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&locMat);
  ierr = MatDestroy(CovP);CHKERRQ(ierr);
  ierr = MatMatMult(locMat,StateXtrans,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CovP);
  ierr = MatDataInfo(CovQ,"(Prediction) CovQ");CHKERRQ(ierr);
  ierr = MatAXPY(CovP,1.0,CovQ,SUBSET_NONZERO_PATTERN);
  ierr = MatDataInfo(CovP,"(Prediction) CovP");CHKERRQ(ierr);

  // free memory
  ierr = MatDestroy(locMat);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                  s.P = s.A * s.P * s.A' + s.Q
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#undef __FUNCT__
#define __FUNCT__ "UncorrelatedKalmanFilter::CovariancePredict"
PetscErrorCode UncorrelatedKalmanFilter::CovariancePredict()
{

  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr=PetscPrintf(PETSC_COMM_WORLD,
                   "Prediction for covariance...\n");CHKERRQ(ierr);
  ierr = MatDataInfo(CovP,"(Prediction) CovP");CHKERRQ(ierr);
  ierr = MatDataInfo(CovQ,"(Prediction) CovQ");CHKERRQ(ierr);
  ierr = MatAXPY(CovP,1.0,CovQ,SUBSET_NONZERO_PATTERN);
  ierr = MatDataInfo(CovP,"(Prediction) CovP");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "DenseKalmanFilter::InitializeState"
/* 
   DirectInitializeState - Forms initial state and covariance

   Input Parameters:
   StateTemp - vector

 */
PetscErrorCode DenseKalmanFilter::InitializeState(Vec StateTemp)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  // initial conditions are for Temp difference
  ierr = VecSet(StateTemp,0.0);CHKERRQ(ierr);

  if(noDenseUpdate) PetscFunctionReturn(0);

  // store dense state matrix and transpose from tmpMatDen data structure
  ierr = MatConvert(tmpMatDen,MATSAME,MAT_INITIAL_MATRIX,&StateX);
  CHKERRQ(ierr);
  ierr = MatTranspose(StateX   ,      MAT_INITIAL_MATRIX,&StateXtrans);
  CHKERRQ(ierr);

  // print info
  ierr = MatDataInfo(tmpMatDen,"tmpMatDen");CHKERRQ(ierr);
  ierr = MatDataInfo(StateX,"State");CHKERRQ(ierr);
  ierr = MatDataInfo(StateXtrans,"StateTrans");CHKERRQ(ierr);

  // echo data
  ierr=PetscPrintf(PETSC_COMM_WORLD,"statecov=%f modelcov=%f \n",
                                     statecov,modelcov); CHKERRQ(ierr);

  // initialize covariance for measurements
  PetscScalar dummy = 1.0; // just need to initialize... will be overwritten
  ierr = InitializeSparseUncorrCov(CovR, dummy); CHKERRQ(ierr);
  ierr = MatDataInfo(CovR,"Init CovR");CHKERRQ(ierr);

  // initialize covariance for model
  ierr = InitializeDenseUncorrCov(CovQ, modelcov); CHKERRQ(ierr);
  ierr = MatDataInfo(CovQ,"Init CovQ");CHKERRQ(ierr);

  // initialize covariance for state
  ierr = InitializeDenseUncorrCov(CovP, statecov);CHKERRQ(ierr);
  ierr = MatDataInfo(CovP,"Init CovP");CHKERRQ(ierr);

  PetscFunctionReturn(0);
} 
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "SparseKalmanFilter::InitializeState"
/* 
   SparseInitializeState - Forms initial state and covariance

   Input Parameters:
   StateTemp - vector

 */
PetscErrorCode SparseKalmanFilter::InitializeState(Vec StateTemp)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  // initial conditions are for Temp difference
  ierr = VecSet(StateTemp,Setup::bodyTemp);CHKERRQ(ierr);

  // echo data
  ierr=PetscPrintf(PETSC_COMM_WORLD,"statecov=%f modelcov=%f \n",
                                     statecov,modelcov); CHKERRQ(ierr);

  // initialize covariance for measurements
  PetscScalar dummy = 1.0; // just need to initialize... will be overwritten
  ierr = InitializeSparseUncorrCov(CovR, dummy); CHKERRQ(ierr);
  ierr = MatDataInfo(CovR,"Init CovR");CHKERRQ(ierr);

  // initialize covariance for model
  ierr = InitializeSparseUncorrCov(CovQ, modelcov); CHKERRQ(ierr);
  ierr = MatDataInfo(CovQ,"Init CovQ");CHKERRQ(ierr);

  // initialize covariance for state
  ierr = InitializeSparseUncorrCov(CovP, statecov);CHKERRQ(ierr);
  ierr = MatDataInfo(CovP,"Init CovP");CHKERRQ(ierr);

  if(noDenseUpdate) PetscFunctionReturn(0);
  // store state and it's transpose
  //ierr = DAGetMatrix(dac,MATAIJ,&StateX);CHKERRQ(ierr);
  PetscInt M,N,m,n;
  ierr = MatGetSize(tmpMatDen,&M,&N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(tmpMatDen,&m,&n);CHKERRQ(ierr);
  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,m,n,M,N,n,PETSC_NULL,
                                       N-n,PETSC_NULL,&StateX);CHKERRQ(ierr);
  ierr = GetSparseMatrixFromTmp(StateX); CHKERRQ(ierr);
  ierr = MatTranspose(StateX   ,MAT_INITIAL_MATRIX,&StateXtrans);
  CHKERRQ(ierr);

  // print info
  ierr = MatDataInfo(tmpMatDen,"tmpMatDen");CHKERRQ(ierr);
  ierr = MatDataInfo(StateX,"State");CHKERRQ(ierr);
  ierr = MatDataInfo(StateXtrans,"StateTrans");CHKERRQ(ierr);

  PetscFunctionReturn(0);
} 
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::InitializeDenseUncorrCov"
PetscErrorCode KalmanFilter::InitializeDenseUncorrCov(Mat &CovMatrix,PetscScalar covariance)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscInt M,N,m,n;
  ierr = MatGetSize(tmpMatDen,&M,&N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(tmpMatDen,&m,&n);CHKERRQ(ierr);
  ierr = MatCreateMPIDense(PETSC_COMM_WORLD,m,n,M,N,
                           PETSC_NULL,&CovMatrix);CHKERRQ(ierr);

  /* Insert Ones and assemble to ensure full matrix is allocated */
  PetscInt Istart,Iend;
  ierr = MatGetOwnershipRange(CovMatrix,&Istart,&Iend);CHKERRQ(ierr);
  for (PetscInt Ii=Istart; Ii<Iend; Ii++) 
        MatSetValue(CovMatrix,Ii,Ii,covariance,INSERT_VALUES);
  ierr = MatAssemblyBegin(CovMatrix,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  CovMatrix,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr=PetscPrintf(PETSC_COMM_WORLD,
                   "(Dense) Initial Covariance Formed...\n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::InitializeSparseUncorrCov"
PetscErrorCode KalmanFilter::InitializeSparseUncorrCov(Mat &CovMatrix,PetscScalar covariance)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscInt       mx,my,mz;
  DAGetInfo(dac,0,&mx,&my,&mz,0,0,0,0,0,0,0);
  ierr = DAGetMatrix(dac,MATAIJ,   &CovMatrix  );CHKERRQ(ierr);

  // assume initially uncorrelated and
  // initialize covariance as a multiple of the identity matrix
  PetscScalar    v[7];
  MatStencil     col[7],row;
  for (int k=ProcStart[2]; k<ProcStart[2]+ProcWidth[2]; k++) {
    for (int j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) {
      for (int i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) {
        row.k = k; row.j = j; row.i = i;
        if (i == 0 || i == mx-1 || j == 0 || j == my-1 || k == 0 || k == mz-1) {
          v[0] = covariance;
          ierr = MatSetValuesStencil(CovMatrix,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
        } else {
        /* interior grid points */
          v[0] = 0.0       ; col[0].k=k-1;col[0].j=j;  col[0].i = i;
          v[1] = 0.0       ; col[1].k=k;  col[1].j=j-1;col[1].i = i;
          v[2] = 0.0       ; col[2].k=k;  col[2].j=j;  col[2].i = i-1;
          v[3] = covariance; col[3].k=row.k;col[3].j=row.j;col[3].i = row.i;
          v[4] = 0.0       ; col[4].k=k;  col[4].j=j;  col[4].i = i+1;
          v[5] = 0.0       ; col[5].k=k;  col[5].j=j+1;col[5].i = i;
          v[6] = 0.0       ; col[6].k=k+1;col[6].j=j;  col[6].i = i;
          ierr = MatSetValuesStencil(CovMatrix,1,&row,7,col,v,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
    }
  }

  /* Assemble matrices  */
  ierr = MatAssemblyBegin(CovMatrix,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  CovMatrix,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr=PetscPrintf(PETSC_COMM_WORLD,
                   "(Sparse) Initial Covariance Formed...\n");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             disp('Compute factorization needed for Kalman gain factor');
             %K = P*inv(P+R);

             disp('State Correction based on observation');
              
             define (P+R) q = z-x ==>  Pq = z-x-Rq
                  x = x + K*(z-x) 
                    = x + P*inv(P+R) *(z-x)
                    = x + P*q
                    = z - R*q 
             
             disp('Covariance Correction based on observation');
              
             define (P+R) S = P ==> -P*S = R*S - P
                  P = P - K*P  
                    = P - P*inv(P+R)*P
                    = P - P*S
                    = P + R*S - P
                    = R*S 

  Perform state update AND setup part of the update for the covariance update
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::StateUpdate"
PetscErrorCode KalmanFilter::StateUpdate(Vec State, Vec Measurement)
{

  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr=PetscPrintf(PETSC_COMM_WORLD,
                   "Update state from measurement data\n");CHKERRQ(ierr);
  
  // compute the kalman gain
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Kalman Gain...\n");CHKERRQ(ierr);

  // SumCov = P + R
  const MatType storagetype;
  ierr = MatGetType(SumCov,&storagetype);CHKERRQ(ierr);

  ierr = MatCopy(CovP,SumCov,DIFFERENT_NONZERO_PATTERN);
  if( !strcasecmp(storagetype,MATMPIDENSE)  || 
      !strcasecmp(storagetype,MATSEQDENSE) ){
    //  MatAYPX not playing w/ dense matrices
    ierr = MatAYPX(SumCov,1.0,CovR,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    //ierr = MatShift(SumCov,meascov);CHKERRQ(ierr);
  } else { // default is mpiiaj 
    ierr = MatAYPX(SumCov,1.0,CovR,SUBSET_NONZERO_PATTERN);CHKERRQ(ierr);
  }

  ierr = MatDataInfo(SumCov,  "(Update) SumCov");CHKERRQ(ierr);

  ierr=PetscPrintf(PETSC_COMM_WORLD,"Update state... \n");CHKERRQ(ierr);
  //globalVec = Measurement - State
  ierr = VecWAXPY(globalVec,-1.0,State,Measurement);CHKERRQ(ierr);

  // print info
  ierr = VecDataInfo(Measurement,"(Update) Measurement");CHKERRQ(ierr);
  ierr = VecDataInfo(globalVec,"(Update) globalVec");CHKERRQ(ierr);

  // SumCov* tmpVec = globalVec
  KSP            kspStateUpdate;         /* Krylov subspace method context */
  ierr = KSPCreate(PETSC_COMM_WORLD,&kspStateUpdate);CHKERRQ(ierr);
  ierr = KSPSetOperators(kspStateUpdate,SumCov,SumCov,
                         DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  // set defaults
  PC pcstate;
  ierr = KSPGetPC(kspStateUpdate,&pcstate); CHKERRQ(ierr);
  ierr = PCSetType(pcstate,PCBJACOBI); CHKERRQ(ierr);
  PetscInt maxiter = Setup::n_roi[0] * Setup::n_roi[1] * Setup::n_roi[2] / 10;
  ierr = KSPSetTolerances(kspStateUpdate,rtol0,PETSC_DEFAULT,PETSC_DEFAULT,
                          maxiter); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(kspStateUpdate);CHKERRQ(ierr);
  ierr = KSPSolve(kspStateUpdate,globalVec,tmpVec);CHKERRQ(ierr);
  ierr = KSPDestroy(kspStateUpdate);CHKERRQ(ierr); 

  // State= State + P * tmpVec  
  //ierr = MatMultAdd(CovP,tmpVec,State,State);CHKERRQ(ierr);
  //                            -OR-  

  // State= Measurement - R * tmpVec
  ierr = VecScale(tmpVec,-1.0);CHKERRQ(ierr);

  ierr = VecDataInfo(tmpVec,"(Update) tmpVec");CHKERRQ(ierr);
  ierr = MatDataInfo(CovR,  "(Update) CovR");CHKERRQ(ierr);
  ierr = VecDataInfo(State, "(Update) State");CHKERRQ(ierr);

  ierr = MatMultAdd(CovR,tmpVec,Measurement,State);CHKERRQ(ierr);

  ierr = VecDataInfo(State,"(Update) State");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
// subroutine to extract image data to a petsc vec
#undef __FUNCT__
#define __FUNCT__ "SparseKalmanFilter::CovarianceUpdate"
PetscErrorCode SparseKalmanFilter::CovarianceUpdate()
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  if(noDenseUpdate) PetscFunctionReturn(0); 

  ierr=PetscPrintf(PETSC_COMM_WORLD,
          "Update covariance from measurement data\n");CHKERRQ(ierr);

  // factor the matrix SumCov 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Factor Matrix SOLVER %s FACTOR %d ORDER %s...\n",
                             solvertype,factortype,ordertype);CHKERRQ(ierr);
  Mat SumCovFact;  // covariance of the state vector estimate
  ierr = MatGetFactor(SumCov,solvertype,
                             factortype,&SumCovFact); CHKERRQ(ierr);
  IS isrow,iscol;
  ierr = MatGetOrdering(SumCov,ordertype,&isrow,&iscol); CHKERRQ(ierr);
  ierr = MatLUFactorSymbolic(SumCovFact,SumCov,isrow,iscol,&factorInfo);CHKERRQ(ierr);
  ierr = MatLUFactorNumeric(SumCovFact,SumCov,&factorInfo);CHKERRQ(ierr);

  // P = SumCovFact * tmpMatDen
  Mat Ptmp;  // covariance of the state vector estimate
  ierr = MatConvert(CovP,MATMPIDENSE,MAT_INITIAL_MATRIX,&Ptmp);

  ierr=PetscPrintf(PETSC_COMM_WORLD,"MatMatSolve... \n");CHKERRQ(ierr);
  ierr = MatMatSolve(SumCovFact,Ptmp,tmpMatDen);CHKERRQ(ierr);
  // P = P * (I - tmpMatDen) = R * tmpMatDen
  PetscTruth  assembled ; 
  ierr = MatAssembled(tmpMatDen,&assembled); CHKERRQ(ierr);
  if(!assembled)
    { 
     ierr = MatAssemblyBegin(tmpMatDen,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
     ierr = MatAssemblyEnd(tmpMatDen,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }

  ierr = MatDataInfo(tmpMatDen,"(Update) tmpMatDen");CHKERRQ(ierr);
  // tmpMatDen = I - tmpMatDen
  //ierr= MatScale(tmpMatDen,-1.0);CHKERRQ(ierr);
  //ierr= MatShift(tmpMatDen, 1.0);CHKERRQ(ierr);
  
  // free memory
  ierr = MatDestroy(Ptmp      );CHKERRQ(ierr);
  ierr = MatDestroy(SumCovFact);CHKERRQ(ierr);

  // get sparse matrix from dense matrix
  Mat tmpMatAIJ;
  //ierr = DAGetMatrix(dac,MATAIJ,&tmpMatAIJ);CHKERRQ(ierr);
  PetscInt M,N,m,n;
  ierr = MatGetSize(tmpMatDen,&M,&N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(tmpMatDen,&m,&n);CHKERRQ(ierr);
  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,m,n,M,N,n,PETSC_NULL,
                          N-n,PETSC_NULL,&tmpMatAIJ);CHKERRQ(ierr);
  ierr = GetSparseMatrixFromTmp(tmpMatAIJ); CHKERRQ(ierr);

  // print info
  ierr = MatDataInfo(tmpMatDen,"(Sparse Update) tmpMatDen");CHKERRQ(ierr);
  ierr = MatDataInfo(tmpMatAIJ,"(Sparse Update) tmpMatAIJ");CHKERRQ(ierr);

  ierr = MatDataInfo(CovP,"(Sparse Update) CovP");CHKERRQ(ierr);
  ierr = MatDestroy(CovP);CHKERRQ(ierr); // MatMatMult will create P again
  ierr = PetscPrintf(PETSC_COMM_WORLD,"MatMatMult... \n");CHKERRQ(ierr);
  ierr = MatMatMult(CovR,tmpMatAIJ,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&CovP);
  ierr = MatDataInfo(CovP,"(Sparse Update) CovP");CHKERRQ(ierr);

  // free memory
  ierr = MatDestroy( tmpMatAIJ );CHKERRQ(ierr); 

  PetscFunctionReturn(0);
}
// update covariance from measurement data
#undef __FUNCT__
#define __FUNCT__ "DenseKalmanFilter::CovarianceUpdate"
PetscErrorCode DenseKalmanFilter::CovarianceUpdate()
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  if(noDenseUpdate) PetscFunctionReturn(0);

  // covariance of the state vector estimate
  ierr = MatDataInfo(CovP,"(Update) CovP");CHKERRQ(ierr);
  ierr=PetscPrintf(PETSC_COMM_WORLD,"MatMatMult... \n");CHKERRQ(ierr);
  ierr = MatMatMult(CovR,tmpMatDen,MAT_REUSE_MATRIX,PETSC_DEFAULT,&CovP);
  ierr = MatDataInfo(CovP,"(Update) CovP");CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
// update covariance from measurement data
#undef __FUNCT__
#define __FUNCT__ "UncorrelatedKalmanFilter::CovarianceUpdate"
PetscErrorCode UncorrelatedKalmanFilter::CovarianceUpdate()
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  Vec diagSumCov,diagCovP;

  // covariance of the state vector estimate
  ierr = MatDataInfo(CovP,"(Update) CovP");CHKERRQ(ierr);

  // create vectors 
  ierr = VecDuplicate(globalVec,&diagSumCov);
  ierr = VecDuplicate(globalVec,&diagCovP);

  // get diagonals
  ierr = MatGetDiagonal(SumCov,diagSumCov);
  ierr = MatGetDiagonal(CovP,diagCovP);

  // 1+\sigma^z/\sigma^x
  ierr = VecPointwiseDivide(diagSumCov,diagSumCov,diagCovP);

  // 1/(1+\sigma^z/\sigma^x)
  ierr = VecReciprocal(diagSumCov);

  // 1-\sigma^x
  ierr = VecScale(diagCovP,-1.0);
  ierr = VecShift(diagCovP, 1.0);

  // ( 1-\sigma^x )/( 1+\sigma^z/\sigma^x )
  ierr = VecPointwiseMult(diagSumCov,diagSumCov,diagCovP);

  // restore to matrix
  ierr = MatDiagonalSet(CovP,diagSumCov,INSERT_VALUES);
  ierr = MatDataInfo(CovP,"(Update) CovP");CHKERRQ(ierr);

  ierr = VecDestroy(diagSumCov);
  ierr = VecDestroy(diagCovP);

  PetscFunctionReturn(0);
}

// subroutine to extract image data to a petsc vec
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::ExtractImageData"
PetscErrorCode KalmanFilter::ExtractImageData(
                                  InputImageType::Pointer tempImage,
                                  InputImageType::RegionType &procRegion,
                                                           Vec DataBuffer)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // Software Guide : BeginLatex
  // setup real, imaginary, base phase, and temperature map iterators
  // The const slice iterator walks the 3D input image, and the non-const
  // linear iterator walks the 2D output image. The iterators are initialized
  // to walk the same linear path through a slice.  Remember that the
  // \emph{second} direction of the slice iterator defines the direction that
  // linear iteration walks within a slice

  InputIteratorType     tempIt(tempImage   , procRegion ); 
  //InputIteratorType   maskIt(maskImage   , maskImage->GetRequestedRegion());
  
  tempIt.SetFirstDirection(  0 );  tempIt.SetSecondDirection( 1 );
  //maskIt.SetFirstDirection(   0 );  maskIt.SetSecondDirection(  1 );
  
  // get pointer to petsc data structures
  PetscScalar   ****MapPixel;
  ierr = VecGetArray4d(DataBuffer,
            ProcWidth[2],ProcWidth[1],ProcWidth[0],nvarplot,
            ProcStart[2],ProcStart[1],ProcStart[0],0,&MapPixel); CHKERRQ(ierr);

  ierr=PetscPrintf(PETSC_COMM_WORLD,"Extracting Data...\n");CHKERRQ(ierr);
  /* loop through parallel data structures 
      tempIt.GetIndex() should contain processor wise indicies
  */
  tempIt.GoToBegin(); 
  while( !tempIt.IsAtEnd() )
    {
    const int kkk = tempIt.GetIndex()[2] - Setup::index[2];
    while ( !tempIt.IsAtEndOfSlice() )
      {
      const int jjj = tempIt.GetIndex()[1] - Setup::index[1];
      while ( !tempIt.IsAtEndOfLine() )
        {
        const int iii = tempIt.GetIndex()[0] - Setup::index[0];
        for(PetscInt ivar = 0; ivar < nvarplot ; ivar++) 
                  MapPixel[kkk][jjj][iii][ivar] = tempIt.Get() ;
        ++tempIt;
        }
      tempIt.NextLine();
      }
    tempIt.NextSlice();
    }
  //for (PetscInt k=ProcStart[2]; k<ProcStart[2]+ProcWidth[2]; k++) 
  //{
  //  for (PetscInt j=ProcStart[1]; j<ProcStart[1]+ProcWidth[1]; j++) 
  //  {
  //    for (PetscInt i=ProcStart[0]; i<ProcStart[0]+ProcWidth[0]; i++) 
  //    {
  //      for(PetscInt ivar = 0; ivar < nvarplot ; ivar++) 
  //                MapPixel[k][j][i][ivar] = tempIt.Get() ;
  //      ++tempIt; // update iterators
  //    }
  //    // get next line
  //    tempIt.NextLine(); 
  //  }
  //  // get next slice
  //  tempIt.NextSlice(); 
  //}
  CHKMEMQ; // check for memory corruption use -malloc_debug to enable

  /* Restore vector */
  ierr = VecRestoreArray4d(DataBuffer,
            ProcWidth[2],ProcWidth[1],ProcWidth[0],nvarplot,
            ProcStart[2],ProcStart[1],ProcStart[0],0,&MapPixel); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// subroutine to extract image data to a petsc vec
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::ApplyDirichletData"
PetscErrorCode KalmanFilter::ApplyDirichletData( Vec X, PetscScalar value)
{
  PetscErrorCode ierr=0;
  PetscFunctionBegin;

  if(this->dirichletNodes.size())
   {
    for (PetscInt Ii = 0 ; Ii < dirichletNodes.size(); Ii++ )
       ierr = VecSetValues (X, 1, &dirichletNodes[Ii], &value, INSERT_VALUES);
    ierr = VecAssemblyBegin(X);
    ierr = VecAssemblyEnd(X);
   }

  PetscFunctionReturn(ierr);
}
// subroutine to extract image data to a petsc vec
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::ExtractCovarianceData"
PetscErrorCode KalmanFilter::ExtractCovarianceData( const int istep,
                                  InputImageType::Pointer snrImage,
                                  InputImageType::Pointer varImage,
                                  InputImageType::RegionType &procRegion)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  Vec diagCovR;
  ierr = VecDuplicate(globalVec,&diagCovR);
  //ierr = VecDuplicate(globalVec,&tmpVec);

  // extract the data and take the max
  ierr = this->ExtractImageData(snrImage,procRegion,diagCovR); 
  //ierr = this->ExtractImageData(varImage,procRegion,tmpVec  ); 
  ierr = VecDataInfo(diagCovR,"(Extract Covariance) SNR");CHKERRQ(ierr);
  //ierr = VecDataInfo(tmpVec  ,"(Extract Covariance) Var");CHKERRQ(ierr);


  PetscTruth  overrideMax=PETSC_FALSE;
  ierr = PetscOptionsGetTruth(PETSC_NULL,"-override_max",&overrideMax,PETSC_NULL); 
  CHKERRQ(ierr);

  if(overrideMax)
    ierr = VecCopy(diagCovR,covMaxVec);
  else
    ierr = VecPointwiseMax(covMaxVec,diagCovR,covMaxVec);

  this->WriteDAImage(istep,"tmapstd",covMaxVec);

  ierr = MatDiagonalSet(CovR,covMaxVec,INSERT_VALUES);
  ierr = MatDataInfo(CovR,"(Extract Covariance) CovR");CHKERRQ(ierr);

  ierr = VecDestroy(diagCovR);

  PetscFunctionReturn(0);
}
// subroutine to write the image to disk
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::WriteDAImage"
PetscErrorCode KalmanFilter::WriteDAImage(const PetscInt istep, 
                                          const char *FileID, Vec Data)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // gather the image buffer on processor zero
  VecScatterBegin(gather,Data,imageVec,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(gather,Data,imageVec,INSERT_VALUES,SCATTER_FORWARD);

  // output various maps
  if(!Setup::rank) // write image from root process
   {
     // allocate memory for the output image
     VecOutputImageType::RegionType region;
     region.SetSize(  ROIsize  );
     VecOutputImageType::Pointer outputImage = VecOutputImageType::New();
     outputImage->SetRegions( region );
     outputImage->Allocate();
     // set image dimensions
     outputImage->SetSpacing(Setup::sp);
     outputImage->SetOrigin( ROIorgn );
       
     // initialize ITK iterators
     typedef itk::ImageSliceIteratorWithIndex< VecOutputImageType > 
                                                                 VecOutIterType;
     VecOutIterType  outputIt( outputImage, outputImage->GetRequestedRegion() );
     outputIt.SetFirstDirection(  0 ); outputIt.SetSecondDirection( 1 );
     outputIt.GoToBegin();

     // get pointer to petsc data structures
     PetscScalar   ****MapPixel;
     VecOutputImageType::PixelType   pixelValue;
     ierr = VecGetArray4d(imageVec,ROIsize[2],ROIsize[1],ROIsize[0],
                                   nvarplot, 0,0,0,0,&MapPixel); CHKERRQ(ierr);

     /*
        loop through parallel data structures
     */
     for (PetscInt k=0; k<ROIsize[2]; k++) 
     {
       for (PetscInt j=0; j<ROIsize[1]; j++) 
       {
         for (PetscInt i=0; i<ROIsize[0]; i++) 
         {
           outputIt.Set( MapPixel[k][j][i][0] ) ; 
           //for(PetscInt ivar = 0; ivar < nvarplot ; ivar++) 
           //            pixelValue[ivar] = MapPixel[k][j][i][ivar];
           //outputIt.Set( pixelValue ) ; 
           ++outputIt; // update iterators
         }
         outputIt.NextLine(); // get next line
       }
       outputIt.NextSlice(); // get next slice
     }
     ierr=VecRestoreArray4d(imageVec,ROIsize[2],ROIsize[1],ROIsize[0],
                                 nvarplot, 0,0,0,0,&MapPixel); CHKERRQ(ierr);

     // set filename 
     std::ostringstream output_filename;
     output_filename << Setup::OutputDir << "/"<< FileID << ".";
     OSSRealzeroright(output_filename,4,0,istep);
     output_filename << ".vtk" ;

     // setup writer
     VecWriterType::Pointer writer = VecWriterType::New();
     writer->SetFileName( output_filename.str() );

     // use custom class for proper data name
     itk::VTKImageVariableNameIO::Pointer ioBasePointer = 
                                          itk::VTKImageVariableNameIO::New();
     std::ostringstream annotatedVariableName;
     annotatedVariableName<< FileID << baseInfo::profileID ;
     ioBasePointer->SetVariableName( annotatedVariableName.str()  ) ;
     writer->SetImageIO( ioBasePointer ); // set variable name

     std::cout << "writing " << output_filename.str() << std::endl;
     writer->SetInput( outputImage );
     try
       {
       writer->Update();
       }
     catch (itk::ExceptionObject &ex)
       {
       std::cout << ex; abort();
       }
   }

  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "KalmanFilter::GetSparseMatrixFromTmp"
PetscErrorCode KalmanFilter::GetSparseMatrixFromTmp(Mat &SparseMat)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // store nonzero structure in AIJ format for later
 
  /* 
  PetscInt       m,N;
  Vec              x;
  PetscScalar    *xx,*data;
  ierr = MatGetArray(tmpMatDen,&xx);CHKERRQ(ierr);
  ierr = MatGetLocalSize(tmpMatDen,&m,PETSC_NULL);CHKERRQ(ierr);  // number local rows 
  ierr = MatGetSize(tmpMatDen,PETSC_NULL,&N);CHKERRQ(ierr);       // total columns in dense matrix
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,m,
                               PETSC_DETERMINE,PETSC_NULL,&x);CHKERRQ(ierr);
  for (PetscInt Jj=0; Jj<N; Jj++) {
    ierr = VecPlaceArray(x,xx + Jj*m);CHKERRQ(ierr);
    ierr = VecGetArray(   x,&data);CHKERRQ(ierr);
    for (PetscInt Ii=0; Ii<m; Ii++) {
      if( std::abs(data[Ii]) > zerotol ) { 
        MatSetValue(SparseMat     ,Ii,Jj,data[Ii],INSERT_VALUES);
      }
    }
    ierr = VecRestoreArray(x,&data);CHKERRQ(ierr);
    ierr = VecResetArray(x);CHKERRQ(ierr);
  }
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = MatRestoreArray(tmpMatDen,&xx);CHKERRQ(ierr);
  */

  const PetscScalar  *vwork;
  const PetscInt     *cwork;
  PetscInt Istart,Iend, nz;

  ierr=PetscPrintf(PETSC_COMM_WORLD,
                   "  Converting...\n");CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(tmpMatDen,&Istart,&Iend);CHKERRQ(ierr);
  for (PetscInt Ii=Istart; Ii<Iend; Ii++) {
    ierr = MatGetRow(tmpMatDen,Ii,&nz,&cwork,&vwork);CHKERRQ(ierr);
    for (PetscInt Jj=0; Jj<nz; Jj++) {
      //ierr = MatSetValues(SparseMat,1,&Ii,nz,cwork,vwork,INSERT_VALUES);CHKERRQ(ierr);
      if( std::abs(vwork[Jj]) > zerotol ) { 
        MatSetValue(SparseMat ,Ii,cwork[Jj],vwork[Jj],INSERT_VALUES);
      }
    }
    ierr = MatRestoreRow(tmpMatDen,Ii,&nz,&cwork,&vwork);CHKERRQ(ierr);
  }

  ierr=PetscPrintf(PETSC_COMM_WORLD,
                   "  Assembling...\n");CHKERRQ(ierr);
  ierr = MatAssemblyBegin(SparseMat ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(  SparseMat ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);


  PetscFunctionReturn(0);
}
