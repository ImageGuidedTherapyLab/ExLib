static char help[] = "Solves -Laplacian u - exp(u) = 0,  0 < x < 1 using GPU\n\n";
/*
   Same as ex47.c except it also uses the GPU to evaluate the function
*/

#include <petscdmda.h>
#include <petscsnes.h>
#include <petsccusp.h>
#include "cusp/detail/device/utils.h"

extern PetscErrorCode ComputeFunction(SNES,Vec,Vec,void*), ComputeJacobian(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
PetscBool  useCUSP = PETSC_FALSE;
PetscLogEvent LogFunction = 0;
__device__ PetscInt *cudaTest;

struct StarStencil
{
  PetscInt       m_rank,m_deviceNum; //device info
  PetscInt       m_xs,m_ys,m_zs,m_xm,m_ym,m_zm; //corners
  PetscScalar    m_hx,m_hy,m_hz;
  PetscScalar    m_x0,m_y0,m_z0;
  PetscScalar    m_density           ;
  PetscScalar    m_specificheat      ;
  PetscScalar    m_deltat            ;
  PetscScalar    m_conduction        ;
  PetscScalar    m_bloodspecificheat ;
  PetscScalar    m_perfusion         ;
  PetscScalar    m_bodytemp          ;
  
  StarStencil(PetscInt rank, PetscInt deviceNum,
              PetscInt  xs,PetscInt  ys,PetscInt  zs, 
              PetscInt  xm,PetscInt  ym,PetscInt  zm,
              PetscScalar hx,PetscScalar hy,PetscScalar hz) : 
               m_rank(rank),m_deviceNum(deviceNum),
               m_xs(xs),m_ys(ys),m_zs(zs), 
               m_xm(xm),m_ym(ym),m_zm(zm), 
               m_hx(hx),m_hy(hy),m_hz(hz) 
               {
                  m_density           = 1.e3;
                  m_specificheat      = 3.8e3;
                  m_deltat            = 1.00;
                  m_conduction        = 0.57;
                  m_bloodspecificheat = 3.4e3;
                  m_perfusion         =  6.0;
                  m_bodytemp          = 37.0;
                  m_x0          = 0.005;
                  m_y0          = 0.005;
                  m_z0          = 0.005;
               }

	template <typename Tuple>
	__host__ __device__
	void operator()(Tuple t)
	{
		/* f = (2*u_i - u_(i+1) - u_(i-1))/h - h*exp(u_i) */
	     thrust::get<0>(t) = 1;
             PetscInt Iz = thrust::get<8>(t)/m_ym/m_xm;
             PetscInt Iy = (thrust::get<8>(t)-Iz*m_ym*m_xm)/m_xm;
             PetscInt Ix = (thrust::get<8>(t)-Iz*m_ym*m_xm- Iy*m_xm);
             PetscScalar sc      = m_hx*m_hz*m_hy;
             PetscScalar hxhzdhy = m_hx*m_hz/m_hy;
             PetscScalar hyhzdhx = m_hy*m_hz/m_hx;
             PetscScalar hxhydhz = m_hx*m_hy/m_hz;
             PetscScalar two     = 2.0;
             PetscScalar u_val       = thrust::get<1>(t) ;//1  u(i  ,j  ,k  )
             // print launch parameters and dbg info
             // printf("rank=%d device=%d blockDim=(%d,%d,%d) gridDim=(%d,%d,%d) warpSize=%d blockIdx=(%d,%d,%d) threadIdx=(%d,%d,%d) size=(%d,%d,%d) globalID=%d index=(%d,%d,%d)\n",m_rank,m_deviceNum,blockDim.x, blockDim.y, blockDim.z, gridDim.x, gridDim.y, gridDim.z, warpSize,blockIdx.x,blockIdx.y,blockIdx.z,threadIdx.x,threadIdx.y,threadIdx.z,m_xm,m_ym,m_zm,thrust::get<8>(t),Ix,Iy,Iz);
             if (
                 Ix > 0  && Ix < m_xm-1
                         &&
                 Iy > 0  && Iy < m_ym-1
                         &&
                 Iz > 0  && Iz < m_zm-1
                ) {
               // decode the tuple
               PetscScalar u_east      = thrust::get<2>(t) ;//2  u(i+1,j  ,k  )
               PetscScalar u_west      = thrust::get<3>(t) ;//3  u(i-1,j  ,k  )
               PetscScalar u_north     = thrust::get<4>(t) ;//4  u(i  ,j+1,k  )
               PetscScalar u_south     = thrust::get<5>(t) ;//5  u(i  ,j-1,k  )
               PetscScalar u_up        = thrust::get<6>(t) ;//6  u(i  ,j  ,k+1)
               PetscScalar u_down      = thrust::get<7>(t) ;//7  u(i  ,j  ,k-1)
               PetscScalar u_xx        = (-u_east  + two*u_val - u_west )*hyhzdhx;
               PetscScalar u_yy        = (-u_north + two*u_val - u_south)*hxhzdhy;
               PetscScalar u_zz        = (-u_up    + two*u_val - u_down )*hxhydhz;
               PetscScalar sqdist      = (m_hx * Ix - m_x0)*(m_hx * Ix - m_x0)
                                       + (m_hy * Iy - m_y0)*(m_hy * Iy - m_y0)
                                       + (m_hz * Iz - m_z0)*(m_hz * Iz - m_z0);
               PetscScalar source      = 1.e4 * exp(5.0/(sqdist +1.0));
               thrust::get<0>(t) = sc * ( source
                              + m_density*m_specificheat/m_deltat* u_val 
                              + m_bloodspecificheat*m_perfusion*(m_bodytemp - 0.5*u_val) ) 
                              + m_conduction/2.0* (u_xx + u_yy + u_zz) ;
             } else { // dirichlet bc everywhere else
               thrust::get<0>(t) = u_val;
             } 
		
	}
};
int main(int argc,char **argv) 
{
  SNES           snes; 
  Vec            x,f;  
  Mat            J;
  PetscErrorCode ierr;
  cudaError      ierrCuda;
  char           *tmp,typeName[256];
  int            myrank;
  PetscBool      flg;

  PetscInitialize(&argc,&argv,(char *)0,help);

  MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
  int deviceNum=myrank;
  {
    int deviceCount;
    CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));
    
    ierr = PetscPrintf(PETSC_COMM_SELF, "!!!!!found %d devices !!!!!\n",deviceCount);CHKERRQ(ierr);
    if (deviceCount == 0) {
      ierr = PetscPrintf(PETSC_COMM_SELF, "!!!!!No devices found!!!!!\n");CHKERRQ(ierr);
      return -1000;
    }

    if (deviceNum >= deviceCount || deviceNum < 0) {
      ierr = PetscPrintf(PETSC_COMM_SELF, "\n!!!!!Invalid GPU number %d given hence default gpu %d will be used !!!!!\n", deviceNum, 0);CHKERRQ(ierr);
      deviceNum = 0;
    }
  }

  ierrCuda =  cudaSetDevice(deviceNum);
  if (ierrCuda != cudaSuccess) {
    ierr = PetscPrintf(PETSC_COMM_SELF, " cuda Error: %s , exiting\n",cudaGetErrorString( ierrCuda));CHKERRQ(ierr);
    return -1;
  }
  ierr = PetscPrintf(PETSC_COMM_SELF, " reseting GPU: \n");CHKERRQ(ierr);
  CUDA_SAFE_CALL(cudaDeviceReset());

  ierr = PetscPrintf(PETSC_COMM_SELF, "Running on...\n\n");CHKERRQ(ierr);
  cudaDeviceProp deviceProp;
  if (cudaGetDeviceProperties(&deviceProp, deviceNum) == cudaSuccess) {
    ierr = PetscPrintf(PETSC_COMM_SELF, " Device %d: %s %d.%d\n", deviceNum, deviceProp.name,deviceProp.major,deviceProp.minor);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF," Global memory available on device in bytes %d\n"                            ,  deviceProp.totalGlobalMem                  );
    ierr = PetscPrintf(PETSC_COMM_SELF," Shared memory available per block in bytes %d\n"                            ,  deviceProp.sharedMemPerBlock               );
    ierr = PetscPrintf(PETSC_COMM_SELF," 32-bit registers available per block %d\n"                                  ,  deviceProp.regsPerBlock                    );
    ierr = PetscPrintf(PETSC_COMM_SELF," Warp size in threads %d\n"                                                  ,  deviceProp.warpSize                        );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum pitch in bytes allowed by memory copies %d\n"                       ,  deviceProp.memPitch                        );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum number of threads per block %d\n"                                   ,  deviceProp.maxThreadsPerBlock              );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum size of each dimension of a block %d\n"                             ,  deviceProp.maxThreadsDim[0]                );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum size of each dimension of a block %d\n"                             ,  deviceProp.maxThreadsDim[1]                );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum size of each dimension of a block %d\n"                             ,  deviceProp.maxThreadsDim[2]                );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum size of each dimension of a grid %d\n"                              ,  deviceProp.maxGridSize[0]                  );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum size of each dimension of a grid %d\n"                              ,  deviceProp.maxGridSize[1]                  );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum size of each dimension of a grid %d\n"                              ,  deviceProp.maxGridSize[2]                  );
    ierr = PetscPrintf(PETSC_COMM_SELF," Clock frequency in kilohertz %d\n"                                          ,  deviceProp.clockRate                       );
    ierr = PetscPrintf(PETSC_COMM_SELF," Constant memory available on device in bytes %d\n"                          ,  deviceProp.totalConstMem                   );
    ierr = PetscPrintf(PETSC_COMM_SELF," Alignment requirement for textures %d\n"                                    ,  deviceProp.textureAlignment                );
    ierr = PetscPrintf(PETSC_COMM_SELF," Number of multiprocessors on device %d\n"                                   ,  deviceProp.multiProcessorCount             );
    ierr = PetscPrintf(PETSC_COMM_SELF," Specified whether there is a run time limit on kernels %d\n"                ,  deviceProp.kernelExecTimeoutEnabled        );
    ierr = PetscPrintf(PETSC_COMM_SELF," Device is integrated as opposed to discrete %d\n"                           ,  deviceProp.integrated                      );
    ierr = PetscPrintf(PETSC_COMM_SELF," Device can map host memory with cudaHostAlloc/cudaHostGetDevicePointer %d\n",  deviceProp.canMapHostMemory                );
    ierr = PetscPrintf(PETSC_COMM_SELF," Compute mode (See ::cudaComputeMode) %d\n"                                  ,  deviceProp.computeMode                     );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum 1D texture size %d\n"                                               ,  deviceProp.maxTexture1D                    );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum 2D texture dimensions %d\n"                                         ,  deviceProp.maxTexture2D[0]                 );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum 2D texture dimensions %d\n"                                         ,  deviceProp.maxTexture2D[1]                 );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum 3D texture dimensions %d\n"                                         ,  deviceProp.maxTexture3D[0]                 );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum 3D texture dimensions %d\n"                                         ,  deviceProp.maxTexture3D[1]                 );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum 3D texture dimensions %d\n"                                         ,  deviceProp.maxTexture3D[2]                 );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum 1D layered texture dimensions %d\n"                                 ,  deviceProp.maxTexture1DLayered[0]          );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum 1D layered texture dimensions %d\n"                                 ,  deviceProp.maxTexture1DLayered[1]          );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum 2D layered texture dimensions %d\n"                                 ,  deviceProp.maxTexture2DLayered[0]          );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum 2D layered texture dimensions %d\n"                                 ,  deviceProp.maxTexture2DLayered[1]          );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum 2D layered texture dimensions %d\n"                                 ,  deviceProp.maxTexture2DLayered[2]          );
    ierr = PetscPrintf(PETSC_COMM_SELF," Alignment requirements for surfaces %d\n"                                   ,  deviceProp.surfaceAlignment                );
    ierr = PetscPrintf(PETSC_COMM_SELF," Device can possibly execute multiple kernels concurrently %d\n"             ,  deviceProp.concurrentKernels               );
    ierr = PetscPrintf(PETSC_COMM_SELF," Device has ECC support enabled %d\n"                                        ,  deviceProp.ECCEnabled                      );
    ierr = PetscPrintf(PETSC_COMM_SELF," PCI bus ID of the device %d\n"                                              ,  deviceProp.pciBusID                        );
    ierr = PetscPrintf(PETSC_COMM_SELF," PCI device ID of the device %d\n"                                           ,  deviceProp.pciDeviceID                     );
    ierr = PetscPrintf(PETSC_COMM_SELF," PCI domain ID of the device %d\n"                                           ,  deviceProp.pciDomainID                     );
    ierr = PetscPrintf(PETSC_COMM_SELF," 1 if device is a Tesla device using TCC driver, 0 otherwise %d\n"           ,  deviceProp.tccDriver                       );
    ierr = PetscPrintf(PETSC_COMM_SELF," Number of asynchronous engines %d\n"                                        ,  deviceProp.asyncEngineCount                );
    ierr = PetscPrintf(PETSC_COMM_SELF," Device shares a unified address space with the host %d\n"                   ,  deviceProp.unifiedAddressing               );
    ierr = PetscPrintf(PETSC_COMM_SELF," Peak memory clock frequency in kilohertz %d\n"                              ,  deviceProp.memoryClockRate                 );
    ierr = PetscPrintf(PETSC_COMM_SELF," Global memory bus width in bits %d\n"                                       ,  deviceProp.memoryBusWidth                  );
    ierr = PetscPrintf(PETSC_COMM_SELF," Size of L2 cache in bytes %d\n"                                             ,  deviceProp.l2CacheSize                     );
    ierr = PetscPrintf(PETSC_COMM_SELF," Maximum resident threads per multiprocessor %d\n"                           ,  deviceProp.maxThreadsPerMultiProcessor     );
  } else {
    ierr = PetscPrintf(PETSC_COMM_SELF, " Unable to determine device %d properties, exiting\n",deviceNum);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, " cuda Error: %s , exiting\n",cudaGetErrorString( ierrCuda));CHKERRQ(ierr);
    return -1;
  }

  PetscLogEventRegister("ComputeFunction",0,&LogFunction); 
  ierr = PetscOptionsGetString(PETSC_NULL,"-da_vec_type",typeName,256,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscStrstr(typeName,"cusp",&tmp);CHKERRQ(ierr);
    if (tmp) useCUSP = PETSC_TRUE;
  }

  size_t sizeIndex = 3 * sizeof(PetscInt);
  CUDA_SAFE_CALL(cudaMalloc((void **) &cudaTest, sizeIndex));   // Allocate array on device

  //ierr = DMDACreate1d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,-8,1,1,PETSC_NULL,&da);CHKERRQ(ierr);
  PetscInt globalSize = 125;
  globalSize = 99;
  DM             da;
  ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,-globalSize,-globalSize,-globalSize,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&da);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da,&x); VecDuplicate(x,&f);CHKERRQ(ierr);
  if (useCUSP)
    {
     ierr = DMCreateMatrix(da,MATAIJCUSP,&J);CHKERRQ(ierr);
    }
  else
    {
     ierr = DMCreateMatrix(da,MATAIJ,&J);CHKERRQ(ierr);
    }

  PetscInt       GlobalDAMx,GlobalDAMy,GlobalDAMz,xs,xm,ys,ym,zs,zm;
  PetscScalar    hx,hy,hz;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&GlobalDAMx,&GlobalDAMy,&GlobalDAMz,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  hx     = 1.0/(PetscReal)(GlobalDAMx-1);
  hy     = 1.0/(PetscReal)(GlobalDAMy-1);
  hz     = 1.0/(PetscReal)(GlobalDAMz-1);
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  StarStencil  stencil_op(0,0,xs,ys,zs,xm,ym,zm,hx,hy,hz);// transformation operator
  ierr = DMSetApplicationContext(da,&stencil_op);CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,f,ComputeFunction,da);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,J,J,ComputeJacobian,da);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = ComputeFunction(snes,x,f,(void *)da);
  ierr = SNESSolve(snes,PETSC_NULL,x);CHKERRQ(ierr);

  ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&f);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  // call device reset to flush buffer
  CUDA_SAFE_CALL(cudaDeviceReset());
  PetscFinalize();
  return 0;
}


PetscErrorCode ComputeFunction(SNES snes,Vec u,Vec f,void *ctx) 
{
  PetscInt       i,j,k;
  PetscInt       ustartshift,uendshift,xoffset,yoffset,zoffset,fstart;
  PetscScalar    ***uu,***ff,hxhzdhy,hyhzdhx,hxhydhz;
  PetscScalar    u_val,u_east,u_west,u_north,u_south,u_up, u_down, u_xx, u_yy,u_zz,sc ,two =2.0;
  DM             da = (DM) ctx; 
  Vec            ulocal;
  PetscErrorCode ierr;
  PetscMPIInt    rank,size;
  MPI_Comm       comm;
  CUSPARRAY      *uarray,*farray;
  PetscLogEventBegin(LogFunction,0,0,0,0); // init libMesh

  ierr = DMGetLocalVector(da,&ulocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,u,INSERT_VALUES,ulocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,u,INSERT_VALUES,ulocal);CHKERRQ(ierr);
  StarStencil  *stencil_op;
  ierr = DMGetApplicationContext(da,(void *)&stencil_op);CHKERRQ(ierr);
  hxhzdhy = stencil_op->m_hx*stencil_op->m_hz/stencil_op->m_hy;
  hyhzdhx = stencil_op->m_hy*stencil_op->m_hz/stencil_op->m_hx;
  hxhydhz = stencil_op->m_hx*stencil_op->m_hy/stencil_op->m_hz;
  sc      = stencil_op->m_hx*stencil_op->m_hy*stencil_op->m_hz*3.0;

  if (useCUSP) {
    ierr = VecCUSPGetArrayRead(ulocal,&uarray);CHKERRQ(ierr);
    ierr = VecCUSPGetArrayWrite(f,&farray);CHKERRQ(ierr);
    ierr = PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
    if (rank) ustartshift = 1; else ustartshift = 0;
    if (rank != size-1) uendshift = 1; else uendshift = 0;
    xoffset = 1;
    yoffset = stencil_op->m_xm;
    zoffset = stencil_op->m_xm*stencil_op->m_ym;
    ierr = VecGetOwnershipRange(f,&fstart,PETSC_NULL);CHKERRQ(ierr);
    try {
      
      // typedef these iterators for shorthand
      thrust::for_each(
		       thrust::make_zip_iterator(
						 thrust::make_tuple(
            farray->begin(),                              //0
            uarray->begin()+ustartshift,                  //1  u(i  ,j  ,k  )
            uarray->begin()+ustartshift + xoffset,        //2  u(i+1,j  ,k  )
            uarray->begin()+ustartshift - xoffset,        //3  u(i-1,j  ,k  )
            uarray->begin()+ustartshift + yoffset,        //4  u(i  ,j+1,k  )
            uarray->begin()+ustartshift - yoffset,        //5  u(i  ,j-1,k  )
            uarray->begin()+ustartshift + zoffset,        //6  u(i  ,j  ,k+1)
            uarray->begin()+ustartshift - zoffset,        //7  u(i  ,j  ,k-1)
            thrust::counting_iterator<int>(fstart)        //8
                                                                    )), 
		       thrust::make_zip_iterator(
						 thrust::make_tuple(
            farray->end(),                            //0
            //farray->begin()+10,                            //0
            uarray->end()+uendshift,                  //1  u(i  ,j  ,k  )
            uarray->end()+uendshift + xoffset,        //2  u(i+1,j  ,k  )
            uarray->end()+uendshift - xoffset,        //3  u(i-1,j  ,k  )
            uarray->end()+uendshift + yoffset,        //4  u(i  ,j+1,k  )
            uarray->end()+uendshift - yoffset,        //5  u(i  ,j-1,k  )
            uarray->end()+uendshift + zoffset,        //6  u(i  ,j  ,k+1)
            uarray->end()+uendshift - zoffset,        //7  u(i  ,j  ,k-1)
            thrust::counting_iterator<int>(fstart) + u->map->n        //8
                                                                    )),
		       *stencil_op);
      
      PetscInt hostTest[3]={-1,-1,-1};
      //CUDA_SAFE_CALL(cudaMemcpy(hostTest, cudaTest,3*sizeof(PetscInt),cudaMemcpyDeviceToHost));
      ierr = PetscPrintf(PETSC_COMM_WORLD, "%d %d %d \n",hostTest[0],hostTest[1],hostTest[2]);CHKERRQ(ierr);
    }
    catch(char* all){
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Thrust is not working\n");CHKERRQ(ierr);
    }
    ierr = VecCUSPRestoreArrayRead(ulocal,&uarray);CHKERRQ(ierr);
    ierr = VecCUSPRestoreArrayWrite(f,&farray);CHKERRQ(ierr);
  } else {
    ierr = DMDAVecGetArray(da,ulocal,&uu);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,f,&ff);CHKERRQ(ierr);
    
    PetscInt       GlobalDAMx,GlobalDAMy,GlobalDAMz;
    ierr = DMDAGetInfo(da,PETSC_IGNORE,&GlobalDAMx,&GlobalDAMy,&GlobalDAMz,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
    /* Compute function over the locally owned part of the grid */
    for (k=stencil_op->m_zs; k<stencil_op->m_zs+stencil_op->m_zm; k++) {
      for (j=stencil_op->m_ys; j<stencil_op->m_ys+stencil_op->m_ym; j++) {
        for (i=stencil_op->m_xs; i<stencil_op->m_xs+stencil_op->m_xm; i++) {
          if (i == 0 || j == 0 || k == 0 || i == GlobalDAMx-1 || j == GlobalDAMy-1 || k == GlobalDAMz-1) {
            ff[k][j][i] = uu[k][j][i];
          } else {
            u_val       = uu[k][j][i];
            u_east      = uu[k][j][i+1];
            u_west      = uu[k][j][i-1];
            u_north     = uu[k][j+1][i];
            u_south     = uu[k][j-1][i];
            u_up        = uu[k+1][j][i];
            u_down      = uu[k-1][j][i];
            u_xx        = (-u_east  + two*u_val - u_west )*hyhzdhx;
            u_yy        = (-u_north + two*u_val - u_south)*hxhzdhy;
            u_zz        = (-u_up    + two*u_val - u_down )*hxhydhz;
            ff[k][j][i]  = u_xx + u_yy + u_zz - sc*PetscExpScalar(u_val);
          }
        }
      }
    }
    ierr = DMDAVecRestoreArray(da,ulocal,&uu);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,f,&ff);CHKERRQ(ierr);
  }
  ierr = DMRestoreLocalVector(da,&ulocal);CHKERRQ(ierr);
  PetscLogEventEnd(LogFunction,0,0,0,0);   // init libMesh
  //VecView(u,0);printf("f\n");
  //VecView(f,0);
  return 0;

}
PetscErrorCode ComputeJacobian(SNES snes,Vec x,Mat *J,Mat *B,MatStructure *flag,void *ctx)
{
  DM             da = (DM) ctx; 
  Vec            xlocal;
  PetscErrorCode ierr;

  ierr = DMGetLocalVector(da,&xlocal);DMGlobalToLocalBegin(da,x,INSERT_VALUES,xlocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,x,INSERT_VALUES,xlocal);CHKERRQ(ierr);

  PetscInt       GlobalDAMx,GlobalDAMy,GlobalDAMz,xs,xm,ys,ym,zs,zm;
  PetscScalar    hx,hy,hz;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&GlobalDAMx,&GlobalDAMy,&GlobalDAMz,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  hx     = 1.0/(PetscReal)(GlobalDAMx-1);
  hy     = 1.0/(PetscReal)(GlobalDAMy-1);
  hz     = 1.0/(PetscReal)(GlobalDAMz-1);
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  PetscScalar    hxhzdhy,hyhzdhx,hxhydhz,sc;
  hxhzdhy = hx*hz/hy;
  hyhzdhx = hy*hz/hx;
  hxhydhz = hx*hy/hz;
  sc      = hx*hy*hz*3.0;

  ierr = MatZeroEntries(*J);CHKERRQ(ierr);
  ierr = MatShift(*J,1.0);CHKERRQ(ierr);

  StarStencil  *stencil_op;
  ierr = DMGetApplicationContext(da,(void *)&stencil_op);CHKERRQ(ierr);

  /* Compute function over the locally owned part of the grid */
  PetscScalar    v[7],two = 2.0;
  MatStencil     col[7],row;
  PetscInt       i,j,k;
  for (k=zs; k<zs+zm; k++) {
    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
        row.k = k; row.j = j; row.i = i;
        if (i > 0 && j > 0 && k > 0 && i < GlobalDAMx-1 && j < GlobalDAMy-1 && k < GlobalDAMz-1) {
          v[0] = -0.5 * stencil_op->m_conduction * hxhydhz; col[0].k=k-1;col[0].j=j;  col[0].i = i;
          v[1] = -0.5 * stencil_op->m_conduction * hxhzdhy; col[1].k=k;  col[1].j=j-1;col[1].i = i;
          v[2] = -0.5 * stencil_op->m_conduction * hyhzdhx; col[2].k=k;  col[2].j=j;  col[2].i = i-1;
          v[3] =  sc*(  stencil_op->m_density*stencil_op->m_specificheat/stencil_op->m_deltat 
                    + 0.5 * stencil_op->m_perfusion * stencil_op->m_bloodspecificheat) 
                    + 1.0 * stencil_op->m_bloodspecificheat * (hyhzdhx+hxhzdhy+hxhydhz);
                           col[3].k=row.k;col[3].j=row.j;col[3].i = row.i;
          v[4] = -0.5 * stencil_op->m_conduction * hyhzdhx; col[4].k=k;  col[4].j=j;  col[4].i = i+1;
          v[5] = -0.5 * stencil_op->m_conduction * hxhzdhy; col[5].k=k;  col[5].j=j+1;col[5].i = i;
          v[6] = -0.5 * stencil_op->m_conduction * hxhydhz; col[6].k=k+1;col[6].j=j;  col[6].i = i;
          ierr = MatSetValuesStencil(*J,1,&row,7,col,v,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
    }
  }

  ierr = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  *flag = SAME_NONZERO_PATTERN;
  ierr = DMRestoreLocalVector(da,&xlocal);CHKERRQ(ierr);
  return 0;}

