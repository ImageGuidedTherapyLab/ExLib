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

int main(int argc,char **argv) 
{
  SNES           snes; 
  Vec            x,f;  
  Mat            J;
  DM             da;
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
  PetscInt globalSize = 4;
  ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,-globalSize,-globalSize,-globalSize,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&da);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da,&x); VecDuplicate(x,&f);CHKERRQ(ierr);
  ierr = DMCreateMatrix(da,MATAIJ,&J);CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,f,ComputeFunction,da);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,J,J,ComputeJacobian,da);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = ComputeFunction(snes,x,f,(void *)da);
  //ierr = SNESSolve(snes,PETSC_NULL,x);CHKERRQ(ierr);

  ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&f);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}

struct StarStencil
{
  PetscInt       m_rank,m_deviceNum; //device info
  PetscInt       m_xs,m_ys,m_zs,m_xm,m_ym,m_zm; //corners
  PetscScalar    m_hxhzdhy,m_hyhzdhx,m_hxhydhz;
  
  StarStencil(PetscInt rank, PetscInt deviceNum,
              PetscInt  xs,PetscInt  ys,PetscInt  zs, 
              PetscInt  xm,PetscInt  ym,PetscInt  zm,
              PetscScalar hxhzdhy,PetscScalar hyhzdhx,PetscScalar hxhydhz) : 
               m_rank(rank),m_deviceNum(deviceNum),
               m_xs(xs),m_ys(ys),m_zs(zs), 
               m_xm(xm),m_ym(ym),m_zm(zm), 
               m_hxhzdhy(m_hxhzdhy),m_hyhzdhx(hyhzdhx),m_hxhydhz(hxhydhz) {}

	template <typename Tuple>
	__host__ __device__
	void operator()(Tuple t)
	{
		/* f = (2*u_i - u_(i+1) - u_(i-1))/h - h*exp(u_i) */
	     thrust::get<0>(t) = 1;
             PetscInt Iz = thrust::get<8>(t)/m_ym/m_xm;
             PetscInt Iy = (thrust::get<8>(t)-Iz*m_ym*m_xm)/m_xm;
             PetscInt Ix = (thrust::get<8>(t)-Iz*m_ym*m_xm- Iy*m_xm);
             // print launch parameters and dbg info
             printf("rank=%d device=%d blockDim=(%d,%d,%d) gridDim=(%d,%d,%d) warpSize=%d blockIdx=(%d,%d,%d) threadIdx=(%d,%d,%d) size=(%d,%d,%d) globalID=%d index=(%d,%d,%d)\n",m_rank,m_deviceNum,blockDim.x, blockDim.y, blockDim.z, gridDim.x, gridDim.y, gridDim.z, warpSize,blockIdx.x,blockIdx.y,blockIdx.z,threadIdx.x,threadIdx.y,threadIdx.z,m_xm,m_ym,m_zm,thrust::get<8>(t),Ix,Iy,Iz);
             if (Ix > 0  && Ix < m_xm-1) {
               thrust::get<0>(t) = (2.0*thrust::get<1>(t) - thrust::get<2>(t) - thrust::get<3>(t)) / m_hxhzdhy - m_hxhydhz*exp(thrust::get<1>(t));
             } else if (Ix == 0) {
               thrust::get<0>(t) = thrust::get<1>(t) / m_hxhzdhy;
             } else if (Ix == m_xm-1) {
               thrust::get<0>(t) = thrust::get<1>(t) / m_hxhzdhy;
             } 
		
	}
};

PetscErrorCode ComputeFunction(SNES snes,Vec u,Vec f,void *ctx) 
{
  PetscInt       i,j,k,GlobalDAMx,GlobalDAMy,GlobalDAMz,xs,xm,ys,ym,zs,zm;
  PetscInt       ustartshift,uendshift,xoffset,yoffset,zoffset,fstart;
  PetscScalar    ***uu,***ff,hx,hy,hz, hxhzdhy,hyhzdhx,hxhydhz;
  PetscScalar    u_val,u_east,u_west,u_north,u_south,u_up, u_down, u_xx, u_yy,u_zz,sc ,two =2.0;
  DM             da = (DM) ctx; 
  Vec            ulocal;
  PetscErrorCode ierr;
  PetscMPIInt    rank,size;
  MPI_Comm       comm;
  CUSPARRAY      *uarray,*farray;
  PetscLogEventBegin(LogFunction,0,0,0,0); // init libMesh

  ierr = DMDAGetInfo(da,PETSC_IGNORE,&GlobalDAMx,&GlobalDAMy,&GlobalDAMz,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  hx     = 1.0/(PetscReal)(GlobalDAMx-1);
  hy     = 1.0/(PetscReal)(GlobalDAMy-1);
  hz     = 1.0/(PetscReal)(GlobalDAMz-1);
  hxhzdhy = hx*hz/hy;
  hyhzdhx = hy*hz/hx;
  hxhydhz = hx*hy/hz;
  sc     = hx*hy*hz*3.0;
  ierr = DMGetLocalVector(da,&ulocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,u,INSERT_VALUES,ulocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,u,INSERT_VALUES,ulocal);CHKERRQ(ierr);

  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  if (useCUSP) {
    StarStencil  stencil_op(0,0,xs,ys,zs,xm,ym,zm,hxhzdhy,hyhzdhx,hxhydhz);// transformation operator
    ierr = VecCUSPGetArrayRead(ulocal,&uarray);CHKERRQ(ierr);
    ierr = VecCUSPGetArrayWrite(f,&farray);CHKERRQ(ierr);
    ierr = PetscObjectGetComm((PetscObject)da,&comm);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
    if (rank) ustartshift = 1; else ustartshift = 0;
    if (rank != size-1) uendshift = 1; else uendshift = 0;
    xoffset = 1;
    yoffset = xm;
    zoffset = xm*ym;
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
		       stencil_op);
      
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
    
    /* Compute function over the locally owned part of the grid */
    for (k=zs; k<zs+zm; k++) {
      for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {
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
  PetscInt       i,Mx,xm,xs; 
  PetscScalar    hx,*xx; 
  Vec            xlocal;
  PetscErrorCode ierr;

  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  hx = 1.0/(PetscReal)(Mx-1);
  ierr = DMGetLocalVector(da,&xlocal);DMGlobalToLocalBegin(da,x,INSERT_VALUES,xlocal);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,x,INSERT_VALUES,xlocal);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,xlocal,&xx);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,PETSC_NULL,PETSC_NULL,&xm,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

  ierr = MatZeroEntries(*J);CHKERRQ(ierr);
  ierr = MatShift(*J,1.0);CHKERRQ(ierr);
  for (i=xs; i<xs+xm; i++) {
    if (i == 0 || i == Mx-1) { 
      ierr = MatSetValue(*J,i,i,1.0/hx,INSERT_VALUES);CHKERRQ(ierr);
    } else {
      ierr = MatSetValue(*J,i,i-1,-1.0/hx,INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatSetValue(*J,i,i,2.0/hx - hx*PetscExpScalar(xx[i]),INSERT_VALUES);CHKERRQ(ierr);
      ierr = MatSetValue(*J,i,i+1,-1.0/hx,INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  *flag = SAME_NONZERO_PATTERN;
  ierr = DMDAVecRestoreArray(da,xlocal,&xx);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&xlocal);CHKERRQ(ierr);
  return 0;}

