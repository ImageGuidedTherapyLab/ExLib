        #include <thrust/version.h>
        #include <cusp/version.h>
        #include <iostream>
        
        int main(void)
        {
            int thrust_major = THRUST_MAJOR_VERSION;
            int thrust_minor = THRUST_MINOR_VERSION;
        
            int cusp_major = CUSP_MAJOR_VERSION;
            int cusp_minor = CUSP_MINOR_VERSION;
            int cusp_subminor = CUSP_SUBMINOR_VERSION;

            int ierr=cudaDeviceReset();
            std::cout << "Thrust v" << thrust_major << "." << thrust_minor << std::endl;
            std::cout << "Cusp   v" << cusp_major << "." << cusp_minor << "." << cusp_subminor << std::endl;
            int numdev = 0;
            ierr=cudaGetDeviceCount(&numdev);
            std::cout << "ierr" << ierr  << std::endl;
            std::cout << "dev " << numdev<< std::endl;
        
            return 0;
        }

