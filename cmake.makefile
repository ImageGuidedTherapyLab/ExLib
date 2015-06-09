.PHONY: tags 

# The following are the definitions for each target individually.
Makefile: CMakeLists.txt
	cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_VERBOSE_MAKEFILE=ON  -DOTB_DIR=$(OTB_DIR)

clean:
	rm -rf CMakeCache.txt Makefile CMakeFiles/ ITKIOFactoryRegistration/ cmake_install.cmake

tags:
	#ctags -R  --langmap=c++:+.cu --langmap=c++:+.cuh $(MATLABROOT) . 
	ctags -R  --langmap=c++:+.cu --langmap=c++:+.cuh --langmap=c++:+.txx --langmap=c++:+.cl $(ITK_SOURCE) .
