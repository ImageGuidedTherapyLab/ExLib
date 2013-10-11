OPENCLDIR=/opt/apps/intel/opencl-1.2-3.0.67279
CXX = icpc

hello: hello.c
	$(CXX) -o hello -I$(OPENCLDIR)/include -g -O0 hello.c  -Wl,-rpath,$(OPENCLDIR)/lib64/ -L$(OPENCLDIR)/lib64/ -lOpenCL
tags: 
	ctags -R $(OPENCL_INCLUDE)/CL/

