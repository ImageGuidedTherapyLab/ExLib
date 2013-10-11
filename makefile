OPENCL_INCLUDE=/opt/apps/khronos/1.0/
hello: hello.c
	g++ -o hello -I$(OPENCL_INCLUDE) -g -O0 hello.c -lOpenCL
tags: 
	ctags -R $(OPENCL_INCLUDE)/CL/

