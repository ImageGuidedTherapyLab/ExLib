# intel
OPENCLINCLUDE=/opt/apps/intel/opencl-1.2-3.0.67279/include 
OPENCLLIB=/opt/apps/intel/opencl-1.2-3.0.67279/lib64
CXX = icpc
# mac
OPENCLINCLUDE=/Applications/Xcode.app/Contents//Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk//System/Library/Frameworks/OpenCL.framework/
OPENCLLIB=/Applications/Xcode.app/Contents//Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk//System/Library/Frameworks/OpenCL.framework/
CXX = gcc
# ubuntu 12.04
OPENCLINCLUDE=/usr/include
OPENCLLIB=/usr/lib
CXX = g++


hello: hello.c
	$(CXX) -o hello -I$(OPENCLINCLUDE) -g -O0 hello.c  -Wl,-rpath,$(OPENCLLIB) -L$(OPENCLLIB) -lOpenCL
tags: 
	ctags -R $(OPENCLINCLUDE)

