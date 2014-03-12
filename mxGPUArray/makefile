#
# file : Makefile (UNIX)
#
#  mexGPUArray   from 
# 
#     $(MATLABROOT)/toolbox/distcomp/gpu/extern/src/mex/mexGPUExample.cu
#
# You can invoke this Makefile using 
#  make -f Makefile MATLABROOT=[directory where MATLAB is installed]
#
# If you do not want to specify the MATLABROOT at the gmake command line 
# every single time, you can define it by uncommenting the line below
# and assigning the correct root of MATLAB (with no trailing '/') on
# the right hand side.
#
# MATLABROOT	:= /opt/apps/MATLAB/R2013a/
#

#
# Defaults
#

MEX=$(MATLABROOT)/bin/mex

# The following are the definitions for each target individually.

mexGPUExample:   mexGPUExample.cu
	$(MEX) -v -f glnxa64/mexopts.sh $^
