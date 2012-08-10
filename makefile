#
# file : Makefile (UNIX)
#
#  Example of a stand-alone makefile for compiler-generated M files
#
#  Copyright 1997-2005 The MathWorks, Inc.
#  $Revision: 1.7.6.4 $  $Date: 2009/05/18 19:50:05 $
#

#
# You can invoke this Makefile using 
#  gmake -f Makefile MATLABROOT=[directory where MATLAB is installed]
#
# If you do not want to specify the MATLABROOT at the gmake command line 
# every single time, you can define it by uncommenting the line below
# and assigning the correct root of MATLAB (with no trailing '/') on
# the right hand side.
#
# MATLABROOT	:=
#

#
# Defaults
#

MCC=$(MATLABROOT)/bin/mcc
MBUILD=$(MATLABROOT)/bin/mbuild

#
# It's important that MRANK gets built before MRANKP, since MRANK has
# a dependency on MRANK.M to build MRANK.C as does MRANKP.  however,
# MRANKP doesn't build MRANK.C only MRANK does.  so it's important
# that MRANK builds first, so that MRANK.C is up to date by the time
# MRANKP builds with it
#
C_TARGETS=mrank driver magicsquare mrankp multargp matrixdriver  collect flames triangle ex_loadsave

CPP_TARGETS=matrixdriverp trianglep

STAND_ALONE_TARGETS=$(C_TARGETS) $(CPP_TARGETS)


# Use the following target to build and test all the executables
bnt_exes: build_all
	@echo "*********** testing mrank ***********"
	./mrank 
	@echo "*********** testing driver ***********"
	./driver 
	@echo "*********** testing flames ***********"
	./flames 
	@echo "*********** testing magicsquare 4 ***********"
	./magicsquare 4
	@echo "*********** testing mrankp ***********"
	./mrankp 
	@echo "*********** testing multargp ***********"
	./multargp 
	@echo "*********** testing matrixdriver ***********"
	./matrixdriver 
	@echo "*********** testing triangle ***********"
	./triangle 
	@echo "*********** testing collect ***********"
	./collect 
	@echo "*********** testing trianglep ***********"
	./trianglep 
	@echo "*********** testing matrixdriverp ***********"
	./matrixdriverp
	@echo "*********** testing ex_loadsave ***********"
	./ex_loadsave



# Use this target to only build all the targets
build_all: $(STAND_ALONE_TARGETS)

# Use this target to only build all the C++ targets
only_cpp: $(CPP_TARGETS)

# Use this target to only build the C targets
only_c: $(C_TARGETS)


# The following are the definitions for each target individually.

mrank:  main.m mrank.m
	$(MCC) -m $^ -o $@
	@rm -f *.o*

driver:  driver.m 
	$(MCC) -m -R -nojvm -R -nodisplay -v $^ -o $@
	#@rm -f *.o*

flames:  flames.m flames.mat
	$(MCC) -m flames.m -a flames.mat
	@rm -f *.o*

magicsquare:  magicsquare.m
	$(MCC) -m $^ -o magicsquare
	@rm -f *.o*

ex_loadsave:  
	$(MCC) -m "./Data_Handling/ex_loadsave.m" -a "./Data_Handling/user_data.mat" -a "./Data_Handling/userdata/extra_data.mat" -a "./externdata/extern_data.mat" 
	@rm -f *.o*

# If you do not wish to create a shared library, use the single line MCC command.
# This command will create the wrapper files and then compile all the source files
# including the wrapper code and the user supplied C/C++ code into a single executable
# binary. No shared library is created in this process.  

triangle:  triangle.c main_for_lib.c main_for_lib.h sierpinski.m
#	$(MCC)  -l -o libtriangle sierpinski.m
#	$(MBUILD) triangle.c main_for_lib.c  -L. -ltriangle
	$(MCC)  -W lib:libtriangle -T link:exe triangle.c main_for_lib.c  sierpinski.m
	@rm -f *.o*


# If you wish to create a shared library first and then link matrixdriver.c with the 
# resulting library, comment the first line of the execution code for the target below
# and uncomment the second and the third line. See the comments for the triangle target
# above for more details.

matrixdriver:  matrixdriver.c addmatrix.m eigmatrix.m multiplymatrix.m
	$(MCC) -W lib:libmatrix -T link:exe $^
#	$(MCC)  -l -o libmatrix addmatrix.m eigmatrix.m multiplymatrix.m
#	$(MBUILD) matrixdriver.c -L. -lmatrix
	@rm -f *.o*


multargp:  multargp.c  main_for_lib.c multarg.m printmatrix.m
	$(MCC) -W lib:libMultpkg -T link:exe $^
	@rm -f *.o*


mrankp:  mrankp.c main_for_lib.c mrank.m printmatrix.m
	$(MCC) -W lib:libPkg -T link:exe $^
	@rm -f *.o*

# The next two targets are for the CPP shared library target example
trianglep:  triangle.cpp main_for_lib.c main_for_lib.h  sierpinski.m
	$(MCC) -B cpplib:libtrianglep sierpinski.m
	$(MBUILD) triangle.cpp main_for_lib.c -L. -ltrianglep -output trianglep
	@rm -f *.o*

matrixdriverp:  matrixdriver.cpp addmatrix.m eigmatrix.m multiplymatrix.m
	$(MCC) -B cpplib:libmatrixp addmatrix.m eigmatrix.m multiplymatrix.m
	$(MBUILD) matrixdriver.cpp -L. -lmatrixp -output matrixdriverp
	@rm -f *.o*


# The following example depicts %#external pragma.
collect:  collect.m measure.c
	$(MCC) -m $^


clean:
	rm -f $(STAND_ALONE_TARGETS)

reallyclean:
	rm -rf $(STAND_ALONE_TARGETS)				\
		mrank.h mrank.c                   		\
		main_main.c main.h main.c         		\
		driver.h driver.c hello_main.c      		\
		flames.h flames.c flames_main.c			\
		magicsquare.h magicsquare.c magicsquare_main.c	\
		multarg.h multarg.c		  		\
		libPkg* libMult* libtriangle* libmatrix*	\
		*.xml *.ctf *_mcr *mcc_component_data.c 		\
		*_main.c *.exports lib*                         \
		collect_one_external.h mccExcludedFIles.log  \
		output	
