
# This directory contains SNES example programs for solving systems of
# nonlinear equations.  

#CPPFLAGS	 = -I/PETSc3/geodynamics/PetscSimulationsViewers/src
CFLAGS           =
FFLAGS		 =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR		 = src/snes/examples/tutorials/
MANSEC           = SNES
EXAMPLESC	 = ex1.c ex2.c ex3.c ex4.c ex5.c ex5s.c ex6.c ex7.c ex8.c \
                ex10.c ex12.c ex14.c ex15.c ex18.c ex19.c ex20.c ex21.c ex22.c \
                ex25.c ex28.c ex30.c \
                ex35.c ex38.c ex42.c ex43.c ex44.c ex45.c ex46.c ex48.c \
                ex51.c ex52.c ex53.c ex54.c ex55.c ex58.c ex59.c ex60.c \
                ex61.c ex61gen.c ex61view.c ex62.c ex63.c ex633D_DB.c ex64.c ex65.c
EXAMPLESF	 = ex1f.F ex5f.F ex5f90.F ex5f90t.F ex5fs.F ex39f90.F ex40f90.F90
EXAMPLESCH	 = ex43-44.h
EXAMPLESFH       = ex5f.h
EXAMPLESMATLAB   = ex5m.m  ex61genm.m ex61m.m
DIRS		 = ex10d

all: FocusUltraSoundModel
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

FocusUltraSoundModel: FocusUltraSoundModel.o  chkopts
	-${CLINKER}  -o FocusUltraSoundModel FocusUltraSoundModel.o ${PETSC_SNES_LIB}
	${RM} FocusUltraSoundModel.o

tags:
	ctags -f petsc_tags -R --langmap=c++:+.inl --langmap=c:+.cu --langmap=python:+.pxi --langmap=c++:+.txx --languages=c,c++,python,fortran $(PETSC_DIR) $(CUDA_DIR) .

ex1: ex1.o  chkopts
	-${CLINKER} -o ex1 ex1.o ${PETSC_SNES_LIB}
	${RM} ex1.o
ex1f: ex1f.o  chkopts
	-${FLINKER} -o ex1f ex1f.o ${PETSC_SNES_LIB}
	${RM} ex1f.o
ex1f90: ex1f90.o  chkopts
	-${FLINKER} -o ex1f90 ex1f90.o ${PETSC_SNES_LIB}
	${RM} ex1f90.o
ex2: ex2.o  chkopts
	-${CLINKER} -o ex2 ex2.o ${PETSC_SNES_LIB}
	${RM} ex2.o
ex3: ex3.o  chkopts
	-${CLINKER} -o ex3 ex3.o ${PETSC_SNES_LIB}
	${RM} ex3.o
ex4: ex4.o  chkopts
	-${CLINKER} -o ex4 ex4.o ${PETSC_SNES_LIB}
	${RM} ex4.o
ex5: ex5.o chkopts
	-${CLINKER} -o ex5 ex5.o ${PETSC_SNES_LIB}
	${RM} ex5.o
ex5f: ex5f.o  chkopts
	-${FLINKER} -o ex5f ex5f.o ${PETSC_SNES_LIB}
	${RM} ex5f.o
#
#  The SGI parallelizing compiler generates incorrect code by treating 
#  the math functions (such as sqrt and exp) as local variables. The 
#  sed below patches this.
#
ex5s: chkopts
	@if [ "${PETSC_ARCH}" != "IRIX64" ]; then echo "Only for PETSC_ARCH of IRIX64"; false ; fi
	-${CC} -pca keep  -WK,-lo=l ${FCONF} ${CFLAGS} -c ex5s.c
	sed "s/, sqrt/ /g"   ex5s.M | sed "s/, exp/ /g"  > ex5s_tmp.c
	-${CC} -mp ${PCC_FLAGS} ${CFLAGS} ${CCPPFLAGS} -c ex5s_tmp.c
	-${FC} -pfa keep -mp -64 ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS} -c ex5fs.F
	-${CLINKER} -mp -o ex5s ex5s_tmp.o ex5fs.o ${PETSC_SNES_LIB}
	${RM} ex5s.o
ex5f90: ex5f90.o  chkopts
	-${FLINKER} -o ex5f90 ex5f90.o ${PETSC_SNES_LIB}
	${RM} ex5f90.o
ex5f90t: ex5f90t.o  chkopts
	-${FLINKER} -o ex5f90t ex5f90t.o ${PETSC_SNES_LIB}
	${RM} ex5f90t.o
ex6: ex6.o  chkopts
	-${CLINKER} -o ex6 ex6.o ${PETSC_SNES_LIB}
	${RM} ex6.o
ex7: ex7.o  chkopts
	-${CLINKER} -o ex7 ex7.o ${PETSC_SNES_LIB}
	${RM} ex7.o
ex8: ex8.o chkopts
	-${CLINKER} -o ex8 ex8.o  ${PETSC_SNES_LIB}
	${RM} ex8.o
ex9: ex9.o chkopts
	-${CLINKER} -o ex9 ex9.o ${PETSC_SNES_LIB}
	${RM} ex9.o
ex10: ex10.o chkopts
	-${CLINKER} -o ex10 ex10.o ${PETSC_SNES_LIB}
	${RM} ex10.o
ex12: ex12.o chkopts
	-${CLINKER} -o ex12 ex12.o ${PETSC_SNES_LIB}
	${RM} ex12.o
ex13: ex13.o chkopts
	-${CLINKER} -o ex13 ex13.o ${PETSC_SNES_LIB}
	${RM} ex13.o
ex14: ex14.o chkopts
	-${CLINKER} -o ex14 ex14.o ${PETSC_SNES_LIB}
	${RM} ex14.o
ex15: ex15.o chkopts
	-${CLINKER} -o ex15 ex15.o ${PETSC_SNES_LIB}
	${RM} ex15.o
ex16: ex16.o chkopts
	-${CLINKER} -o ex16 ex16.o ${PETSC_SNES_LIB}
	${RM} ex16.o
ex17: ex17.o chkopts
	-${CLINKER} -o ex17 ex17.o ${PETSC_SNES_LIB}
	${RM} ex17.o
ex18: ex18.o chkopts
	-${CLINKER} -o ex18 ex18.o ${PETSC_SNES_LIB}
	${RM} ex18.o
ex19:  ex19.o chkopts
	-${CLINKER} -o ex19  ex19.o ${PETSC_SNES_LIB}
	${RM} ex19.o
ex19tu:  ex19tu.o chkopts
	-${CLINKER} -o ex19tu  ex19tu.o ${PETSC_SNES_LIB}
	${RM} ex19tu.o
ex20: ex20.o chkopts
	-${CLINKER} -o ex20 ex20.o ${PETSC_SNES_LIB}
	${RM} ex20.o
ex21: ex21.o chkopts
	-${CLINKER} -o ex21 ex21.o ${PETSC_SNES_LIB}
	${RM} ex21.o
ex22: ex22.o chkopts
	-${CLINKER} -o ex22 ex22.o ${PETSC_SNES_LIB}
	${RM} ex22.o
ex23: ex23.o chkopts
	-${CLINKER} -o ex23 ex23.o ${PETSC_SNES_LIB}
	${RM} ex23.o
ex24: ex24.o chkopts
	-${CLINKER} -o ex24 ex24.o ${PETSC_SNES_LIB}
	${RM} ex24.o
ex25: ex25.o chkopts
	-${CLINKER} -o ex25 ex25.o ${PETSC_SNES_LIB}
	${RM} ex25.o
ex26: ex26.o chkopts
	-${CLINKER} -o ex26 ex26.o ${PETSC_SNES_LIB}
	${RM} ex26.o
ex27: ex27.o chkopts
	-${CLINKER} -o ex27 ex27.o ${PETSC_SNES_LIB}
	${RM} ex27.o
ex28: ex28.o chkopts
	-${CLINKER} -o ex28 ex28.o ${PETSC_SNES_LIB}
	${RM} ex28.o
ex29: ex29.o chkopts
	-${CLINKER} -o ex29 ex29.o ${PETSC_SNES_LIB}
	${RM} ex29.o
ex30: ex30.o chkopts
	-${CLINKER} -o ex30 ex30.o ${PETSC_SNES_LIB}
	${RM} ex30.o
ex31: ex31.o chkopts
	-${CLINKER} -o ex31 ex31.o ${PETSC_SNES_LIB}
	${RM} ex31.o
ex32: ex32.o chkopts
	-${CLINKER} -o ex32 ex32.o ${PETSC_SNES_LIB}
	${RM} ex32.o
ex35: ex35.o chkopts
	-${CLINKER} -o ex35 ex35.o ${PETSC_SNES_LIB}
	${RM} ex35.o
ex38: ex38.o  chkopts
	-${CLINKER} -o ex38 ex38.o ${PETSC_SNES_LIB}
	${RM} ex38.o
ex39f90: ex39f90.o  chkopts
	-${FLINKER} -o ex39f90 ex39f90.o ${PETSC_SNES_LIB}
	${RM} ex39f90.o
ex40f90: ex40f90.o  chkopts
	-${FLINKER} -o ex40f90 ex40f90.o ${PETSC_SNES_LIB}
	${RM} ex40f90.o
ex41: ex41.o  chkopts
	-${CLINKER} -o ex41 ex41.o ${PETSC_SNES_LIB}
	${RM} ex41.o
ex42: ex42.o  chkopts
	-${CLINKER} -o ex42 ex42.o ${PETSC_SNES_LIB}
	${RM} ex42.o
ex43: ex43.o  chkopts
	-${CLINKER} -o ex43 ex43.o ${PETSC_SNES_LIB}
	${RM} ex43.o
ex44: ex44.o  chkopts
	-${CLINKER} -o ex44 ex44.o ${PETSC_SNES_LIB}
	${RM} ex44.o
ex45: ex45.o  chkopts
	-${CLINKER} -o ex45 ex45.o ${PETSC_SNES_LIB}
	${RM} ex45.o
ex46: ex46.o  chkopts
	-${CLINKER} -o ex46 ex46.o ${PETSC_SNES_LIB}
	${RM} ex46.o
ex48: ex48.o  chkopts
	-${CLINKER} -o ex48 ex48.o ${PETSC_SNES_LIB}
	${RM} ex48.o
ex51: ex51.o  chkopts
	-${CLINKER} -o ex51 ex51.o ${PETSC_SNES_LIB}
	${RM} ex51.o
ex52: ex52.o ex52_integrateElement.o chkopts
	${PYTHON} ${PETSC_DIR}/bin/pythonscripts/PetscGenerateFEMQuadrature.py 2 1 1 1 laplacian ex52.h
	-${CLINKER} -o ex52 ex52.o ex52_integrateElement.o  ${PETSC_SNES_LIB}
	${RM} ex52.o ex52_integrateElement.o
ex53: ex53.o  chkopts
	-${CLINKER} -o ex53 ex53.o ${PETSC_SNES_LIB}
	${RM} ex53.o
ex54: ex54.o  chkopts
	-${CLINKER} -o ex54 ex54.o ${PETSC_SNES_LIB}
	${RM} ex54.o
ex55: ex55.o  chkopts
	-${CLINKER} -o ex55 ex55.o ${PETSC_SNES_LIB}
	${RM} ex55.o
ex56: ex56.o  chkopts
	-${CLINKER} -o ex56 ex56.o ${PETSC_SNES_LIB}
	${RM} ex56.o
ex57: ex57.o  chkopts
	-${CLINKER} -o ex57 ex57.o ${PETSC_SNES_LIB}
	${RM} ex57.o
ex58: ex58.o  chkopts
	-${CLINKER} -o ex58 ex58.o ${PETSC_SNES_LIB}
	${RM} ex58.o
ex59: ex59.o  chkopts
	-${CLINKER} -o ex59 ex59.o ${PETSC_SNES_LIB}
	${RM} ex59.o
ex60: ex60.o  chkopts
	-${CLINKER} -o ex60 ex60.o ${PETSC_SNES_LIB}
	${RM} ex60.o
ex61: ex61.o  chkopts
	-${CLINKER} -o ex61 ex61.o ${PETSC_SNES_LIB}
	${RM} ex61.o
ex61gen: ex61gen.o  chkopts
	-${CLINKER} -o ex61gen ex61gen.o ${PETSC_SNES_LIB}
	${RM} ex61gen.o
ex61view: ex61view.o  chkopts
	-${CLINKER} -o ex61view ex61view.o ${PETSC_SNES_LIB}
	${RM} ex61view.o
ex62: ex62.o  chkopts
	-${CLINKER} -o ex62 ex62.o ${PETSC_SNES_LIB}
	${RM} ex62.o
ex63: ex63.o  chkopts
	-${CLINKER} -o ex63 ex63.o ${PETSC_SNES_LIB}
	${RM} ex63.o
ex64: ex64.o  chkopts
	-${CLINKER} -o ex64 ex64.o ${PETSC_SNES_LIB}
	${RM} ex64.o
ex65: ex65.o  chkopts
	-${CLINKER} -o ex65 ex65.o ${PETSC_SNES_LIB}
	${RM} ex65.o
ex653D: ex653D.o chkopts
	-${CLINKER} -o ex653D ex653D.o ${PETSC_SNES_LIB}
	${RM} ex653D.o
ex65dm: ex65dm.o  chkopts
	-${CLINKER} -o ex65dm ex65dm.o ${PETSC_SNES_LIB}
	${RM} ex65dm.o
ex633D_DB: ex633D_DB.o  chkopts
	-${CLINKER} -o ex633D_DB ex633D_DB.o ${PETSC_SNES_LIB}
	${RM} ex633D_DB.o
#--------------------------------------------------------------------------
runex1:
	-@${MPIEXEC} -n 1 ./ex1 -ksp_gmres_cgs_refinement_type refine_always -snes_monitor_short > ex1_1.tmp 2>&1;   \
	   if (${DIFF} output/ex1_1.out ex1_1.tmp) then true; \
	   else echo  ${PWD} "\nPossible problem with with ex1_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex1_1.tmp
runex1_X:
	-@${MPIEXEC} -n 1 ./ex1 -ksp_monitor_short -ksp_type gmres -ksp_gmres_krylov_monitor -snes_monitor_short > ex1_X.tmp 2>&1;   \
	   if (${DIFF} output/ex1_X.out ex1_X.tmp) then true; \
	   else echo  ${PWD} "\nPossible problem with with ex1_X, diffs above \n========================================="; fi; \
	   ${RM} -f ex1_X.tmp
runex1_2:
	-@${MPIEXEC} -n 2 ./ex1 -snes_mf_operator -snes_fd -snes_monitor_short > ex1_2.tmp 2>&1;   \
	   if (${DIFF} output/ex1_2.out ex1_2.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex1_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex1_2.tmp
runex1f:
	-@${MPIEXEC} -n 1 ./ex1f -ksp_gmres_cgs_refinement_type refine_always -snes_monitor_short > ex1f_1.tmp 2>&1;   \
	   if (${DIFF} output/ex1f_1.out ex1f_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex1f_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex1f_1.tmp
runex2:
	-@${MPIEXEC} -n 1 ./ex2 -nox -snes_monitor_cancel -snes_monitor_short -snes_view -pc_type jacobi -ksp_gmres_cgs_refinement_type refine_always >ex2_1.tmp 2>&1;	\
	   if (${DIFF} output/ex2_1.out ex2_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex2_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex2_1.tmp
runex2_2:
	-@${MPIEXEC} -n 1 ./ex2 -nox -snes_monitor_cancel -snes_monitor_short -snes_type tr -snes_view -pc_type jacobi -ksp_gmres_cgs_refinement_type refine_always > ex2_2.tmp 2>&1; \
	   if (${DIFF} output/ex2_2.out ex2_2.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex2_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex2_2.tmp
runex2_3:
	-@${MPIEXEC} -n 1 ./ex2 -nox -snes_monitor_cancel -snes_monitor_short -malloc no -snes_view -pc_type jacobi -ksp_gmres_cgs_refinement_type refine_always >ex2_3.tmp 2>&1;	\
	   if (${DIFF} output/ex2_1.out ex2_3.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex2_3, diffs above \n========================================="; fi; \
	   ${RM} -f ex2_3.tmp
runex3:
	-@${MPIEXEC} -n 1 ./ex3 -nox -snes_monitor_cancel -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex3_1.tmp 2>&1;	  \
	   if (${DIFF} output/ex3_1.out ex3_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex3_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex3_1.tmp
runex3_2:
	-@${MPIEXEC} -n 3 ./ex3 -nox -pc_type asm -mat_type mpiaij -snes_monitor_cancel -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex3_2.tmp 2>&1; \
	   if (${DIFF} output/ex3_2.out ex3_2.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex3_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex3_2.tmp
runex3_3:
	-@${MPIEXEC} -n 2 ./ex3 -nox -snes_monitor_cancel -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex3_3.tmp 2>&1; \
	   if (${DIFF} output/ex3_3.out ex3_3.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex3_3, diffs above \n========================================="; fi; \
	   ${RM} -f ex3_3.tmp
runex3_4:
	-@${MPIEXEC} -n 1 ./ex3 -nox -pre_check_iterates -post_check_iterates > ex3_4.tmp 2>&1; \
	   if (${DIFF} output/ex3_4.out ex3_4.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex3_4, diffs above \n========================================="; fi; \
	   ${RM} -f ex3_4.tmp
runex4:
	-@${MPIEXEC} -n 1 ./ex4 -snes_monitor_short -snes_grid_sequence 1 -pc_type mg -pc_mg_type full -mg_levels_ksp_type gmres -mg_levels_ksp_max_it 1 -mg_levels_pc_type ilu -ksp_type fgmres -snes_converged_reason > ex4_1.tmp 2>&1;	  \
	   if (${DIFF} output/ex4_1.out ex4_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex4_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex4_1.tmp
runex5:
	-@${MPIEXEC} -n 1 ./ex5 -pc_type mg -ksp_monitor_short  -snes_view -pc_mg_levels 3 -pc_mg_galerkin -da_grid_x 17 -da_grid_y 17 -mg_levels_ksp_monitor_short -mg_levels_ksp_norm_type unpreconditioned -snes_monitor_short -mg_levels_ksp_chebyshev_estimate_eigenvalues 0.5,1.1 -mg_levels_pc_type sor -pc_mg_type full > ex5_1.tmp 2>&1; \
	   if (${DIFF} output/ex5_1.out ex5_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_1.tmp
runex5_2:
	-@${MPIEXEC} -n 1 ./ex5 -pc_type mg -ksp_monitor_short  -snes_view -pc_mg_galerkin -snes_grid_sequence 3 -mg_levels_ksp_monitor_short -mg_levels_ksp_norm_type unpreconditioned -snes_monitor_short -mg_levels_ksp_chebyshev_estimate_eigenvalues 0.5,1.1 -mg_levels_pc_type sor -pc_mg_type full > ex5_2.tmp 2>&1; \
	   if (${DIFF} output/ex5_2.out ex5_2.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_2.tmp
runex5_3:
	-@${MPIEXEC} -n 2 ./ex5 -snes_grid_sequence 2 -snes_mf_operator -snes_converged_reason -snes_view -pc_type mg > ex5_3.tmp 2>&1; \
	   if (${DIFF} output/ex5_3.out ex5_3.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5_3, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_3.tmp
runex5_4:
	-@${MPIEXEC} -n 2 ./ex5 -snes_grid_sequence 2 -snes_monitor_short -ksp_monitor_short -ksp_converged_reason -snes_converged_reason -snes_view -pc_type mg > ex5_4.tmp 2>&1; \
	   if (${DIFF} output/ex5_4.out ex5_4.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5_4, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_4.tmp

#COMPOSABLE SOLVER DEMOS: OPTIONS
CSD_BASIC_COMMAND_LINE = -@${MPIEXEC} -n 1 ./ex5 -da_grid_x 81 -da_grid_y 81 -snes_monitor_short -snes_max_it 50 -par 6.0
N_SMOOTHS = 3
N_RESTART = 10

runex5_5_ngmres:
	-@${CSD_BASIC_COMMAND_LINE} -snes_type ngmres -snes_ngmres_m ${N_RESTART} > ex5_5_ngmres.tmp 2>&1; \
	   if (${DIFF} output/ex5_5_ngmres.out ex5_5_ngmres.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5_5_ngmres, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_5_ngmres.tmp

runex5_5_ngmres_nrichardson:
	-@${CSD_BASIC_COMMAND_LINE} -snes_type ngmres -snes_ngmres_m ${N_RESTART} -npc_snes_type nrichardson -npc_snes_max_it ${N_SMOOTHS} \
        > ex5_5_ngmres_richardson.tmp 2>&1; \
	   if (${DIFF} output/ex5_5_ngmres_richardson.out ex5_5_ngmres_richardson.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5_5_ngmres_richardson, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_5_ngmres_richardson.tmp

runex5_5_ncg:
	-@${CSD_BASIC_COMMAND_LINE} -snes_type ncg -snes_ncg_type fr \
        > ex5_5_ncg.tmp 2>&1; \
	   if (${DIFF} output/ex5_5_ncg.out ex5_5_ncg.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5_5_ncg, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_5_ncg.tmp

runex5_5_nrichardson:
	-@${CSD_BASIC_COMMAND_LINE} -snes_type nrichardson \
        > ex5_5_nrichardson.tmp 2>&1; \
	   if (${DIFF} output/ex5_5_nrichardson.out ex5_5_nrichardson.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5_5_nrichardson, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_5_nrichardson.tmp

runex5_5_ngmres_ngs:
	-@${CSD_BASIC_COMMAND_LINE} -snes_type ngmres -npc_snes_type gs -npc_snes_max_it 1 \
           > ex5_5_ngmres_ngs.tmp 2>&1; \
	   if (${DIFF} output/ex5_5_ngmres_ngs.out ex5_5_ngmres_ngs.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5_5_ngmres_ngs, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_5_ngmres_ngs.tmp

runex5_5_qn:
	-@${CSD_BASIC_COMMAND_LINE} -snes_type qn -snes_linesearch_type cp -snes_qn_m ${N_RESTART} \
        > ex5_5_qn.tmp 2>&1; \
	   if (${DIFF} output/ex5_5_qn.out ex5_5_qn.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5_5_qn, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_5_qn.tmp

runex5_5_ls:
	-@${CSD_BASIC_COMMAND_LINE} -snes_type ls \
        > ex5_5_ls.tmp 2>&1; \
	   if (${DIFF} output/ex5_5_ls.out ex5_5_ls.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5_5_ls, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_5_ls.tmp

runex5_5_fas:
	-@${MPIEXEC} -n 1 ./ex5 -fas_coarse_snes_max_it 1 -fas_coarse_pc_type lu -fas_coarse_ksp_type preonly -snes_rtol 1.e-12 -snes_monitor_short -snes_type fas -fas_coarse_ksp_type richardson -da_refine 6 > ex5_5_fas.tmp 2>&1; \
	   if (${DIFF} output/ex5_5_fas.out ex5_5_fas.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5_5_fas, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_5_fas.tmp

runex5_5_fas_additive:
	-@${MPIEXEC} -n 1 ./ex5 -fas_coarse_snes_max_it 1 -fas_coarse_pc_type lu -fas_coarse_ksp_type preonly -snes_rtol 1.e-12 -snes_monitor_short -snes_type fas -fas_coarse_ksp_type richardson -da_refine 6 -snes_fas_type additive -snes_max_it 50 > ex5_5_fas_additive.tmp 2>&1; \
	   if (${DIFF} output/ex5_5_fas_additive.out ex5_5_fas_additive.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5_5_fas_additive, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_5_fas_additive.tmp

runex5_5_ngmres_fas:
	-@${MPIEXEC} -n 1 ./ex5 -snes_type ngmres -npc_fas_coarse_snes_max_it 1 -npc_fas_coarse_snes_type ls -npc_fas_coarse_pc_type lu -npc_fas_coarse_ksp_type preonly -snes_ngmres_m 10 -snes_rtol 1.e-12 -snes_monitor_short -npc_snes_max_it 1 -npc_snes_type fas -npc_fas_coarse_ksp_type richardson -da_refine 6 > ex5_5_ngmres_fas.tmp 2>&1; \
	   if (${DIFF} output/ex5_5_ngmres_fas.out ex5_5_ngmres_fas.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5_5_ngmres_fas, diffs above \n========================================="; fi; \
	   ${RM} -f ex5_5_ngmres_fas.tmp

runex5_6:
	-@${MPIEXEC} -n 4 ./ex5 -snes_converged_reason -ksp_converged_reason -da_grid_x 129 -da_grid_y 129  -pc_type mg -pc_mg_levels 8 -mg_levels_ksp_type chebyshev -mg_levels_ksp_chebyshev_estimate_eigenvalues 0,0.5,0,1.1 -mg_levels_ksp_max_it 2 > ex5_6.tmp 2>&1; \
	   ${DIFF} output/ex5_6.out ex5_6.tmp || echo  ${PWD} "\nPossible problem with ex5_6, diffs above \n========================================="; \
	   ${RM} -f ex5_6.tmp

runex5f:
	-@${MPIEXEC} -n 4 ./ex5f -snes_mf -da_processors_x 4 -da_processors_y 1 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5f_1.tmp 2>&1; \
	   if (${DIFF} output/ex5f_1.out ex5f_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5f_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex5f_1.tmp
runex5f_2:
	-@${MPIEXEC} -n 4 ./ex5f  -da_processors_x 2 -da_processors_y 2 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5f_2.tmp 2>&1; \
	   if (${DIFF} output/ex5f_2.out ex5f_2.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5f_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex5f_2.tmp
runex5f_3:
	-@${MPIEXEC} -n 3 ./ex5f -snes_fd  -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5f_3.tmp 2>&1;\
	   if (${DIFF} output/ex5f_3.out ex5f_3.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5f_3, diffs above \n========================================="; fi; \
	   ${RM} -f ex5f_3.tmp
runex5f_4:
	-@${MPIEXEC} -n 2 ./ex5f -adifor_jacobian -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5f_4.tmp 2>&1;\
	   if (${DIFF} output/ex5f_4.out ex5f_4.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5f_4, diffs above \n========================================="; fi; \
	   ${RM} -f ex5f_4.tmp
runex5f_5:
	-@${MPIEXEC} -n 2 ./ex5f -adiformf_jacobian  -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5f_5.tmp 2>&1;\
	   if (${DIFF} output/ex5f_5.out ex5f_5.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5f_5, diffs above \n========================================="; fi; \
	   ${RM} -f ex5f_5.tmp
testex5f: ex5f.PETSc
	@if [ "${PETSC_WITH_BATCH}" != "" ]; then \
           echo "Running with batch filesystem; to test run src/snes/examples/tutorials/ex5f with" ; \
           echo "your systems batch system"; \
        elif [ "${MPIEXEC}" = "/bin/false" ]; then \
           echo "*mpiexec not found*. Please run src/snes/examples/tutorials/ex5f manually"; \
        elif [ -f ex5f ]; then \
	   ${MPIEXEC} -n 1 ./ex5f > ex5f_1.tmp 2>&1; \
	   if (${DIFF} output/ex5f_1.testout ex5f_1.tmp > /dev/null 2>&1) then \
           echo "Fortran example src/snes/examples/tutorials/ex5f run successfully with 1 MPI process"; \
	   else echo "Possible error running Fortran example src/snes/examples/tutorials/ex5f with 1 MPI process"; \
           echo "See http://www.mcs.anl.gov/petsc/documentation/faq.html";\
           cat ex5f_1.tmp; fi;  \
         ${RM} -f ex5f_1.tmp ;\
         ${MAKE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} ex5f.rm; fi
runex5f90:
	-@${MPIEXEC} -n 4 ./ex5f90 -snes_mf -da_processors_x 4 -da_processors_y 1 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5f90_1.tmp 2>&1; \
	   if (${DIFF} output/ex5f90_1.out ex5f90_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5f90_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex5f90_1.tmp
runex5f90_2:
	-@${MPIEXEC} -n 4 ./ex5f90 -da_processors_x 2 -da_processors_y 2 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5f90_2.tmp 2>&1; \
	   if (${DIFF} output/ex5f90_2.out ex5f90_2.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5f90_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex5f90_2.tmp
runex5f90_3:
	-@${MPIEXEC} -n 3 ./ex5f90 -snes_fd  -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5f90_3.tmp 2>&1;\
	   if (${DIFF} output/ex5f90_3.out ex5f90_3.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5f90_3, diffs above \n========================================="; fi; \
	   ${RM} -f ex5f90_3.tmp
runex5f90_4:
	-@${MPIEXEC} -n 3 ./ex5f90 -snes_mf_operator  -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5f90_4.tmp 2>&1;\
	   if (${DIFF} output/ex5f90_4.out ex5f90_4.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5f90_4, diffs above \n========================================="; fi; \
	   ${RM} -f ex5f90_4.tmp
runex5f90_5:
	-@${MPIEXEC} -n 1 ./ex5f90  > ex5f90_5.tmp 2>&1;\
	   if (${DIFF} output/ex5f90_5.out ex5f90_5.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5f90_5, diffs above \n========================================="; fi; \
	   ${RM} -f ex5f90_5.tmp
runex5f90t:
	-@${MPIEXEC} -n 4 ./ex5f90t -snes_mf -da_processors_x 4 -da_processors_y 1 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex5f90t_1.tmp 2>&1; \
	   if (${DIFF} output/ex5f90_1.out ex5f90t_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex5f90t_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex5f90t_1.tmp
runex6:
	-@${MPIEXEC} -n 1 ./ex6 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex6_1.tmp 2>&1;   \
	   if (${DIFF} output/ex6_1.out ex6_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex6_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex6_1.tmp
runex6_2:
	-@${MPIEXEC} -n 1 ./ex6 -snes_mf -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex6_2.tmp 2>&1;   \
	   if (${DIFF} output/ex6_2.out ex6_2.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex6_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex6_2.tmp
runex8:
	-@${MPIEXEC} -n 2 ./ex8  -mx 16 -my 16 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always \
                    -lidvelocity 0. -grashof 1000. > ex8_1.tmp 2>&1; \
           if (${DIFF} output/ex8_1.out ex8_1.tmp) then true; \
           else  echo  ${PWD} "\nPossible problem with with ex8_1, diffs above \n========================================="; fi; \
           ${RM} -f ex8_1.tmp
runex8_2:
	-@${MPIEXEC} -n 2 ./ex8 -snes_monitor_solution -ksp_gmres_cgs_refinement_type refine_always > /dev/null 2>&1 
runex8_3:
	-@${MPIEXEC} -n 2 ./ex8 -is_coloring_view > ex8_3.tmp 2>&1;\
           if (${DIFF} output/ex8_3.out ex8_3.tmp) then true; \
           else  echo  ${PWD} "\nPossible problem with with ex8_3, diffs above \n========================================="; fi; \
           ${RM} -f ex8_3.tmp
runex14:
	-@${MPIEXEC} -n 4 ./ex14 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex14_1.tmp 2>&1; \
	   if (${DIFF} output/ex14_1.out ex14_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex14_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex14_1.tmp
runex14_2:
	-@${MPIEXEC} -n 4 ./ex14 -snes_monitor_short -ksp_gmres_cgs_refinement_type refine_always  > ex14_2.tmp 2>&1; \
	   if (${DIFF} output/ex14_2.out ex14_2.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex14_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex14_2.tmp
runex15:
	-@${MPIEXEC} -n 2 ./ex15 -snes_monitor_short -da_grid_x 20 -da_grid_y 20 -p 1.3 -lambda 1 -jtype 4 > ex15_1.tmp 2>&1; \
	   ${DIFF} output/ex15_1.out ex15_1.tmp || echo  ${PWD} "\nPossible problem with with ex15_1, diffs above \n========================================="; \
	   ${RM} -f ex15_1.tmp
runex15_2:
	-@${MPIEXEC} -n 2 ./ex15 -snes_monitor_short -da_grid_x 20 -da_grid_y 20 -p 1.3 -lambda 1 -jtype 2 -precheck 1 > ex15_2.tmp 2>&1; \
	   ${DIFF} output/ex15_2.out ex15_2.tmp || echo  ${PWD} "\nPossible problem with with ex15_2, diffs above \n========================================="; \
	   ${RM} -f ex15_2.tmp
runex15_3:
	-@${MPIEXEC} -n 2 ./ex15 -snes_monitor_short -da_grid_x 20 -da_grid_y 20 -p 1.3 -lambda 1 -jtype 2 -picard -precheck 1 > ex15_3.tmp 2>&1; \
	   ${DIFF} output/ex15_3.out ex15_3.tmp || echo  ${PWD} "\nPossible problem with with ex15_3, diffs above \n========================================="; \
	   ${RM} -f ex15_3.tmp
runex19:
	-@${MPIEXEC} -n 2 ./ex19 -da_refine 3 -snes_monitor_short -pc_type mg -ksp_type fgmres -pc_mg_type full > ex19_1.tmp 2>&1; \
	   if (${DIFF} output/ex19_1.out ex19_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_1.tmp
runex19_3: #test pc_redundant
	-@${MPIEXEC} -n 4 ./ex19 -da_refine 3 -snes_monitor_short -pc_type redundant -mat_type mpiaij -redundant_ksp_type preonly -redundant_pc_factor_mat_solver_package mumps -pc_redundant_number 2 > ex19_3.tmp 2>&1; \
	   if (${DIFF} output/ex19_3.out ex19_3.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_3, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_3.tmp
runex19_4: #test pc_redundant
	-@${MPIEXEC} -n 12 ./ex19 -da_refine 3 -snes_monitor_short -pc_type redundant -mat_type mpiaij -redundant_ksp_type preonly -redundant_pc_factor_mat_solver_package mumps -pc_redundant_number 5 > ex19_4.tmp 2>&1; \
	   if (${DIFF} output/ex19_3.out ex19_4.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_4, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_4.tmp
runex19_5: #test different scatters
	-@for A in " " -vecscatter_rsend -vecscatter_ssend -vecscatter_alltoall "-vecscatter_alltoall -vecscatter_nopack" -vecscatter_window; do \
           for B in " " -vecscatter_merge ; do \
             ${MPIEXEC} -n 4 ./ex19 -da_refine 3 -ksp_type fgmres -pc_type mg -pc_mg_type full $$A $$B > ex19_5.tmp 2>&1; \
	     if (${DIFF} output/ex19_5.out ex19_5.tmp) then true; \
	     else  echo  ${PWD} "\nPossible problem with with ex19_5 " $$A $$B " diffs above \n========================================="; fi; \
           done;\
         done;\
	 ${RM} -f ex19_5.tmp
# fieldsplit preconditioner tests
runex19_6: 
	-@${MPIEXEC} -n 1 ./ex19 -snes_monitor_short -ksp_monitor_short -pc_type fieldsplit -snes_view -ksp_type fgmres -da_refine 1  > ex19_6.tmp 2>&1; \
	   if (${DIFF} output/ex19_6.out ex19_6.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_6, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_6.tmp
runex19_7: 
	-@${MPIEXEC} -n 3 ./ex19 -snes_monitor_short -ksp_monitor_short -pc_type fieldsplit -snes_view -da_refine 1 -ksp_type fgmres  > ex19_7.tmp 2>&1; \
	   if (${DIFF} output/ex19_7.out ex19_7.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_7, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_7.tmp
runex19_8: 
	-@${MPIEXEC} -n 1 ./ex19 -snes_monitor_short -ksp_monitor_short -pc_type fieldsplit -pc_fieldsplit_block_size 2 -pc_fieldsplit_0_fields 0,1 -pc_fieldsplit_1_fields 0,1 -pc_fieldsplit_type multiplicative -snes_view   -fieldsplit_pc_type lu -da_refine 1 -ksp_type fgmres > ex19_8.tmp 2>&1; \
	   if (${DIFF} output/ex19_8.out ex19_8.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_8, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_8.tmp
runex19_9: 
	-@${MPIEXEC} -n 3 ./ex19 -snes_monitor_short -ksp_monitor_short -pc_type fieldsplit -pc_fieldsplit_type multiplicative -snes_view -da_refine 1 -ksp_type fgmres  > ex19_9.tmp 2>&1; \
	   if (${DIFF} output/ex19_9.out ex19_9.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_9, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_9.tmp
runex19_10: 
	-@${MPIEXEC} -n 3 ./ex19 -snes_monitor_short -ksp_monitor_short -pc_type fieldsplit -pc_fieldsplit_type symmetric_multiplicative -snes_view  -da_refine 1 -ksp_type fgmres > ex19_10.tmp 2>&1; \
	   if (${DIFF} output/ex19_10.out ex19_10.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_10, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_10.tmp
runex19_11: #test pc_redundant
	-@${MPIEXEC} -n 4 ./ex19  -snes_monitor_short -pc_type redundant -mat_type mpiaij -redundant_pc_factor_mat_solver_package pastix -pc_redundant_number 2 -da_refine 4 -ksp_type fgmres > ex19_11.tmp 2>&1; \
	   if (${DIFF} output/ex19_11.out ex19_11.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_11, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_11.tmp
runex19_12: #test pc_redundant
	-@${MPIEXEC} -n 12 ./ex19  -snes_monitor_short -pc_type redundant -mat_type mpiaij -redundant_pc_factor_mat_solver_package pastix -pc_redundant_number 5  -da_refine 4 -ksp_type fgmres  > ex19_12.tmp 2>&1; \
	   if (${DIFF} output/ex19_12.out ex19_12.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_12, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_12.tmp
runex19_13: #test pc_fieldsplit with -snes_mf_operator
	-@${MPIEXEC} -n 3 ./ex19 -snes_monitor_short -ksp_monitor_short -pc_type fieldsplit -pc_fieldsplit_type multiplicative -snes_view  -da_refine 1 -ksp_type fgmres  -snes_mf_operator > ex19_13.tmp 2>&1; \
	   if (${DIFF} output/ex19_13.out ex19_13.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_13, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_13.tmp
runex19_14:
	-@${MPIEXEC} -n 2 ./ex19 -snes_monitor_short -pc_type mg -dm_mat_type baij -mg_coarse_pc_type bjacobi -da_refine 3 -ksp_type fgmres > ex19_14.tmp 2>&1; \
	   if (${DIFF} output/ex19_14.out ex19_14.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_14, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_14.tmp
runex19_superlu:
	-@${MPIEXEC} -n 1 ./ex19 -da_grid_x 20 -da_grid_y 20 -pc_type lu -pc_factor_mat_solver_package superlu > ex19.tmp 2>&1; \
	   if (${DIFF} output/ex19_superlu.out ex19.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_superlu, diffs above \n========================================="; fi; \
	   ${RM} -f ex19.tmp
runex19_superlu_equil: #This test fails - I'll check it. Hong
	-@${MPIEXEC} -n 1 ./ex19  -da_grid_x 20 -da_grid_y 20 -{snes,ksp}_monitor_short -pc_type lu -pc_factor_mat_solver_package superlu -mat_superlu_equil > ex19.tmp 2>&1; \
	   if (${DIFF} output/ex19_superlu_equil.out ex19.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_superlu_equil, diffs above \n========================================="; fi; \
	   ${RM} -f ex19.tmp
runex19_superlu_dist:
	-@${MPIEXEC} -n 1 ./ex19 -da_grid_x 20 -da_grid_y 20 -pc_type lu -pc_factor_mat_solver_package superlu_dist > ex19.tmp 2>&1; \
	   if (${DIFF} output/ex19_superlu.out ex19.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_superlu_dist, diffs above \n========================================="; fi; \
	   ${RM} -f ex19.tmp
runex19_superlu_dist_2:
	-@${MPIEXEC} -n 2 ./ex19  -da_grid_x 20 -da_grid_y 20 -pc_type lu -pc_factor_mat_solver_package superlu_dist > ex19.tmp 2>&1; \
	   if (${DIFF} output/ex19_superlu.out ex19.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_superlu_dist, diffs above \n========================================="; fi; \
	   ${RM} -f ex19.tmp
runex19_fieldsplit_2:
	-@${MPIEXEC} -n 1 ./ex19 -ksp_type fgmres -pc_type fieldsplit -pc_fieldsplit_block_size 4 -pc_fieldsplit_type additive -pc_fieldsplit_0_fields 0,1,2 -pc_fieldsplit_1_fields 3 -snes_monitor_short -ksp_monitor_short  > ex19_6.tmp 2>&1; \
	   if (${DIFF} output/ex19_fieldsplit_2.out ex19_6.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_fieldsplit_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_6.tmp
runex19_fieldsplit_3:
	-@${MPIEXEC} -n 1 ./ex19 -ksp_type fgmres -pc_type fieldsplit -pc_fieldsplit_block_size 4 -pc_fieldsplit_type additive -pc_fieldsplit_0_fields 0,1,2 -pc_fieldsplit_1_fields 3 -fieldsplit_0_pc_type lu -fieldsplit_1_pc_type lu -snes_monitor_short -ksp_monitor_short  > ex19_6.tmp 2>&1; \
	   if (${DIFF} output/ex19_fieldsplit_3.out ex19_6.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_fieldsplit_3, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_6.tmp
runex19_fieldsplit_4:
	-@${MPIEXEC} -n 1 ./ex19 -ksp_type fgmres -pc_type fieldsplit -pc_fieldsplit_block_size 4 -pc_fieldsplit_type SCHUR -pc_fieldsplit_0_fields 0,1,2 -pc_fieldsplit_1_fields 3 -fieldsplit_0_pc_type lu -fieldsplit_1_pc_type lu -snes_monitor_short -ksp_monitor_short  > ex19_6.tmp 2>&1; \
	   if (${DIFF} output/ex19_fieldsplit_4.out ex19_6.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_fieldsplit_4, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_6.tmp
runex19_fieldsplit_mumps:
	-@${MPIEXEC} -n 2 ./ex19 -pc_type fieldsplit -pc_fieldsplit_block_size 4 -pc_fieldsplit_type SCHUR -pc_fieldsplit_0_fields 0,1,2 -pc_fieldsplit_1_fields 3 -fieldsplit_0_pc_type lu -fieldsplit_1_pc_type lu -snes_monitor_short -ksp_monitor_short  -fieldsplit_0_pc_factor_mat_solver_package mumps -fieldsplit_1_pc_factor_mat_solver_package mumps > ex19_6.tmp 2>&1; \
	   if (${DIFF} output/ex19_fieldsplit_5.out ex19_6.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_fieldsplit_fieldsplit_mumps, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_6.tmp
runex19_fieldsplit_hypre:
	-@${MPIEXEC} -n 2 ./ex19  -pc_type fieldsplit -pc_fieldsplit_block_size 4 -pc_fieldsplit_type SCHUR -pc_fieldsplit_0_fields 0,1,2 -pc_fieldsplit_1_fields 3 -fieldsplit_0_pc_type lu -fieldsplit_0_pc_factor_mat_solver_package mumps -fieldsplit_1_pc_type hypre -fieldsplit_1_pc_hypre_type boomeramg -snes_monitor_short -ksp_monitor_short  > ex19_6.tmp 2>&1; \
	   if (${DIFF} output/ex19_fieldsplit_hypre.out ex19_6.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_fieldsplit_hypre, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_6.tmp
runex19_composite_fieldsplit: #similar to runex19_fieldsplit_2
	-@${MPIEXEC} -n 1 ./ex19  -ksp_type fgmres -pc_type composite -pc_composite_type MULTIPLICATIVE -pc_composite_pcs fieldsplit,none -sub_0_pc_fieldsplit_block_size 4 -sub_0_pc_fieldsplit_type additive -sub_0_pc_fieldsplit_0_fields 0,1,2 -sub_0_pc_fieldsplit_1_fields 3 -snes_monitor_short -ksp_monitor_short  > ex19_6.tmp 2>&1; \
	   if (${DIFF} output/ex19_composite_fieldsplit.out ex19_6.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_composite_fieldsplit, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_6.tmp
runex19_composite_fieldsplit_bjacobi:
	-@${MPIEXEC} -n 1 ./ex19 -ksp_type fgmres -pc_type composite -pc_composite_type MULTIPLICATIVE -pc_composite_pcs fieldsplit,bjacobi -sub_0_pc_fieldsplit_block_size 4 -sub_0_pc_fieldsplit_type additive -sub_0_pc_fieldsplit_0_fields 0,1,2 -sub_0_pc_fieldsplit_1_fields 3 -sub_1_pc_bjacobi_blocks 16 -sub_1_sub_pc_type lu -snes_monitor_short -ksp_monitor_short  > ex19_6.tmp 2>&1; \
	   if (${DIFF} output/ex19_composite_fieldsplit_bjacobi.out ex19_6.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_composite_fieldsplit_bjacobi, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_6.tmp
runex19_composite_fieldsplit_bjacobi_2:
	-@${MPIEXEC} -n 4 ./ex19 -ksp_type fgmres -pc_type composite -pc_composite_type MULTIPLICATIVE -pc_composite_pcs fieldsplit,bjacobi -sub_0_pc_fieldsplit_block_size 4 -sub_0_pc_fieldsplit_type additive -sub_0_pc_fieldsplit_0_fields 0,1,2 -sub_0_pc_fieldsplit_1_fields 3 -sub_1_pc_bjacobi_blocks 16 -sub_1_sub_pc_type lu -snes_monitor_short -ksp_monitor_short  > ex19_6.tmp 2>&1; \
	   if (${DIFF} output/ex19_composite_fieldsplit_bjacobi_2.out ex19_6.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_composite_fieldsplit_bjacobi_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_6.tmp

runex19_fas: #test ex19 with FAS and pointwise GS smoother
	-@${MPIEXEC} -n 1 ./ex19 -da_refine 4 -snes_monitor_short -snes_type fas \
        -fas_levels_snes_type gs -grashof 4e4 -snes_fas_smoothup 6 -snes_fas_smoothdown 6 -lidvelocity 100 > ex19_fas.tmp 2>&1; \
	   ${DIFF} output/ex19_fas.out ex19_fas.tmp || echo  ${PWD} "\nPossible problem with ex19_fas, diffs above \n========================================="; \
           ${RM} -f ex19_fas.tmp

runex19_ngmres_fas: #test ex19 with NGMRES preconditioned by FAS with pointwise GS smoother
	-@${MPIEXEC} -n 1 ./ex19 -da_refine 4 -snes_monitor_short -snes_type ngmres \
        -npc_snes_type fas -npc_fas_levels_snes_type gs -npc_snes_max_it 1 -npc_snes_fas_smoothup 6 -npc_snes_fas_smoothdown 6  \
        -lidvelocity 100 -grashof 4e4 \
         > ex19_ngmres_fas.tmp 2>&1; \
	   ${DIFF} output/ex19_ngmres_fas.out ex19_ngmres_fas.tmp || echo  ${PWD} "\nPossible problem with ex19_ngmres_fas, diffs above \n========================================="; \
           ${RM} -f ex19_ngmres_fas.tmp

runex19_ngmres_fas_ms: #test ex19 with NGMRES preconditioned by FAS with multi-stage smoother
	-@${MPIEXEC} -n 2 ./ex19 -snes_grid_sequence 2 -lidvelocity 200 -grashof 1e4 -snes_monitor_short -snes_view -snes_converged_reason -snes_type ngmres -npc_snes_type fas -npc_fas_coarse_snes_type ls -npc_fas_coarse_ksp_type preonly -npc_snes_max_it 1 -npc_fas_snes_type ms -npc_fas_snes_max_it 1 -npc_fas_ksp_type preonly -npc_fas_pc_type none -npc_fas_snes_ms_type m62 -npc_fas_snes_max_it 1 -npc_fas_snes_ms_damping 0.22 > ex19_ngmres_fas_ms.tmp 2>&1; \
	   ${DIFF} output/ex19_ngmres_fas_ms.out ex19_ngmres_fas_ms.tmp || echo  ${PWD} "\nPossible problem with ex19_ngmres_fas_ms, diffs above \n========================================="; \
           ${RM} -f ex19_ngmres_fas_ms.tmp

runex19_bjacobi: #test hierarchical pc
	-@${MPIEXEC} -n 4 ./ex19 -da_refine 4 -ksp_type fgmres -pc_type bjacobi -pc_bjacobi_blocks 2 -sub_ksp_type gmres -sub_ksp_max_it 2 -sub_pc_type bjacobi -sub_sub_ksp_type preonly -sub_sub_pc_type ilu -snes_monitor_short > ex19_1.tmp 2>&1; \
	   if (${DIFF} output/ex19_bjacobi.out ex19_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_bjacobi, diffs above \n========================================="; fi; \
	   ${RM} -f ex19_1.tmp

runex18:
	-@${MPIEXEC} -n 1 ./ex18 -pc_type mg -ksp_type fgmres -da_refine 2 -pc_mg_galerkin -snes_view > ex18_1.tmp 2>&1; \
	   if (${DIFF} output/ex18_1.out ex18_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex18, diffs above \n========================================="; fi; \
	   ${RM} -f ex18_1.tmp


testex19: ex19.PETSc
	-@if [ "${PETSC_WITH_BATCH}" != "" ]; then \
           echo "Running with batch filesystem; to test run src/snes/examples/tutorials/ex19 with" ; \
           echo "your systems batch system"; \
        elif [ "${MPIEXEC}" = "/bin/false" ]; then \
           echo "*mpiexec not found*. Please run src/snes/examples/tutorials/ex19 manually"; \
	elif [ -f ex19 ]; then \
           ${MPIEXEC} -n 1 ./ex19 -da_refine 3 -pc_type mg -ksp_type fgmres  > ex19_1.tmp 2>&1; \
	   if (${DIFF} output/ex19_1.testout ex19_1.tmp > /dev/null 2>&1) then \
           echo "C/C++ example src/snes/examples/tutorials/ex19 run successfully with 1 MPI process"; \
	   else echo "Possible error running C/C++ src/snes/examples/tutorials/ex19 with 1 MPI process"; \
           echo "See http://www.mcs.anl.gov/petsc/documentation/faq.html";\
           cat ex19_1.tmp; fi; \
	if [ "${MPIEXEC}" != "${PETSC_DIR}/bin/mpiexec.uni" ]; then \
           ${MPIEXEC} -n 2 ./ex19 -da_refine 3 -pc_type mg -ksp_type fgmres  > ex19_1.tmp 2>&1; \
	   if (${DIFF} output/ex19_1.testout ex19_1.tmp > /dev/null 2>&1) then \
           echo "C/C++ example src/snes/examples/tutorials/ex19 run successfully with 2 MPI processes"; \
	   else echo "Possible error running C/C++ src/snes/examples/tutorials/ex19 with 2 MPI processes"; \
           echo "See http://www.mcs.anl.gov/petsc/documentation/faq.html";\
           cat ex19_1.tmp; fi; fi; \
        ${RM} -f ex19_1.tmp; \
        ${MAKE} PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} ex19.rm ; fi

runex19_cusp:
	-@${MPIEXEC} -n 1 ./ex19 -dm_vec_type cusp -dm_mat_type aijcusp -pc_type none -ksp_type fgmres -snes_monitor_short -snes_rtol 1.e-5 > ex19.tmp 2>&1; \
	   if (${DIFF} output/ex19_cusp.out ex19.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex19_cusp, diffs above \n========================================="; fi; \
	   ${RM} -f ex19.tmp
runex20:
	-@${MPIEXEC} -n 4 ./ex20 -snes_view -snes_monitor_short -pc_mg_type full -ksp_type fgmres -pc_type mg -snes_view -pc_mg_levels 2 -pc_mg_galerkin  > ex20_1.tmp 2>&1; \
	   if (${DIFF} output/ex20_1.out ex20_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex20_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex20_1.tmp
runex21:
	-@${MPIEXEC} -n 4 ./ex21 -snes_linesearch_monitor -snes_mf -snes_monitor_short -nox -ksp_monitor -snes_converged_reason > ex21_1.tmp 2>&1; \
	   if (${DIFF} output/ex21_1.out ex21_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex21_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex21_1.tmp
runex22:
	-@${MPIEXEC} -n 2 ./ex22 -da_grid_x 10 -snes_monitor_short -ksp_monitor_short -snes_view > ex22_1.tmp 2>&1; \
	   if (${DIFF} output/ex22_1.out ex22_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex22_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex22_1.tmp

runex25:
	-@${MPIEXEC} -n 1 ./ex25 -pc_type mg -da_refine 1  -ksp_type fgmres  > ex25_1.tmp 2>&1;	  \
	   if (${DIFF} output/ex25_1.out ex25_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex25_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex25_1.tmp

runex25_2:
	-@${MPIEXEC} -n 2 ./ex25 -pc_type mg -da_refine 1  -ksp_type fgmres > ex25_2.tmp 2>&1;	  \
	   if (${DIFF} output/ex25_2.out ex25_2.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex25_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex25_2.tmp

runex28_0:
	-@${MPIEXEC} -n 3 ./ex28 -da_grid_x 20 -snes_converged_reason -snes_monitor_short -problem_type 0 > ex28_0.tmp 2>&1; \
	  ${DIFF} output/ex28_0.out ex28_0.tmp || echo  ${PWD} "\nPossible problem with ex28_0, diffs above \n========================================="; ${RM} -f ex28_0.tmp
runex28_1:
	-@${MPIEXEC} -n 3 ./ex28 -da_grid_x 20 -snes_converged_reason -snes_monitor_short -problem_type 1 > ex28_1.tmp 2>&1; \
	  ${DIFF} output/ex28_1.out ex28_1.tmp || echo  ${PWD} "\nPossible problem with ex28_1, diffs above \n========================================="; ${RM} -f ex28_1.tmp
runex28_2:
	-@${MPIEXEC} -n 3 ./ex28 -da_grid_x 20 -snes_converged_reason -snes_monitor_short -problem_type 2 > ex28_2.tmp 2>&1; \
	  ${DIFF} output/ex28_2.out ex28_2.tmp || echo  ${PWD} "\nPossible problem with ex28_2, diffs above \n========================================="; ${RM} -f ex28_2.tmp
runex28_3:
	-@for mtype in aij nest ; do \
	    ${MPIEXEC} -n 3 ./ex28 -da_grid_x 20 -snes_converged_reason -snes_monitor_short -ksp_monitor_short -problem_type 2 -snes_mf_operator \
	    -pack_dm_mat_type $$mtype -pc_type fieldsplit -pc_fieldsplit_type additive -fieldsplit_u_ksp_type gmres -fieldsplit_k_pc_type jacobi > ex28_3.tmp 2>&1; \
	    ${DIFF} output/ex28_3.out ex28_3.tmp || echo  ${PWD} "\nPossible problem with ex28_3 mtype=$${mtype}; diffs above \n========================================="; ${RM} -f ex28_3.tmp; \
	  done
runex28_4:
	-@${MPIEXEC} -n 6 ./ex28 -da_grid_x 257 -snes_converged_reason -snes_monitor_short -ksp_monitor_short -problem_type 2 -snes_mf_operator -pack_dm_mat_type aij -pc_type fieldsplit -pc_fieldsplit_type multiplicative -fieldsplit_u_ksp_type gmres -fieldsplit_u_ksp_pc_side right -fieldsplit_u_pc_type mg -fieldsplit_u_pc_mg_levels 4 -fieldsplit_u_mg_levels_ksp_type richardson -fieldsplit_u_mg_levels_ksp_max_it 1 -fieldsplit_u_mg_levels_pc_type sor -fieldsplit_u_pc_mg_galerkin -fieldsplit_u_ksp_converged_reason -fieldsplit_k_pc_type jacobi > ex28_4.tmp 2>&1; \
	  ${DIFF} output/ex28_4.out ex28_4.tmp || echo  ${PWD} "\nPossible problem with ex28_4 ; diffs above \n========================================="; ${RM} -f ex28_4.tmp;
runex30:
	-@${MPIEXEC} -n 1 ./ex30  > ex30_1.tmp 2>&1;	  \
	   if (${DIFF} output/ex30_1.out ex30_1.tmp) then true; \
	   else  echo  ${PWD} "\nPossible problem with with ex30_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex30_1.tmp

runex40f90:
	-@${MPIEXEC} -n 1 ./ex40f90 -snes_monitor_short -snes_view  -da_refine 1 -pc_type mg  -pc_mg_type full -ksp_type fgmres -pc_mg_galerkin  > ex40f90_1.tmp 2>&1;	  \
	   if (${DIFF} output/ex40f90.out ex40f90_1.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex40f90_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex40f90_1.tmp

runex42:
	-@${MPIEXEC} -n 1 ./ex42 -snes_monitor_short -snes_max_it 1000 -snes_rtol 1.e-14  > ex42_1.tmp 2>&1;	  \
	   if (${DIFF} output/ex42_1.out ex42_1.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex42_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex42_1.tmp

runex44:
	-@${MPIEXEC} -n 1 ./ex44 -n 10 -snes_mf -user_precond > ex44.tmp 2>&1; \
	   if (${DIFF} output/ex44_0.out ex44.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex44, diffs above \n========================================="; fi; \
	   ${RM} -f ex44.tmp
runex44_2:
	-@${MPIEXEC} -n 1 ./ex44 -n 10 -snes_mf > ex44.tmp 2>&1; \
	   if (${DIFF} output/ex44_0.out ex44.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex44_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex44.tmp

runex45:
	-@${MPIEXEC} -n 1 ./ex45 -n 10 > ex45.tmp 2>&1; \
	   if (${DIFF} output/ex45_0.out ex45.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex45, diffs above \n========================================="; fi; \
	   ${RM} -f ex45.tmp
runex45_2:
	-@${MPIEXEC} -n 1 ./ex45 -n 10 -snes_jacobian_default -pc_type jacobi > ex45.tmp 2>&1; \
	   if (${DIFF} output/ex45_0.out ex45.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex45_2, diffs above \n========================================="; fi; \
	   ${RM} -f ex45.tmp
runex45_3:
	-@${MPIEXEC} -n 1 ./ex45 -n 10 -snes_jacobian_default -fd_jacobian_coloring -my_jacobian_struct > ex45.tmp 2>&1; \
	   if (${DIFF} output/ex45_0.out ex45.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex45_3, diffs above \n========================================="; fi; \
	   ${RM} -f ex45.tmp
runex46:
	-@${MPIEXEC} -n 1 ./ex46 -snes_view -snes_monitor_short -da_refine 1 -pc_type mg -ksp_type fgmres -pc_mg_type full -mg_levels_ksp_chebyshev_estimate_eigenvalues 0.5,1.1  > ex46.tmp 2>&1; \
	   if (${DIFF} output/ex46_1.out ex46.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex46, diffs above \n========================================="; fi; \
	   ${RM} -f ex46.tmp

runex47cu:
	-@${MPIEXEC} -n 1 ./ex47cu -snes_monitor_short -dm_vec_type cusp  > ex47cu_1.tmp 2>&1;	  \
	   if (${DIFF} output/ex47cu_1.out ex47cu_1.tmp) then true; \
	   else echo ${PWD} ; echo "Possible problem with with ex47cu_1, diffs above \n========================================="; fi; \
	   ${RM} -f ex47cu_1.tmp

runex48:
	-@${MPIEXEC} -n 1 ./ex48 -M 6 -P 4 -da_refine 1 -snes_monitor_short -snes_converged_reason -ksp_monitor_short -ksp_converged_reason -thi_mat_type sbaij -ksp_type fgmres -pc_type mg -pc_mg_type full -mg_levels_ksp_type gmres -mg_levels_ksp_max_it 1 -mg_levels_pc_type icc > ex48.tmp 2>&1; \
	   ${DIFF} output/ex48_1.out ex48.tmp || echo  ${PWD} "\nPossible problem with ex48_1, diffs above \n========================================="; \
	   ${RM} -f ex48.tmp
runex48_2:
	-@${MPIEXEC} -n 2 ./ex48 -M 6 -P 4 -thi_hom z -snes_monitor_short -snes_converged_reason -ksp_monitor_short -ksp_converged_reason -thi_mat_type sbaij -ksp_type fgmres -pc_type mg -pc_mg_type full -mg_levels_ksp_type gmres -mg_levels_ksp_max_it 1 -mg_levels_pc_type asm -mg_levels_pc_asm_blocks 6 -mg_levels_0_pc_type redundant -snes_grid_sequence 1 -mat_partitioning_type current > ex48.tmp 2>&1; \
	   ${DIFF} output/ex48_2.out ex48.tmp || echo  ${PWD} "\nPossible problem with ex48_2, diffs above \n========================================="; \
	   ${RM} -f ex48.tmp
runex48_3:
	-@${MPIEXEC} -n 3 ./ex48 -M 7 -P 4 -thi_hom z -da_refine 1 -snes_monitor_short -snes_converged_reason -ksp_monitor_short -ksp_converged_reason -thi_mat_type baij -ksp_type fgmres -pc_type mg -pc_mg_type full -mg_levels_pc_type asm -mg_levels_pc_asm_blocks 9 -mg_levels_ksp_type gmres -mg_levels_ksp_max_it 1 -mat_partitioning_type current > ex48.tmp 2>&1; \
	   ${DIFF} output/ex48_3.out ex48.tmp || echo  ${PWD} "\nPossible problem with ex48_3, diffs above \n========================================="; \
	   ${RM} -f ex48.tmp
runex48_4:
	-@${MPIEXEC} -n 6 ./ex48 -M 4 -P 2 -da_refine_hierarchy_x 1,1,3 -da_refine_hierarchy_y 2,2,1 -da_refine_hierarchy_z 2,2,1 -snes_grid_sequence 3 -ksp_converged_reason -ksp_type fgmres -ksp_rtol 1e-2 -pc_type mg -mg_levels_ksp_type gmres -mg_levels_ksp_max_it 1 -mg_levels_pc_type bjacobi -mg_levels_1_sub_pc_type cholesky -pc_mg_type multiplicative -snes_converged_reason -snes_stol 1e-12 -thi_L 80e3 -thi_alpha 0.05 -thi_friction_m 1 -thi_hom x -snes_view -mg_levels_0_pc_type redundant -mg_levels_0_ksp_type preonly > ex48.tmp 2>&1; \
	   ${DIFF} output/ex48_4.out ex48.tmp || echo  ${PWD} "\nPossible problem with ex48_4, diffs above \n========================================="; \
	   ${RM} -f ex48.tmp
runex48_5:
	-@for mtype in aij baij sbaij; do \
           ${MPIEXEC} -n 6 ./ex48 -M 12 -P 5 -snes_monitor_short -ksp_converged_reason -pc_type asm -pc_asm_type restrict -dm_mat_type $$mtype > ex48_5.tmp 2>&1; \
	   ${DIFF} output/ex48_5.out ex48_5.tmp || echo  ${PWD} "\nPossible problem with ex48_5 mtype=$${mtype}, diffs above \n========================================="; \
	   ${RM} -f ex48.tmp; \
         done

runex58:
	-@${MPIEXEC} -n 1 ./ex58 -pc_type mg -ksp_monitor_short -pc_mg_galerkin -da_refine 5   -snes_vi_monitor -pc_mg_type full -snes_rtol 1.e-12 -snes_max_it 100 -snes_converged_reason > ex58.tmp 2>&1; \
	   ${DIFF} output/ex58_1.out ex58.tmp || echo  ${PWD} "\nPossible problem with ex58_1, diffs above \n========================================="; \
	   ${RM} -f ex58.tmp
runex58_2:
	-@${MPIEXEC} -n 1 ./ex58 -snes_type viss -pc_type mg -ksp_monitor_short -pc_mg_galerkin -da_refine 5   -snes_vi_monitor -pc_mg_type full -snes_rtol 1.e-12 -snes_max_it 100 -snes_converged_reason > ex58.tmp 2>&1; \
	   ${DIFF} output/ex58_2.out ex58.tmp || echo  ${PWD} "\nPossible problem with ex58_2, diffs above \n========================================="; \
	   ${RM} -f ex58.tmp

TESTEXAMPLES_C		       = ex1.PETSc runex1 ex1.rm ex2.PETSc runex2  runex2_3 ex2.rm ex3.PETSc runex3 \
                                 runex3_2 runex3_3 runex3_4 ex3.rm ex4.PETSc runex4 ex4.rm ex5.PETSc runex5 runex5_2 runex5_3 runex5_4 \
                                 runex5_5_ngmres runex5_5_ngmres_nrichardson runex5_5_ncg runex5_5_nrichardson \
                                 runex5_5_ngmres_ngs runex5_5_qn runex5_5_ls \
                                 runex5_5_fas runex5_5_ngmres_fas runex5_5_fas_additive \
                                 runex5_6 ex5.rm\
                                 ex6.PETSc runex6 runex6_2 ex6.rm\
                                 ex14.PETSc runex14 runex14_2 ex14.rm ex15.PETSc runex15 runex15_3 ex15.rm ex18.PETSc runex18 ex18.rm \
                                 ex19.PETSc runex19 runex19_5 \
                                 runex19_6 runex19_fieldsplit_2 runex19_fieldsplit_3 runex19_fieldsplit_4 \
                                 runex19_composite_fieldsplit runex19_composite_fieldsplit_bjacobi runex19_composite_fieldsplit_bjacobi_2\
                                 runex19_7 runex19_8 runex19_9 \
                                 runex19_10 runex19_14   runex19_fas runex19_bjacobi ex19.rm \
                                 ex22.PETSc runex22 ex22.rm ex18.PETSc ex18.rm ex20.PETSc runex20 ex20.rm ex21.PETSc runex21 ex21.rm \
                                 ex25.PETSc runex25 runex25_2 \
                                 ex25.rm ex28.PETSc runex28_0 runex28_1 runex28_2 runex28_3 ex28.rm ex30.PETSc runex30 ex30.rm \
                                 ex42.PETSc runex42 ex42.rm ex43.PETSc ex43.rm \
                                 ex44.PETSc runex44 runex44_2 ex44.rm ex45.PETSc runex45 runex45_2 runex45_3 ex45.rm ex46.PETSc runex46 ex46.rm \
                                 ex48.PETSc runex48 runex48_2 runex48_3 runex48_4 runex48_5 ex48.rm \
                                 ex58.PETSc runex58 runex58_2 ex58.rm \
                                 ex54.PETSc ex54.rm ex59.PETSc ex59.rm
TESTEXAMPLES_C_X	       = ex1.PETSc runex1_X
TESTEXAMPLES_FORTRAN	       = ex1f.PETSc runex1f ex1f.rm ex40f90.PETSc runex40f90 ex40f90.rm
TESTEXAMPLES_C_NOCOMPLEX       = ex30.PETSc ex30.rm
TESTEXAMPLES_FORTRAN_NOCOMPLEX = ex5f.PETSc runex5f runex5f_3 ex5f.rm
TESTEXAMPLES_FORTRAN_MPIUNI    = ex1f.PETSc runex1f ex1f.rm
TESTEXAMPLES_C_X_MPIUNI      = ex1.PETSc runex1 ex1.rm ex2.PETSc runex2 ex2.rm ex3.PETSc runex3 ex3.rm ex6.PETSc \
                                 runex6 runex6_2 ex6.rm  ex14.PETSc ex14.rm 
TESTEXAMPLES_F90	       = ex5f90.PETSc runex5f90 runex5f90_2 runex5f90_3 runex5f90_4 runex5f90_5 ex5f90.rm
TESTEXAMPLES_13		       = 
TESTEXAMPLES_MATLAB_ENGINE     = 
TESTEXAMPLES_AMS	       = 
TESTEXAMPLES_ADIFOR	       = ex5f.PETSc runex5f_2 ex5f.rm
TESTEXAMPLES_MUMPS             = ex19.PETSc runex19_4 runex19_fieldsplit_mumps runex19_3 runex19_4  ex19.rm 
TESTEXAMPLES_SUPERLU           = ex19.PETSc runex19_superlu ex19.rm
TESTEXAMPLES_SUPERLU_DIST      = ex19.PETSc runex19_superlu_dist runex19_superlu_dist_2 ex19.rm  ex19.PETSc  ex19.rm
TESTEXAMPLES_PASTIX            = ex19.PETSc runex19_11 runex19_12 ex19.rm  ex19.PETSc runex19_11 runex19_12 ex19.rm
TESTEXAMPLES_CUSP              = ex19.PETSc runex19_cusp ex19.rm ex47cu.PETSc runex47cu ex47cu.rm 

include ${PETSC_DIR}/conf/test
