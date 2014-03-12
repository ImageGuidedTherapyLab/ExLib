#include <stdlib.h>

#define NUM_QUADRATURE_POINTS_0 1

/* Quadrature points
   - (x1,y1,x2,y2,...) */
static PetscReal points_0[3] = {
  -0.5,
  -0.5,
  -0.5};

/* Quadrature weights
   - (v1,v2,...) */
static PetscReal weights_0[1] = {1.33333333333};

#define SPATIAL_DIM_0 3

#define NUM_BASIS_FUNCTIONS_0 4

#define NUM_BASIS_COMPONENTS_0 1

/* Number of degrees of freedom for each dimension */
static int numDof_0[4] = {
  1,
  0,
  0,
  0};

/* Nodal basis function evaluations
    - basis function is fastest varying, then point */
static PetscReal Basis_0[4] = {
  0.25,
  0.25,
  0.25,
  0.25};

/* Nodal basis function derivative evaluations,
    - derivative direction fastest varying, then basis function, then point */
static PetscReal BasisDerivatives_0[12] = {
  -0.5,
  -0.5,
  -0.5,
  0.5,
  0.0,
  1.38777878078e-17,
  0.0,
  0.5,
  0.0,
  0.0,
  0.0,
  0.5};
