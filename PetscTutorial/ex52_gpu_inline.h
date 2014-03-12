#include <stdlib.h>

const int numQuadraturePoints_0 = 1;

/* Quadrature points
   - (x1,y1,x2,y2,...) */
const PetscReal points_0[3] = {
  -0.5,
  -0.5,
  -0.5};

/* Quadrature weights
   - (v1,v2,...) */
const PetscReal weights_0[1] = {1.33333333333};

const int numBasisFunctions_0 = 4;

const int numBasisComponents_0 = 1;

/* Nodal basis function evaluations
    - basis function is fastest varying, then point */
const PetscReal Basis_0[4] = {
  0.25,
  0.25,
  0.25,
  0.25};

/* Nodal basis function derivative evaluations,
    - derivative direction fastest varying, then basis function, then point */
const float3 BasisDerivatives_0[4] = {
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

#define f1_func f1_laplacian

#define f1_coef_func f1_laplacian_coef

/* Number of concurrent blocks */
const int N_bl = 1;
