# Solving implicit nonlinear equations using (staggered) finite differences, automatic differentiation and sparsity pattern detection

The purpose of these examples is to show how we can use a few julia tools to solve nonlinear equations using staggered finite differences and a few tools. The only thing a user need to provide is a residual routine, from which a sparse matrix is automatically extracted.

We will go through the exercise in the following order 

### exercise_1.jl
Solves a 1D steady state heat diffusion with linear coefficients

### exercise_2.jl
2D diffusion with linear coefficients. Also introduces helper functions to deal with multiple fields

### exercise_3.jl
2 coupled 1D steady state diffusion equations with nonlinear coefficients. 

### exercise_4.jl
1D viscoelastic porosiuty waves for viscous and elastic endmembers. Also shows how to generate an `*.mp4` movie using GLMakie

### exercise_stokes.jl
2D staggered grid Stokes equation using a velocity-pressure formulation for a falling sphere setup. Whereas the example is for a linear viscous rheology, all is setup to deal with nonlinear rheologies as well. 


