using OrdinaryDiffEq
using Trixi
using Plots

###############################################################################
# Semidiscretization of the compressible Euler equations

# Fluid parameters
const gamma = 5 / 3
const prandtl_number = 0.72

# Parameters for compressible von-Karman vortex street
const Re = 500
const Ma = 0.5f0
const D = 1 # Diameter of the cylinder as in the mesh file

# Parameters that can be freely chosen
const v_in = 1
const p_in = 1

# Parameters that follow from Reynolds and Mach number + adiabatic index gamma
const mu = v_in * D / Re

const c = v_in / Ma
const p_over_rho = c^2 / gamma
const rho_in = p_in / p_over_rho

# Equations for this configuration
equations = CompressibleEulerEquations2D(gamma)
equations_parabolic = CompressibleNavierStokesDiffusion2D(equations, mu = mu,
                                                          Prandtl = prandtl_number,
                                                          gradient_variables = GradientVariablesPrimitive())

# Freestream configuration
@inline function initial_condition(x, t, equations::CompressibleEulerEquations2D)
    rho = rho_in
    v1 = v_in
    v2 = 0.0
    p = p_in

    prim = SVector(rho, v1, v2, p)
    return prim2cons(prim, equations)
end

# Symmetric mesh which is refined around the cylinder and in the wake region
mesh_file = joinpath(@__DIR__, "out/sample.inp")
mesh = P4estMesh{2}(mesh_file)

bc_freestream = BoundaryConditionDirichlet(initial_condition)

# Boundary names follow from the mesh file.
# Since this mesh is been generated using the symmetry feature of 
# HOHQMesh.jl (https://trixi-framework.github.io/HOHQMesh.jl/stable/tutorials/symmetric_mesh/)
# the mirrored boundaries are named with a "_R" suffix.
boundary_conditions = Dict(:cylinder => boundary_condition_slip_wall, # top half of the cylinder
                           :cylinder_R => boundary_condition_slip_wall, # bottom half of the cylinder
                           :top => bc_freestream,
                           :top_R => bc_freestream, # aka bottom
                           :outflow => bc_freestream,
                           :outflow_R => bc_freestream,
                           :inflow => bc_freestream,
                           :inflow_R => bc_freestream)

# Parabolic boundary conditions                            
velocity_bc_free = NoSlip((x, t, equations) -> SVector(v_in, 0))
# Use adiabatic also on the boundaries to "copy" temperature from the domain
heat_bc_free = Adiabatic((x, t, equations) -> 0)
boundary_condition_free = BoundaryConditionNavierStokesWall(velocity_bc_free, heat_bc_free)

velocity_bc_cylinder = NoSlip((x, t, equations) -> SVector(0, 0))
heat_bc_cylinder = Adiabatic((x, t, equations) -> 0)
boundary_condition_cylinder = BoundaryConditionNavierStokesWall(velocity_bc_cylinder,
                                                                heat_bc_cylinder)

boundary_conditions_para = Dict(:cylinder => boundary_condition_cylinder, # top half of the cylinder
                                :cylinder_R => boundary_condition_cylinder, # bottom half of the cylinder
                                :top => boundary_condition_free,
                                :top_R => boundary_condition_free, # aka bottom
                                :outflow => boundary_condition_free,
                                :outflow_R => boundary_condition_free,
                                :inflow => boundary_condition_free,
                                :inflow_R => boundary_condition_free)
# Standard DGSEM sufficient here
solver = DGSEM(polydeg = 3, surface_flux = flux_hll)

semi = SemidiscretizationHyperbolicParabolic(mesh, (equations, equations_parabolic),
                                             initial_condition, solver,
                                             boundary_conditions = (boundary_conditions,
                                                                    boundary_conditions_para))

###############################################################################
# Setup an ODE problem
tspan = (0, 100)
ode = semidiscretize(semi, tspan)

# Callbacks
summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)
visualization = VisualizationCallback(interval = 100)
callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution, visualization)

###############################################################################
# run the simulation

time_int_tol = 1e-7
sol = solve(ode,
            # Moderate number of threads (e.g. 4) advisable to speed things up
            RDPK3SpFSAL49(thread = OrdinaryDiffEq.True());
            abstol = time_int_tol, reltol = time_int_tol,
            ode_default_options()..., callback = callbacks)

summary_callback() # print the timer summary