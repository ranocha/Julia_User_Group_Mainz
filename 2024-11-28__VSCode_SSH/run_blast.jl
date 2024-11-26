using Tenkai
using Plots
plotlyjs()
# Submodules
Eq = Tenkai.EqEuler1D

#------------------------------------------------------------------------------
xmin, xmax = 0.0, 1.0

boundary_condition = (reflect, reflect)
γ = 1.4
final_time = 0.038

initial_value = Eq.blast
exact_solution = Eq.exact_blast # dummy function
boundary_value = Eq.exact_blast # dummy function

degree = 4
solver = "lwfr"
solution_points = "gl"
correction_function = "radau"
numerical_flux = Eq.rusanov
bound_limit = "yes"
bflux = evaluate

nx = 100
cfl = 0.0
bounds = ([-Inf], [Inf]) # Not used in Euler
tvbM = 300.0
save_iter_interval = 0
save_time_interval = 0.1 * final_time
animate = true # Factor on save_iter_interval or save_time_interval
compute_error_interval = 0

# blend parameters
indicator_model = "gassner"
debug_blend = false
cfl_safety_factor = 0.95
pure_fv = false
#------------------------------------------------------------------------------
grid_size = nx
domain = [xmin, xmax]
problem = Problem(domain, initial_value, boundary_value,
                  boundary_condition, final_time, exact_solution)
equation = Eq.get_equation(γ)
limiter = setup_limiter_blend(blend_type = mh_blend(equation),
                              # indicating_variables = Eq.rho_p_indicator!,
                              indicating_variables = Eq.rho_p_indicator!,
                              reconstruction_variables = conservative_reconstruction,
                              indicator_model = indicator_model,
                              debug_blend = debug_blend,
                              pure_fv = pure_fv,
                              numflux = Eq.rusanov)
# limiter = setup_limiter_tvb(equation; tvbM = tvbM)
# limiter = setup_limiter_hierarchical(alpha = 1.0,
#                                      reconstruction = characteristic_reconstruction)
scheme = Scheme(solver, degree, solution_points, correction_function,
                numerical_flux, bound_limit, limiter, bflux)
param = Parameters(grid_size, cfl, bounds, save_iter_interval,
                   save_time_interval, compute_error_interval;
                   animate = animate, cfl_safety_factor = cfl_safety_factor,
                   time_scheme = "SSPRK33")
#------------------------------------------------------------------------------
sol = Tenkai.solve(equation, problem, scheme, param);

sol["plot_data"].p_ua
