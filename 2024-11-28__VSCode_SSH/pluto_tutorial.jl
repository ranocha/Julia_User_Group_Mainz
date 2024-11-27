### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ b0f6a076-ace9-11ef-2134-31b8a2d573b5
begin
	import Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	using Plots
end

# ╔═╡ de231021-e1de-438c-bee8-765c06cb9e12
using DelimitedFiles

# ╔═╡ 21cc7f77-23e7-4d10-8682-400154aa8f28
using GZip

# ╔═╡ deb410dd-4427-48da-9d88-1b45d890d9b5
using Tenkai

# ╔═╡ a2bc3d92-2651-473e-b314-40b6ee9287a2
begin
exact_data_raw = joinpath(Tenkai.data_dir, "blast.dat.gz");
exact_data = GZip.open(exact_data_raw);
exact_blast = readdlm(exact_data);
end;

# ╔═╡ 35a4a72d-4803-4d05-87dd-a84c06d86e0a
let
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

nx = 400
cfl = 0.0
bounds = ([-Inf], [Inf]) # Not used in Euler
tvbM = 300.0
save_iter_interval = 0
save_time_interval = 0.0 * final_time
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
                              indicating_variables = Eq.rho_p_indicator!,
                              reconstruction_variables = conservative_reconstruction,
                              indicator_model = indicator_model,
                              debug_blend = debug_blend,
                              pure_fv = pure_fv,
                              numflux = Eq.rusanov)
scheme = Scheme(solver, degree, solution_points, correction_function,
                numerical_flux, bound_limit, limiter, bflux)
param = Parameters(grid_size, cfl, bounds, save_iter_interval,
                   save_time_interval, compute_error_interval;
                   animate = animate, cfl_safety_factor = cfl_safety_factor,
                   time_scheme = "SSPRK33")
#------------------------------------------------------------------------------
sol = Tenkai.solve(equation, problem, scheme, param);

global x, u = sol["grid"].xc, sol["ua"][1,1:end-1]

sol["plot_data"].p_ua

end

# ╔═╡ f22301bb-1c37-4e8b-9964-2ac6e846da7d
plot(x, u, label = "Tenkai.jl")

# ╔═╡ 2e50bc03-76ac-41f5-bc63-ffe3824690d8
let
plotlyjs()
global p = scatter(x, u, label = "Tenkai.jl")
plot!(exact_blast[:,1], exact_blast[:,2], label = "Reference")
end

# ╔═╡ 3e58a59c-5267-43f5-bd6b-90626d022135
let
anim = Animation()
a, b = 230, 380
x_res = x[a:b]
p_ = scatter(x_res, u[a:b], label = "Tenkai.jl")
for i in 1:10
	y = u[a-5*i:b-5*i]
	p_[1][1][:y] .= y # [subplot_index][curve_index][y_series]
	frame(anim, p_)
end
gif(anim, "soln.mp4", fps = 1)
end

# ╔═╡ Cell order:
# ╠═b0f6a076-ace9-11ef-2134-31b8a2d573b5
# ╠═de231021-e1de-438c-bee8-765c06cb9e12
# ╠═21cc7f77-23e7-4d10-8682-400154aa8f28
# ╠═deb410dd-4427-48da-9d88-1b45d890d9b5
# ╠═a2bc3d92-2651-473e-b314-40b6ee9287a2
# ╠═35a4a72d-4803-4d05-87dd-a84c06d86e0a
# ╠═f22301bb-1c37-4e8b-9964-2ac6e846da7d
# ╠═2e50bc03-76ac-41f5-bc63-ffe3824690d8
# ╠═3e58a59c-5267-43f5-bd6b-90626d022135
