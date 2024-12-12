
update_density!(ρ, phase, rheology) = update_field!(ρ, phase, rheology, compute_density)
update_viscosity!(η, phase, rheology) = update_field!(η, phase, rheology, compute_viscosity)

function non_lin_D!(D, T; dims=1)
    D .= lin_int(1.0 .* T,dims) .^2
end

function non_lin_D_2D(D, T)
    D = zeros(eltype(D), size(D))
    D = 1.0 .* T .^2
end

function non_lin_D_2D!(D, T, D0, n)
    D .= D0 .* T .^n
end

function update_field!(A, phase, rheology, f::F) where F
    args = (;)
    for i in eachindex(A)
        A[i] = f(rheology, phase[i], args)
    end
    return nothing
end

average_2D(x) = (x[2:end,2:end] + x[1:end-1,2:end] + x[2:end,1:end-1] + x[1:end-1,1:end-1])/4

function lin_int(A, dim)
    if dim == 1
        A = (A[1:(end - 1), :] .+ A[2:end, :]) ./ 2.0
    elseif dim == 2
        A = (A[:, 1:(end - 1)] .+ A[:, 2:end]) ./ 2.0
    end
    return A
end

function calculate_average_row(A, dim)
    avg_A = zeros(size(A,dim))
    for i in eachindex(avg_A)
        if dim == 1
            avg_A[i] = mean(A[i, :])
        elseif dim == 2
            avg_A[i] = mean(A[:, i])
        end
    end
    return avg_A
end

function harmonic_avg(η1, η2)
    return 2.0 * (η1 * η2) / (η1 + η2)
end

function smooth_viscosity_harmonic(η, Nx, Nz)
    for i in 2:Nx-1
        for j in 2:Nz-1
            η[i, j] = harmonic_avg(
                harmonic_avg(η[i, j], η[i+1, j]),
                harmonic_avg(η[i, j], η[i, j+1])
            )
        end
    end
    return η
end

function smooth_viscosity(η, Nx, Nz)
    μ_smooth = copy(η)
    for i in 2:Nx-1
        for j in 2:Nz-1
            μ_smooth[i, j] = 0.2 * (η[i, j] + η[i+1, j] + η[i-1, j] + η[i, j+1] + η[i, j-1])
        end
    end
    μ_smooth[:,1] .= μ_smooth[:,2]
    μ_smooth[:,end] .= μ_smooth[:,end-1]
    μ_smooth[1,:] .= μ_smooth[2,:]
    μ_smooth[end,:] .= μ_smooth[end-1,:]
    return μ_smooth
end

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function calculate_η(μ, ϕ, C, R, Peff, λ)
    η = zeros(eltype(Peff),size(μ))
    η .= μ ./ ϕ .* C .* (1.0 .+ 0.5 .* (1.0 ./ R - 1.0) .* (1.0 .+ tanh.(- Peff ./ λ)))
    return η
end

function update_stresses!(τxx, τzz, τxz, εxx,εzz,εxz, P, Params, rheology, Δt)
    I1   = CartesianIndex(1,0)
    I2   = CartesianIndex(1,1)
    I3   = CartesianIndex(0,1)
    τxz .= 0.0
    for I in CartesianIndices(εxx)
        # old stresses
        τ_o = (
            Params.τxx_old[I],
            Params.τzz_old[I],
            (Params.τxz_old[I], Params.τxz_old[I+I1], Params.τxz_old[I+I2], Params.τxz_old[I+I3]),
        )
        ε_phase = (εxx[I], εzz[I], (εxz[I], εxz[I+I1], εxz[I+I2], εxz[I+I3]))
        phase_i = (Params.phasec[I], Params.phasec[I], (Params.phasev[I], Params.phasev[I+I1], Params.phasev[I+I2], Params.phasev[I+I3]))
        
        args = (dt=Δt, P = P[I])
        τij, τII = compute_τij(rheology, ε_phase, args, τ_o, phase_i)
        
        τxx[I   ]  = τij[1]
        τzz[I   ]  = τij[2]
        # we average shear stresses @ vertices
        τxz[I   ] += τij[3]*0.25;
        τxz[I+I1] += τij[3]*0.25;
        τxz[I+I2] += τij[3]*0.25;
        τxz[I+I3] += τij[3]*0.25;
    end
    return nothing
end

function update_stresses_visc!(
    τxx, τzz, τxz, 
    εxx, εzz, εxz, 
    ηv,
    P
)
    Nx, Nz = size(εxx)

    # Reset shear stresses at the vertices before updating
    τxz .= 0.0  

    ηc = average_2D(ηv)
    # Update normal stresses τxx and τzz
    for i in 1:Nx
        for j in 1:Nz
            local_viscosity = ηc[i,j]
            τxx[i,j] = 2.0 * local_viscosity * εxx[i,j] - (1.0 / 3.0) * P[i,j]
            τzz[i,j] = 2.0 * local_viscosity * εzz[i,j] - (1.0 / 3.0) * P[i,j]
        end
    end

    # Update shear stress τxz at the vertices with averaged viscosity
    for i in 1:Nx+1
        for j in 1:Nz+1
            τxz[i,j] = 2 * ηv[i,j] * εxz[i,j]
        end
    end

    # Handle the boundary conditions explicitly for shear stress (free-slip BC)
    τxz[1, :] .= 0.0
    τxz[:, 1] .= 0.0
    τxz[end, :] .= 0.0
    τxz[:, end] .= 0.0

    return nothing
end

function initialize_helper_fields(Nx, Nz)
    # Helper fields
    εxx     =  zeros(Nx,Nz)
    εzz     =  zeros(Nx,Nz)
    εII     =  zeros(Nx,Nz)
    εxz     =  zeros(Nx+1,Nz+1)
    τxx     =  zeros(Nx,Nz)
    τzz     =  zeros(Nx,Nz)
    τII     =  zeros(Nx,Nz)
    τxz     =  zeros(Nx+1,Nz+1)
    τxx_old =  zeros(Nx,Nz)
    τzz_old =  zeros(Nx,Nz)
    τxz_old =  zeros(Nx+1,Nz+1)
    phasec  =  ones(Int64, Nx,Nz)      # material phase tag @ center
    phasev  =  ones(Int64, Nx+1,Nz+1)  # material phase tag @ vertices
    qDx     =  zeros(Nx-1,Nz)
    qDz     =  zeros(Nx,Nz-1)
    P_eff_old   =  zeros(Nx,Nz)
    return εxx, εzz, εII, εxz, τxx, τzz, τII, τxz, τxx_old, τzz_old, τxz_old, phasec, phasev, qDx, qDz, P_eff_old
end

"""
    add_box!(field, val, box_height, box_width, x_center, z_center, xc, zc)

Add a box of the value `val` to the field `field` with dimensions 
`box_height` and `box_width` centered at (`x_center`, `z_center`)
where `xc` and `zc` are the coordinates of the center points of the grid.
"""
function add_box!(field, val, box_height, box_width, x_center, z_center, xc, zc)
    # Iterate over the center points
    for (i, xv) in enumerate(xc), (j, zv) in enumerate(zc)
        # Define the box centered at (x_center, z_center)
        if abs(xv - x_center) < box_width / 2 && abs(zv - z_center) < box_height / 2
            field[i, j] = val
        end
    end
end

function limit_porosity!(ϕ)
    clamp!(ϕ, 1e-3, 1.0)
end