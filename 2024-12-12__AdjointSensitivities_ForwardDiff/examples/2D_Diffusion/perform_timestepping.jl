
"""
    perform_timestepping(Params, U, F, Jac, colors, Δt;  tmax=1e3, max_timesteps=1e6, verbose=false, maxit=1000)

Performs the timestepping loop for a nonlinear system of equations where:
    - `Params` is a NamedTuple containing the physical parameters
    - `U` is a vector of abstract arrays containing the initial guess of the fields
    - `F` is the residual vector of the system
    - `Jac` is the jacobian matrix of the system
    - `colors` is the coloring matrix of the jacobian
    - `Δt` is the timestep

"""
function perform_timestepping(Params, U, F, Jac, colors, Δt, Res_closed!, N;  tmax=1e3, max_timesteps=1e6, verbose=true, maxit=1000, save_fig=false, fig_name="example_1")
    time = 0.0
    it = 0

    time_vec = Float64[]
    time     = 0;
    Usol     = copy(U)
    U0       = copy(U)
    dT_dt_vec = []

    println("Starting simulation")
    
    dT_dt = 1e3

    # #initialize figure
    qD_mag = sqrt.(lin_int(Params.qDx.^2, 2) .+ lin_int(Params.qDz.^2, 1))

    fig = Figure(; size=(1600, 1600))
    ax1 = Axis(fig[1,1], xlabel="x", ylabel="z]", title="Temperature [°C]")
    l1  = heatmap!(ax1, Params.xc, Params.zc, U[1], colormap = :lajolla)
    contour!(ax1, Params.xc, Params.zc, Params.Radc, levels=[0.005], color=:blue)
    Colorbar(fig[1,2], l1)
    ax2 = Axis(fig[1,3], xlabel="x", ylabel="flux", title=" Magnitude of Heat flux [W/m²]")
    l2  = heatmap!(ax2, Params.x[2:end-1], Params.z[2:end-1], qD_mag, colormap = :blues)
    Colorbar(fig[1,4], l2)
    ax3 = Axis(fig[2,1], xlabel="x", ylabel="D", title=" Magnitude of Diffusion coefficient [m²/s]")
    l3  = heatmap!(ax3, Params.xc, Params.zc, Params.D, colormap = :reds)
    Colorbar(fig[2,2], l3)
    display(fig)
    
    # Time-stepping loop
    while dT_dt > 1e-12 #time<=tmax && it<max_timesteps
        converged = false
        maxit = 1000;
        copyto!(U0, U)

        while !converged
            Usol, converged, its = nonlinear_solution(F, U, Jac, colors, Res_closed!, N, maxit=maxit, tol=1e-6, verbose=verbose);
            if !converged
                copyto!(U, U0)
                Δt0 = Δt
                Δt  = Δt / 2
                println("  Did NOT converge in $its steps; decreasing Δt from $Δt0 to $(Δt)")
                if Δt < 1e-20
                    println("timestep too small... ¡ABORT!")
                    break
                end
            elseif its<10  #&& Δt<max_Δt
                println("  Converged in $its steps; Δt=$Δt")
                # Δt *= 1.2
            end
        end
        
        time += Δt
        it   += 1
        
        # dT_dt = sum(((Usol[1] .- U0[1]).^2))
        dT_dt = maximum(abs.(Usol[1] .- U0[1]))
        push!(time_vec, time)
        push!(dT_dt_vec, dT_dt)

        println("  Time=$time, Δt=$Δt, dT/dt=$dT_dt")

        # Update old parameters
        Params.T_old .= copy(Usol[1])
        copyto!(U, Usol)
        
        # Visualisation
        if  it % 1 == 0#dT_dt < 1e-6

            qD_mag = sqrt.(lin_int(Params.qDx.^2, 2) .+ lin_int(Params.qDz.^2, 1))

            # #initialize figure
            fig = Figure(; size=(1600, 1600))
            ax1 = Axis(fig[1,1], xlabel="x", ylabel="z]", title="Temperature [°C]")
            l1  = heatmap!(ax1, Params.xc, Params.zc, U[1], colormap = :lajolla)
            contour!(ax1, Params.xc, Params.zc, Params.Radc, levels=[0.005], color=:blue)
            Colorbar(fig[1,2], l1)
            ax2 = Axis(fig[1,3], xlabel="x", ylabel="flux", title=" Magnitude of Heat flux [W/m²]")
            l2  = heatmap!(ax2, Params.x[2:end-1], Params.z[2:end-1], qD_mag, colormap = :blues)
            Colorbar(fig[1,4], l2)
            ax3 = Axis(fig[2,1], xlabel="x", ylabel="D", title=" Magnitude of Diffusion coefficient [m²/s]")
            l3  = heatmap!(ax3, Params.xc, Params.zc, Params.D, colormap = :reds)
            Colorbar(fig[2,2], l3)
            ax4 = Axis(fig[2,3], xlabel="Time [s]", ylabel="dT/dt", title="dT/dt")
            l4  = lines!(ax4, time_vec, log10.(dT_dt_vec), color=:black, linewidth=2)
            display(fig)
            if save_fig == true
                if !isdir("figures")
                    mkdir("figures")
                end
                save("figures/" * fig_name * "_$(Int(it/1)).png", fig)
            end
        end
        
        println("Timestep $it/$max_timesteps, time=$time, Δt=$Δt")
    end
    return time_vec, U
end