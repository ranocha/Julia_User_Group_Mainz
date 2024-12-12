
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
    
    dT_dt = 1

    # #initialize figure
    fig = Figure(; size=(1600, 1600))
    ax1 = Axis(fig[1,1], xlabel="x", ylabel="Temperature [°C]", title="Temperature [°C]")
    l1  = lines!(ax1, Params.xc, U[1], color=:red, linewidth=2)
    ax2 = Axis(fig[1,2], xlabel="x", ylabel="flux", title="Heat flux [W/m²]")
    l2  = lines!(ax2, Params.x[2:end-1], Params.qDx, color=:blue, linewidth=2)
    ax3 = Axis(fig[2,1], xlabel="x", ylabel="D", title="Diffusion coefficient [m²/s]")
    l3  = lines!(ax3, Params.x[2:end-1], Params.D, color=:green, linewidth=2)
    display(fig)
    
    # Time-stepping loop
    while dT_dt > 1e-5 #time<=tmax && it<max_timesteps
        converged = false
        maxit = 1000;
        copyto!(U0, U)

        while !converged
            Usol, converged, its = nonlinear_solution(F, U, Jac, colors, Res_closed!, N, maxit=maxit, tol=1e-6, verbose=verbose);
            if !converged
                copyto!(U, U0)
                Δt0 = Δt
                Δt  = Δt / 2
                println("  Did NOT converged in $its steps; decreasing Δt from $Δt0 to $(Δt)")
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
        
        dT_dt = maximum(abs.(Usol[1] .- U0[1]))
        push!(time_vec, time)
        push!(dT_dt_vec, dT_dt)

        println("  Time=$time, Δt=$Δt, dT/dt=$dT_dt")

        # Update old parameters
        Params.T_old .= copy(Usol[1])
        copyto!(U, Usol)
        
        # Visualisation
        if it % 1 == 0

            # #initialize figure
            fig = Figure(; size=(1600, 1600))
            ax1 = Axis(fig[1,1], xlabel="x", ylabel="Temperature [°C]", title="Temperature [°C]")
            l1  = lines!(ax1, Params.xc, U[1], color=:red, linewidth=2)
            ax2 = Axis(fig[1,2], xlabel="x", ylabel="flux", title="Heat flux [W/m²]")
            l2  = lines!(ax2, Params.x[2:end-1], Params.qDx, color=:blue, linewidth=2)
            ax3 = Axis(fig[2,1], xlabel="x", ylabel="D", title="Diffusion coefficient [m²/s]")
            l3  = lines!(ax3, Params.x[2:end-1], Params.D, color=:green, linewidth=2)
            ax4 = Axis(fig[2,2], xlabel="Time [s]", ylabel="dT/dt", title="dT/dt")
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