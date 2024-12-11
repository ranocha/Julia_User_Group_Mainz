
function compute_λ(J_forward, ∇_F)
    return - J_forward' \ ∇_F
end

function dR_dp(Params_array, U, U_Tuple, Δ, N, N_params, BC, Params, Δt)

    Res_adjoint_closed! = (R_adjoint, Params_up) -> Res_adjoint!(R_adjoint, Params_up, U_Tuple, Δ, N, N_params, BC, Params, Δt)

    Params_up = vecarray_2_vec(Params_array)

    # Ensure R_adjoint is correctly sized based on the forward problem residual
    R_adjoint = similar(vecarray_2_vec(U))
    
    # Calculate the adjoint residual for the given Params_array
    Res_adjoint_closed!(R_adjoint, Params_array)
    
    # Compute the Jacobian of the adjoint residual with respect to the parameters
    J_adjoint = ForwardDiff.jacobian(Res_adjoint_closed!, R_adjoint, Params_up)

    return J_adjoint
end

function adjoint_sensitivities(J_forward, U, Δ, N, N_params, BC, Params, Δt; IDX=100, IDX1=nothing)

    println("Computing sensitivities for node(s) $IDX")

    U_Tuple = (U=U[1],)
    Params_array = [Params.D]
    R_adjoint = vecarray_2_vec(U)

    # Compute the gradient of the objective function w.r.t. U
    # in the case that we want the sensitivites this is just a mask with 1's and 0's

    ∇_F = zeros(size(R_adjoint))
    if typeof(IDX) == Int
        ∇_F[IDX] = 1.0  # Sensitivity with respect to the node #
    else
        ∇_F[IDX] .= 1.0  # Sensitivity with respect to the node #
    end

    if IDX1 != nothing
        ∇_F[IDX1] .= 1.0  # Sensitivity with respect to the node #
    end

    # Solve the adjoint problem J^T * lambda = -grad_F
    λ = compute_λ(J_forward, ∇_F)

    # Compute dR/dp (Jacobian of residual w.r.t. parameters)
    J_adjoint = dR_dp(Params_array, U, U_Tuple, Δ, N, N_params, BC, Params, Δt)

    # Compute the sensitivity dF/dp
    sensitivities = -λ' * J_adjoint

    return sensitivities, λ, J_adjoint
end