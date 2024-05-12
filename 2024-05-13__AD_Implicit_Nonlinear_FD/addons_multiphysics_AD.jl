# This contains a few routines that makes it easier to deal with multiphysics problems,
# by allowing a user to provide a tuple of arrays as input to a residual function, rather than a single array.
import Symbolics: jacobian_sparsity


"""
    Symbolics.jacobian_sparsity(func::Function, input::NTuple{NN,AbstractVector{T}}; kwargs...) where {NN<:Number, T<:Number}

Allows creating a sparsity pattern for the jacobian of a function `func` with respect to the input `input`, for the case that we have multiple input fields 
"""
function jacobian_sparsity(func::Function, input::NTuple; kwargs...)
    
    # Create a vector with variable names for the different fields
    indices = eachindex.(input) # transfer to linear indices
    vars = map(Symbolics.variable, indices[1])
    vars_tuple = (vars,)
    le   = length(vars)
    for iN=2:length(indices)
        vars_local = map(Symbolics.variable, indices[iN].+le)
        vars = [vars; vars_local ]
        vars_tuple = (vars_tuple..., vars_local)
        le += length(vars_local)
    end
   
    expr = func(vars_tuple; kwargs...)
    return Symbolics.jacobian_sparsity(expr, vars)
end

"""
    jacobian_sparsity(func::Function, output::AbstractArray, input::NTuple; kwargs...) #where {NN<:Number, T<:Number}
    
Computes sparsity pattern for the jacobian of a function `func` with respect to the input `input`, for the case that we have multiple input fields 
"""
function jacobian_sparsity(func::Function, output::AbstractArray{T}, input::Vector{<:AbstractArray{T}}; kwargs...) where {T<:Number}
    
    # Create a vector with variable names for the different fields
    indices = eachindex.(input) # transfer to linear indices
    vars = map(Symbolics.variable, indices[1])
    vars_vecarray = [vars]
    le   = length(vars)
    for iN=2:length(indices)
        vars_local = map(Symbolics.variable, indices[iN].+le)
        vars = [vars; vars_local ]
        push!(vars_vecarray, vars_local)
        le += length(vars_local)
    end
    
    expr = zero(vars)
    func(expr, vars; kwargs...)
    return Symbolics.jacobian_sparsity(expr, vars)
end


"""
    vecarray_2_vec(vecarray::Vector{AbstractArray{T}})

Creates a vector from a vector of arrays
"""
function vecarray_2_vec(vecarray::Vector{<:AbstractArray{T}}) where T
    vec = vecarray[1][:];

    for i=2:length(vecarray)
        vec = [vec; vecarray[i][:]]
    end
    
    return vec
end

"""
    vecarray_2_vec!(v::Vector, vecarray::Vector{AbstractArray{T}}, N::NTuple{NN,Int64})

Creates a vector from a tuple of arrays
"""
function vecarray_2_vec!(v::Vector{T}, vecarray::Vector{AbstractArray{T}}, N::NTuple{NN,Int64}) where {NN,T<:Number}
    ind = cumsum(prod.(N));
    v[1:ind[1]] = vecarray[1][:];
    st = ind[1]
    for i=2:NN
        v[st+1:ind[i]] = vecarray[i][:];
    end
    
    return nothing
end

"""
    vecarray = vec_2_vecarray(vec::Vector, N::NTuple)
Creates a vector of arrays from a single vector, provided the dimensions of the arrays are given as a tuple `N`
"""
function vec_2_vecarray(v::Vector{T}, N::NTuple) where T<:Number
    ind = cumsum(prod.(N));
    st = 0;
    vecarray = [reshape(v[st+1:ind[1]],N[1])]
    
    st = ind[1]
    for i=2:length(N)
        push!(vecarray, reshape(v[st+1:ind[i]],N[i]))
        st = ind[i]
    end
    
    return vecarray
end


"""
    vec_2_vecarray!(vecarray::Vector{Vector{T}}, v::Vector{T}) where {T<:Number}

In-place routine to transfer a vector `v` to a vector of arrays `array_vec`
"""
function vec_2_vecarray!(vecarray::Vector{Vector{T1}}, v::Vector{T}) where {T<:Number,T1<:Number}

    j = 1
    for i=1:length(vecarray)
        for I in eachindex(vecarray[i])
            vecarray[i][I] .= v[j];
            j += 1
        end
    end
    
    return nothing
end

"""
    LineSearch(func::Function, F, x, δx; α = [0.01 0.05 0.1 0.25 0.5 0.75 1.0])
    
Linesearch algorithm provided that the function `func` is an in-place function called as `func!(F, x)` and 
the line search parameters are provided in the structure `LS`

(Provided by Thibault Duretz)
"""
function LineSearch(func::Function, F, x, δx;  α = [0.01 0.05 0.1 0.25 0.5 0.75 1.0])
    Fnorm = zero(α)
    N = length(x)
    for i=1:length(α)
        func(F, x .+ α[i].*δx)
        Fnorm[i] = norm(F)/N
    end
    v, i_opt = findmin(Fnorm)
    return α[i_opt], Fnorm[i_opt]
end


"""
    Usol = nonlinear_solution(Fup::Vector, U::Vector{<:AbstractArray}, J, colors, tol=1e-8, itmax=100)

Computes a nonlinear solution using a Newton method with line search.
`U` needs to be a vector of abstract arrays, which contains the initial guess of every field 
`J` is the sparse jacobian matrix, and `colors` the coloring matrix, usually computed with `matrix_colors(J)`
"""
function nonlinear_solution(Fup::Vector, U::Vector{<:AbstractArray}, J, colors, tol=1e-8, itmax=100)

    r   = zero(Fup)
    Uvec = vecarray_2_vec(U)
    err = 1e3; it=0;
    while err>tol && it<itmax
        # compute residual
        Res_closed!(r,Uvec) 

        # compute jacobian in an in-place manner
        forwarddiff_color_jacobian!(J, Res_closed!, Uvec, colorvec = colors)

        # solve linear system:
        du = J\-r
        
        # use line search to find the optimal step size
        α, err =  LineSearch(Res_closed!, r, Uvec, du);
        
        Uvec = Uvec + α*du  # update solution

        it +=1;
        println("Iteration $it: error = $err")
    end
    Usol = vec_2_vecarray(Uvec, N)

    return Usol
end
