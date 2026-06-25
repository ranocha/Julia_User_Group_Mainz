function without_bounds(A, i)
    @inbounds A[i] = 42
end

A = zeros(1024)
without_bounds(A, 1024)
without_bounds(A, 1025)
