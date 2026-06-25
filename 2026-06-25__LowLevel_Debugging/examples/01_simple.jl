function example()
	r = Ref(0)
	@ccall jl_breakpoint(r::Any)::Cvoid
	for i in 1:10
            r[] += i
        end
	return r
end
