#!/usr/bin/env -S julia
using Pkg

if !contains(ENV["PATH"], DEPOT_PATH[1]*"/bin")
    @warn "The environment variable PATH does not contain $(DEPOT_PATH[1]*"/bin"). Consider modifying your shell startup."
    @info "Adding $(DEPOT_PATH[1])/bin to \$PATH" ENV["PATH"]
    ENV["PATH"] = "$(DEPOT_PATH[1])/bin:$(ENV["PATH"])"
end

if isnothing(Sys.which("juliac"))
    @info "Installing JuliaC"
    Pkg.Apps.add("JuliaC")
end

if !isdir("multiply_numbers")
    @info "Generating multiply_numbers package"
    Pkg.generate("multiply_numbers")
end

code_file::String = "multiply_numbers/src/multiply_numbers.jl"
if !isfile(code_file) || !contains("@main", read(code_file, String))
    @info "Writing $code_file"
    program = """
    module multiply_numbers

    function perform_calculation(numbers::Vector{Int})
        return prod(numbers)
    end

    Base.@ccallable function c_perform_calculation(numbers::Ptr{Int}, len::Int)::Int
        vec = unsafe_wrap(Vector{Int}, numbers, len)
        return perform_calculation(vec)
    end

    function (@main)(args::Vector{String})
        numbers = parse.(Int, args)
        out = perform_calculation(numbers)  # call function to multiply numbers
        println(Core.stdout, out)
        return 0
    end

    end # module multiply_numbers
    """
    println(program)
    write("multiply_numbers/src/multiply_numbers.jl", program)
end

@info "Compiling multiply_numbers.exe with juliac multiply_numbers --output-exe multiply_numbers.exe --bundle build --trim"
run(`juliac multiply_numbers --output-exe multiply_numbers.exe --bundle build --trim`)

@info "Compiling multiply_numbers.lib with juliac multiply_numbers --output-lib multiply_numbers.dylib --bundle build --trim"
run(`juliac multiply_numbers --output-lib multiply_numbers.dylib --bundle build --trim `)


@info "The size of build/bin/multiply_numbers.exe is $(filesize("build/bin/multiply_numbers.exe")/1024^2) MiB"

@info "Running build/bin/multiply_numbers.exe 5 9 10"
run(`build/bin/multiply_numbers.exe 5 9 10`)