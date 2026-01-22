### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 7c87c1a0-bf7d-450c-9387-4459d84fa7f4
using Clang, Pkg, CEnum 

# ╔═╡ 35abdaca-0aed-4e6e-905d-7fd2bd0ca730
Pkg.add(path="https://github.com/boriskaus/BinaryPlayground_jll.jl")

# ╔═╡ 715f8b7d-440a-424a-bcba-73a53868c7e0
using BinaryPlayground_jll

# ╔═╡ 70c1181a-f626-11f0-2e8b-0ba2fee3fd75
md"""
# Interaction of Julia with C/Fortran binaries and libraries

Julia User Group Meeting, 22.01.2026, *Boris Kaus*
"""

# ╔═╡ 7c8322d0-afdd-4571-bf69-b7fddd7dfd02
md"""
# Outline

1. Call executables 
"""

# ╔═╡ 3efc336a-79a4-47f6-848b-8f888ecda31b
md"""
## Background

In the ideal world, all code would be written in Julia so we never have to bother with anything else
"""

# ╔═╡ 1cbc54bd-5af5-4061-b9b2-a1b0b695d5a9
md"""
## Background

But we don't live in an ideal world:

- Many people have invested years of work in their C/Fortran codes and don't want to start fully from zero. Convincing them to switch to Julia is easier if they can keep using these packages

- Many existing packages are the results of decades of work:
- Examples: [PETSc](https://petsc.org/release/) 

- Some packages have become quasi-standard in scientific computing: 
- Examples: LAPACK, BLAS, NetCDF, P4est, HDF5

- Good news: It's pretty easy to combine existing C/Fortran/Python code with julia.

- Python: [PyCall](https://github.com/JuliaPy/PyCall.jl) and [PythonCall](https://github.com/JuliaPy/PythonCall.jl)

- C/(Fortran): Let's have a look at that.
"""

# ╔═╡ b383dcf3-4c00-461c-a6b7-65bc265e5d28
md"""
## A simple C example code

[https://github.com/boriskaus/julia\_binary\_interaction\_playground](https://github.com/boriskaus/julia_binary_interaction_playground)

This small repository demonstrates a tiny C program that implements three
kinds of summation functions (scalars, vectors, and structs) and a simple
CLI to exercise them.

Files of interest
- `binary_playground.c` - implementation and `main()`.
- `binary_playground.h` - public declarations and `MyStruct` definition.
- `Makefile` - builds the executable and a dynamic library.
"""

# ╔═╡ 1edc0718-52c6-4e58-a6f3-08b9e30c3424
md"""
Build the executable:

```sh
$ make
```
"""

# ╔═╡ e0b4837a-5ac6-4626-b9e3-0197ab97c5ca
cd(@__DIR__) do
    run(`make clean`)
	run(`make all`)
end

# ╔═╡ 4d7abd1f-6575-4052-8f37-b4fc54f8145c
md"""
Build the dynamic library (macOS `.dylib` or Linux `.so`):

```sh
$ make lib
```
"""

# ╔═╡ 805f16e4-fceb-4e2a-8f84-1ebf11a8b66b
cd(@__DIR__) do
	run(`make lib`)
end

# ╔═╡ 8be485f0-ce22-4345-a6cf-dc709e4cdec7
md"""
## A simple C example code
Running this can be done in 3 ways:

Run the scalar example and supply three floats:

```sh
$ ./binary_playground --mode scalar --a 1.2 --b 3.4 --c 5.6
```

Run the built-in vector example:

```sh
$ ./binary_playground --mode vector
```

Run the built-in struct example:

```sh
$ ./binary_playground --mode struct
```
"""

# ╔═╡ e891b765-896b-4cbd-aff3-22a5dce00ae4
run(`./binary_playground --mode scalar --a 1.3`)

# ╔═╡ e0a92d18-8b43-438d-a615-e09e0bce67ab
md"""
## A simple C example code

Sum scalars:
```C
/*
 * sum_scalars
 * ---------------
 * Compute the sum of three scalar floats and return the result.
 */
float sum_scalars(float a, float b, float c) {
    return a + b + c;
}
```
"""

# ╔═╡ 1ae20986-d4f7-4a8c-a75c-111348625cd3
md"""
## A simple C example code

Add 2 vectors to the 1th one:
```C
/*
 * sum_vectors
 * ---------------
 * In-place element-wise sum of three vectors. The result is written into
 * the first vector (`v1`). All arrays are assumed to have at least `len`
 * elements.
 */
void sum_vectors(float *v1, const float *v2, const float *v3, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        v1[i] = v1[i] + v2[i] + v3[i];
    }
}
```
"""

# ╔═╡ e789e51b-f506-4d64-affa-fb445ace12dd
md"""
## A simple C example code

Sum two structs:
```C
/*
 * sum_structs
 * ---------------
 * Add corresponding numeric fields of two `MyStruct` instances and return
 * a new `MyStruct` with the aggregated values.
 */
MyStruct sum_structs(const MyStruct *s1, const MyStruct *s2) {
    MyStruct out;
    out.x = s1->x + s2->x;
    out.y = s1->y + s2->y;
    out.z = s1->z + s2->z;
    return out;
}
```

"""

# ╔═╡ 393e2aee-82d3-48a5-8fdb-1297ac495528
md"""
## Call a binary from the shell
We now have a compiled binary `binary_playground`.
This can be called from the command-line as:
```sh
$ ./binary_playground 
sum_scalars(1.000000, 2.000000, 3.000000) = 6.000000
```
"""

# ╔═╡ b0d551c6-cbd1-466f-9f0b-06d7309d8335
md"""
Or:
```sh
$ ./binary_playground --a 2.0 --b 3.3
sum_scalars(2.000000, 3.300000, 3.000000) = 8.300000
```
"""

# ╔═╡ de0ef940-e6b4-4d5c-8961-5c1ac9661aea
md"""
## Call a binary from julia
"""

# ╔═╡ 48d17943-df96-473b-98c0-06db4df39c53
cmd = `./binary_playground --a 2.8 --b 3.1`

# ╔═╡ d6e09c25-14fc-4f8d-997f-fd98ca1a1ac1
run(cmd)

# ╔═╡ 98cef729-222e-4289-a234-5ccdb84f62d9
let
	@show typeof(cmd)
 	@show cmd.env
end

# ╔═╡ 5d8f29fc-cdaa-42c4-a908-7aa587681638
foo(;a=1.0,b=2.0) = run(`./binary_playground --a $a --b $b`)

# ╔═╡ 98c7919a-12aa-40c1-8fa7-0e917db0f20c
foo(a=3.8)

# ╔═╡ f79d611c-289f-40c1-a347-92769c50c2e8
md"""
## Call the scalar function in the dynamic library using ccall
"""

# ╔═╡ 8b057ff6-67d4-4529-b16b-d2a04bce255f
lib = joinpath(@__DIR__, "libbinary_playground.dylib")   # adjust path if needed

# ╔═╡ fa22d3a3-be30-4953-9271-a5e3c9bc7d1e
# Call the sum_scalars function
res = ccall((:sum_scalars, lib), Cdouble, (Cdouble, Cdouble, Cdouble),
            1.3, 2.4, 3.0)

# ╔═╡ 34b2c52e-754f-4ed5-a7b8-f2bcf6abdbde
print(res)

# ╔═╡ cdf229bd-1705-407b-abf6-7beb5f863abd
md"""
```C
/*
 * sum_scalars
 * ---------------
 * Compute the sum of three scalar double and return the result.
 */
double sum_scalars(double a, double b, double c) {
    return a + b + c;
}
```
"""

# ╔═╡ e01ec17e-c296-4cd9-b537-1682c1a702b8
md"""
## Do the same with the vector function
"""

# ╔═╡ d5c2edd0-2319-414b-945c-4510e2026048
let # Call the sum_vectors function
	len = 5
	v1 = [1.0, 2.0, 3.0, 4.0, 5.0]
	v2 = [10.0, 20.0, 30.0, 40.0, 50.0]
	v3 = [100.0, 200.0, 300.0, 400.0, 500.0]    
	ccall((:sum_vectors, lib), Cvoid,
	      (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Csize_t),
	      v1, v2, v3, len)
	println(v1)   # v1 is now modified to contain the sums
end


# ╔═╡ 0c91adcd-870e-4be6-96cc-01b097e35824
md"""
```C
/*
 * sum_vectors
 * ---------------
 * In-place element-wise sum of three vectors. The result is written into
 * the first vector (`v1`). All arrays are assumed to have at least `len`
 * elements.
 */
void sum_vectors(double *v1, const double *v2, const double *v3, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        v1[i] = v1[i] + v2[i] + v3[i];
    }
}
```
"""

# ╔═╡ f17f819e-bac3-4d1e-acfb-2114a2a81e13
md"""
# And with a `struct`
"""

# ╔═╡ e99e96c7-9541-4559-a434-d388e52fb54b
struct MyStruct
	x::Cint
	y::Cdouble
	z::Cdouble
end 

# ╔═╡ e84ee54c-21f5-423d-b6da-8ccf8d9b3d03
let # Call the sum_vectors function
	#struct MyStruct
	#    x::Cint
	#    y::Cdouble
	#    z::Cdouble
	#end 
	s1 = MyStruct(1, 2.0, 3.0)
	s2 = MyStruct(10, 20.0, 30.0)
	s3 = ccall((:sum_structs, lib), MyStruct,
	            (Ref{MyStruct}, Ref{MyStruct}),
	            s1, s2)
end

# ╔═╡ b32100f6-cba9-4845-868e-256870c21248
md"""
```C
/*
 * sum_structs
 * ---------------
 * Add corresponding numeric fields of two `MyStruct` instances and return
 * a new `MyStruct` with the aggregated values.
 */
MyStruct sum_structs(const MyStruct *s1, const MyStruct *s2) {
    MyStruct out;
    out.x = s1->x + s2->x;
    out.y = s1->y + s2->y;
    out.z = s1->z + s2->z;
    return out;
}
```
"""

# ╔═╡ da130d6c-3171-461c-9466-7298f98e788d
md"""
## Using Clang.jl to autogenerate wrappers

Issue: manually generating wrappers for large libraries can be very time-consuming.

- Is there an automatic way to do this?
> sure, Clang.jl 
"""

# ╔═╡ 0ae37634-7677-473e-aba3-9e3c3aa0eb64
md"""
## And how do we do that?

You need to generate 3 files:

`generator.toml`:
```
[general]
library_name = "lib_BinaryPlayground"
output_file_path = "./BinaryPlayground_library.jl"
module_name = "LibBinaryPlayground"
prologue_file_path = "prologue.jl"
```
`generator.jl`:
```julia
using Clang.Generators
using Clang.LibClang.Clang_jll
using Pkg
using Pkg.Artifacts

cd(@__DIR__)

dir = pwd();  

# wrapper generator options
options = load_options(joinpath(@__DIR__, "generator.toml"))

# add compiler flags, e.g. "-DXXXXXXXXX"
args = get_default_args()

# Process all header files in the directory:
include_dir = dir
header_files = [joinpath(dir, header) for header in readdir(include_dir) if endswith(header, ".h")]

# create context
ctx = create_context(header_files, args, options)

# run generator
build!(ctx)
```
`prologue.jl`:
```julia
# START OF PROLOGUE
if isfile("libbinary_playground.dylib")
    const lib_BinaryPlayground = joinpath(pwd(),"libbinary_playground.dylib")
    println("Using locally compiled version of libbinary_playground.dylib")
else
    warning("libbinary_playground.dylib not found in current directory.")
end
```
"""

# ╔═╡ a808f26c-7e50-4734-b22e-a201ea697409
md"""
## Running Clang.jl
"""

# ╔═╡ d5511c75-f07e-45c3-84ca-c529a49a0ff8
md"""
and run it with:
"""

# ╔═╡ 188176cc-3096-4607-9451-d117f062e85a
#include("generator.jl")

# ╔═╡ 293f0857-be62-4d85-82d2-af684cbef39d
md"""
## This will do some magic 

And generates the file `BinaryPlayground_library.jl`:

```julia
module LibBinaryPlayground

using CEnum: CEnum, @cenum

# START OF PROLOGUE
if isfile("libbinary_playground.dylib")
    const lib_BinaryPlayground = joinpath(pwd(),"libbinary_playground.dylib")
    println("Using locally compiled version of libbinary_playground.dylib")
else
    warning("libbinary_playground.dylib not found in current directory.")
end

struct MyStruct
    x::Cint
    y::Cdouble
    z::Cdouble
end

function sum_scalars(a, b, c)
    ccall((:sum_scalars, lib_BinaryPlayground), Cdouble, (Cdouble, Cdouble, Cdouble), a, b, c)
end

function sum_vectors(v1, v2, v3, len)
    ccall((:sum_vectors, lib_BinaryPlayground), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Csize_t), v1, v2, v3, len)
end

function sum_structs(s1, s2)
    ccall((:sum_structs, lib_BinaryPlayground), MyStruct, (Ptr{MyStruct}, Ptr{MyStruct}), s1, s2)
end

end # module
```
"""

# ╔═╡ 42e76098-9f9e-4367-8c67-820797add3e8
#include("BinaryPlayground_library.jl")

# ╔═╡ 9dce55b0-469e-4443-a5e4-30ea45f2850c
md"""
and can call the function with:
"""

# ╔═╡ a0a70e81-7aa0-4cce-a578-1f16639f08a5
LibBinaryPlayground.sum_scalars(1.43,10,3.0)

# ╔═╡ a0da3f4e-3a7c-4284-a191-e9e22a8e3140


# ╔═╡ 053cd563-2b63-4361-aff3-1a9d5959419f
md"""
# Precompile binaries across platforms using BinaryBuilder.jl
[binarybuilder.org](https://binarybuilder.org/)
"""

# ╔═╡ a9ac73b3-e2ef-4b63-a14f-4979fa78a489
md"""
## What's BinaryBuilder.jl?
[https://github.com/JuliaPackaging/BinaryBuilder.jl](https://github.com/JuliaPackaging/BinaryBuilder.jl)
- Provides a stripped-down linux version that allows you to compile codes across platforms (windows, mac, linux). 
- Cross-compiles the packages.
- You can build executables and shared libraries
- These are made available in the Julia ecosystem as `*_jll` packages (e.g., `LaMEM_jll`, `PETSc_jll`, `GMT_jll`)
- There is a central github repo ([https://github.com/JuliaBinaryWrappers/](https://github.com/JuliaBinaryWrappers/)) to store binaries 
- Julia will download the correct one depending on your operating system, julia version etc. etc.
- You can also use them *outside* the julia ecosystem (e.g., `Perple_X_jll`)

Remarks:
- You can only distribute libraries that have an open-source license
- You can create local builds to store in your home directory
- It seems you need to use julia 1.7 to run BinaryBuilder
- If you want to make it available in the Julia package manager, you need to upload the build script to [https://github.com/JuliaPackaging/Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil) 
- Binaries will be build and hosted for free
- `Clang.jl` works very well with `*_jll` packages
"""

# ╔═╡ 1437c807-d9cd-4eeb-86ab-c92f6e98f5bd
md"""
## Lets set that up for our case

1. Use Julia 1.7!
2. Use linux (at least I have issues on Mac)!
3. We will use the github page of this package: [https://github.com/boriskaus/julia\_binary\_interaction\_playground](https://github.com/boriskaus/julia_binary_interaction_playground)
4. And use the wizard of BinaryBuilder. 
```julia
julia> BinaryBuilder.run_wizard()
```
5. This generated the file `build_tarballs.jl`:
```julia
# Note that this script can accept some limited command-line arguments, run
# `julia build_tarballs.jl --help` to see a usage message.
using BinaryBuilder, Pkg

name = "BinaryPlayground"
version = v"0.1.0"

# Collection of sources required to complete build
sources = [
    GitSource("https://github.com/boriskaus/julia_binary_interaction_playground.git", 
                "1df0d115a3eefeb9d1fbfb5e4efa9fb73d436d28")
]

# Bash recipe for building across all platforms
script = raw\"\"\"
cd $WORKSPACE/srcdir
cd julia_binary_interaction_playground/

make all
make lib

install -Dvm 755 libbinary_playground.${dlext} "${libdir}/libbinary_playground.${dlext}"
install -Dvm 755 binary_playground${exeext} "${bindir}/binary_playground${exeext}"

install_license LICENSE
\"\"\"

# These are the platforms we will build for by default, unless further
# platforms are passed in on the command line
platforms = [
    Platform("i686", "linux"; libc = "glibc"),
    Platform("x86_64", "linux"; libc = "glibc"),
    Platform("aarch64", "linux"; libc = "glibc"),
    Platform("armv6l", "linux"; call_abi = "eabihf", libc = "glibc"),
    Platform("armv7l", "linux"; call_abi = "eabihf", libc = "glibc"),
    Platform("powerpc64le", "linux"; libc = "glibc"),
    Platform("riscv64", "linux"; libc = "glibc"),
    Platform("i686", "linux"; libc = "musl"),
    Platform("x86_64", "linux"; libc = "musl"),
    Platform("aarch64", "linux"; libc = "musl"),
    Platform("armv6l", "linux"; call_abi = "eabihf", libc = "musl"),
    Platform("armv7l", "linux"; call_abi = "eabihf", libc = "musl"),
    Platform("x86_64", "macos"; ),
    Platform("aarch64", "macos"; ),
    Platform("x86_64", "freebsd"; ),
    Platform("aarch64", "freebsd"; )
]


# The products that we will ensure are always built
products = [
    LibraryProduct("libbinary_playground", :libbinary_playground),
    ExecutableProduct("binary_playground", :binary_playground)
]

# Dependencies that must be installed before this package can be built
dependencies = Dependency[
]

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies; julia_compat="1.6", preferred_gcc_version = v"12.1.0")
```
"""

# ╔═╡ c81de3e4-8f1d-43ad-a8f5-f88cc3f83909
md"""
## Build it for several versions and deploy it to your own github repository

This is extremely useful, as it allows you to test the binaries locally, before making it available in the julia registry:

```julia
julia1.7 build_tarballs.jl --debug --verbose --deploy="boriskaus/BinaryPlayground_jll.jl" 
```

This will push the results to:

[https://github.com/boriskaus/BinaryPlayground_jll.jl](https://github.com/boriskaus/BinaryPlayground_jll.jl)

Remarks:
- If you rebuild this, you'll have to delete the release first on your github page.
- You can create a `Personal access token` for your github repo, so you don't have to create a password all the time.
"""

# ╔═╡ 0eb99dd9-5ab4-4c13-8155-904260257a32
md"""
## Using the `BinaryPlayground_jll` from your local github

```julia
julia>]
pkg> add https://github.com/boriskaus/BinaryPlayground_jll.jl
julia> using BinaryPlayground_jll
```
"""

# ╔═╡ 673b8909-0822-4dd3-b155-ba702e3c88b6
run(`$(BinaryPlayground_jll.binary_playground())`)

# ╔═╡ 73d5b622-a9de-436d-8374-d272526dce0f


# ╔═╡ 6ad490c3-2924-48f2-83aa-7cd02b232126


# ╔═╡ 40e219fa-5840-4312-8d6a-e3a7bb5798d3
md"""
## Making codes available in Yggdrasil 

Yggdrasil is a central Julia repository with buildscripts for all available `_jll` packages.

[https://github.com/JuliaPackaging/Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil)

Have a look what is available!

If you want your package to be added, you'll need to create a pull request there.

- Be prepared to do some extra work
- In case of problems, use the BinaryBuilder slack channel
- You will need to use other binary dependencies as much as possible
- Doing this can be quite easy (`MAGEMin_jll`) or very complicated/time-consuming (`PETSc_jll`)
"""

# ╔═╡ 2d233a53-1fd7-4e95-b664-06ec17c85225
md"""
## Some tips: calling two executables from _jll files 

In some cases you want to run 2 executables at the same time. 

Example: `mpi` and `LaMEM`.

##### Problem: 
Each of them comes with their own environment.

##### Solution: 
add LaMEM environment to `mpirun`: 
```julia
using LaMEM_jll, MPICH_jll

const MPI_LIBPATH = LaMEM_jll.MPICH_jll.LIBPATH
key = LaMEM_jll.JLLWrappers.JLLWrappers.LIBPATH_env
mpirun = addenv(mpiexec, key=>join((LaMEM_jll.LIBPATH[], MPI_LIBPATH[]), pathsep));

# create command-line object
cmd = `$(mpirun) -n $cores_compute $(LaMEM_jll.LaMEM().exec) -ParamFile $(ParamFile) $args`
```

"""

# ╔═╡ f4704dee-7d77-4ec2-be8e-190095d0ae9c
md"""
## What about C++ codes?

- C++ dynamic libraries don't expose the functions in the dynamic library
- You'll need to create C wrappers for the function calls
- [I haven't tried this yet; needs to be done to fully wrap `LaMEM`]
"""

# ╔═╡ 4c4e7b5a-fa3d-45f3-a82c-2749cb701ba9
md"""
# Examples of codes we have wrapped
"""

# ╔═╡ 9c3fc812-cd4b-48e3-bbaa-6bb723e74c69
md"""
## Some examples:
- [PETSc\_jll](https://github.com/JuliaBinaryWrappers/PETSc_jll.jl) | MPI + PETSc binaries uploaded in Yggdrasil. C library and full wrapper (semi-automatically) created in [PETSc.jl](https://github.com/JuliaParallel/PETSc.jl/tree/main)
- [MAGEMin\_jll](https://github.com/JuliaBinaryWrappers/MAGEMin_jll.jl) | C code to solve thermodynamics of rocks and melt. Fully wrapped using Clang. Distributed in [MAGEMin_C.jl](https://github.com/ComputationalThermodynamics/MAGEMin_C.jl), with every function available (in principle).  
- [LaMEM\_jll](https://github.com/JuliaBinaryWrappers/LaMEM_jll.jl) | MPI-parallel C(++) code to simulate 3D geodynamic processes; Binaries distributed using Yggdrasil. We only run the executable and the [LaMEM.jl](https://github.com/JuliaGeodynamics/LaMEM.jl) package creates an input script.
- [FastScapelib\_jll](https://github.com/JuliaBinaryWrappers/Fastscapelib_jll.jl) | Erosion/surface processes code written in Fortran. Interface written manually but covers the full code. See [FastScape.jl](https://github.com/boriskaus/FastScape.jl)
- [PhreeqcRM\_jll](https://github.com/JuliaBinaryWrappers/PhreeqcRM_jll.jl) - Aqueous chemistry code written in ?. Distributed in Yggdrasil. Julia interface in PhreeqcRM.jl (currently private).
- OpenSWPC | 3D MPI-parallel code to simulate wave propagation. Precompiled binary executable for Mac/Linux available in the github of Boris.  OpenSWPC.jl calls these binaries and provides an interface to the GeophysicalModelGenerator.jl
- [Perple\_X\_jll](https://github.com/JuliaBinaryWrappers/Perple_X_jll.jl) - precompiled binaries for the thermodynamics code Perple_X (Fortran).
- [GMT\_jll](https://github.com/JuliaBinaryWrappers/GMT_jll.jl). Precompiled binarties for the Generic Mapping Tools, a widely used plotting toolbox in the geosciences, with julia interface available in [GMT.jl](https://github.com/GenericMappingTools/GMT.jl)
- GemPy.jl | Python package to perform uncertainty modelling. Julia interface using PythonCall (currently private).
"""

# ╔═╡ 009809be-12b8-4aad-bb13-f654762e39e9
md"""
## What about `PETSc.jl`?

- [PETSc.jl](https://github.com/JuliaParallel/PETSc.jl) 0.4.0 just released 
- Wraps nearly the full PETSc library
- Over 3000 functions, 250k lines of code.
- We started working on that in 2021
- Not using Clang to create wrappers, but created some code to automatically generate function interfaces
- Provides docstrings to every function (translated from the C docs)
- Uses a python library that lists all functions, input/oputput arguments, structs and tuypes in PETSC
- will talk about that some other day
"""

# ╔═╡ d2ab2cd1-f659-47f0-87c9-f6ae2fd89912
md"""
# Going the other way: Create executables and dynamic libraries from Julia

- We have written some amazing piece of Julia code ([RheologyCalculator.jl](https://github.com/albert-de-montserrat/RheologyCalculator.jl), for example).

- But our colleagues have an ancient (efficient) C code.

- It remains challenging to convince them to fully switch to julia

- What can be done?

"""

# ╔═╡ 6725aa97-c688-4ed9-8f4d-38b3aa493019
md"""
## `juliac` compiler!

- the `juliac` compiler allows compiling any piece of julia code into a standalone executable  
- the new `--trim` option released in julia 1.12 reduces the size of the executable
- there was a very nice, compact, example on the Julia discourse by `mkitti`:
[https://discourse.julialang.org/t/self-contained-juliac-demonstration/134366](https://discourse.julialang.org/t/self-contained-juliac-demonstration/134366)

Below, I have modified the example `multiply_numbers.jl` to:
1. Include a separate function to perform the actual calculation (`perform_calculation`)
2. Create a C-interface for this function, visible from within the dynamic library
3. Compile a dynamic library as well

```julia
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
    program = \"\"\"
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
    \"\"\"
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
```

"""

# ╔═╡ 1d276df9-2c0f-4c91-b191-05d9d616e368
md"""
## How to embed this dynamic library in a c-code?

Once compiled, you can check that the function is available in the dynamic library with (on mac):
```sh
mac21-003:2026-01-22_CallingBinaryCodes kausb$ nm -g build/lib/multiply_numbers.dylib
                 U _bzero
000000000000905c T _c_perform_calculation
000000000003117c T _get_jl_RTLD_DEFAULT_handle_addr
                 U _ijl_apply_generic
```

A simple C-code that calls this is `test_c_call_multiply_numbers.c`:
```C
/*
 * Test program to call the c_perform_calculation function from the multiply_numbers.dylib library.
 *
 * Instructions:
 * 1. Ensure the multiply_numbers.dylib library is built in build/lib/ (run julia multiply_numbers.jl).
 * 2. Compile this program: gcc test_c_call_multiply_numbers.c -o test_c_call_multiply_numbers
 * 3. Run the program: ./test_c_call_multiply_numbers
 *    Expected output: Result: 450
 *
 * This demonstrates calling a Julia @ccallable function from C using dlopen/dlsym.
 */

#include <stdio.h>
#include <dlfcn.h>

int main() {
    void* lib = dlopen("build/lib/multiply_numbers.dylib", RTLD_LAZY);
    if (!lib) {
        fprintf(stderr, "Error loading library: %s\n", dlerror());
        return 1;
    }
    
    long long (*c_perform_calculation)(long long*, int) = dlsym(lib, "c_perform_calculation");
    if (!c_perform_calculation) {
        fprintf(stderr, "Error finding symbol: %s\n", dlerror());
        dlclose(lib);
        return 1;
    }
    
    long long numbers[] = {5, 9, 10};
    long long result = c_perform_calculation(numbers, 3);
    printf("Result: %lld\n", result);  // Should print 450
    
    dlclose(lib);
    return 0;
}
```

"""

# ╔═╡ ebcd2e41-5bbd-4f23-a84c-7f01fd1c8170
md"""
## Compile and run this:

```sh
mac21-003:2026-01-22_CallingBinaryCodes kausb$ gcc test_c_call_multiply_numbers.c -o test_c_call_multiply_numbers
```
```sh
mac21-003:2026-01-22_CallingBinaryCodes kausb$ ./test_c_call_multiply_numbers 
Result: 450
```
"""

# ╔═╡ 159f5104-e267-4cd9-a350-1f565a25a1dd
md"""
## Let's modify the julia code

We will only regenerate the dynamic library from julia, but we will NOT recompile `test_c_call_multiply_numbers`!

```julia
function perform_calculation(numbers::Vector{Int})
    return prod(numbers)*2
end
```
```julia
[ Info: Running build/bin/multiply_numbers.exe 5 9 10
900
Process(`build/bin/multiply_numbers.exe 5 9 10`, ProcessExited(0))
```

Without recompiling `test_c_call_multiply_numbers` we now get:
```sh
mac21-003:2026-01-22_CallingBinaryCodes kausb$ ./test_c_call_multiply_numbers 
Result: 900
```
"""

# ╔═╡ b15c961f-fc62-4a38-ad58-fc29f2947491
md"""
# Summary: you have now seen how to link julia and C codes!
"""

# ╔═╡ Cell order:
# ╟─70c1181a-f626-11f0-2e8b-0ba2fee3fd75
# ╠═7c8322d0-afdd-4571-bf69-b7fddd7dfd02
# ╟─3efc336a-79a4-47f6-848b-8f888ecda31b
# ╟─1cbc54bd-5af5-4061-b9b2-a1b0b695d5a9
# ╟─b383dcf3-4c00-461c-a6b7-65bc265e5d28
# ╟─1edc0718-52c6-4e58-a6f3-08b9e30c3424
# ╟─e0b4837a-5ac6-4626-b9e3-0197ab97c5ca
# ╟─4d7abd1f-6575-4052-8f37-b4fc54f8145c
# ╟─805f16e4-fceb-4e2a-8f84-1ebf11a8b66b
# ╟─8be485f0-ce22-4345-a6cf-dc709e4cdec7
# ╠═e891b765-896b-4cbd-aff3-22a5dce00ae4
# ╟─e0a92d18-8b43-438d-a615-e09e0bce67ab
# ╟─1ae20986-d4f7-4a8c-a75c-111348625cd3
# ╟─e789e51b-f506-4d64-affa-fb445ace12dd
# ╟─393e2aee-82d3-48a5-8fdb-1297ac495528
# ╟─b0d551c6-cbd1-466f-9f0b-06d7309d8335
# ╟─de0ef940-e6b4-4d5c-8961-5c1ac9661aea
# ╠═48d17943-df96-473b-98c0-06db4df39c53
# ╠═d6e09c25-14fc-4f8d-997f-fd98ca1a1ac1
# ╠═98cef729-222e-4289-a234-5ccdb84f62d9
# ╠═5d8f29fc-cdaa-42c4-a908-7aa587681638
# ╠═98c7919a-12aa-40c1-8fa7-0e917db0f20c
# ╟─f79d611c-289f-40c1-a347-92769c50c2e8
# ╠═8b057ff6-67d4-4529-b16b-d2a04bce255f
# ╠═fa22d3a3-be30-4953-9271-a5e3c9bc7d1e
# ╠═34b2c52e-754f-4ed5-a7b8-f2bcf6abdbde
# ╟─cdf229bd-1705-407b-abf6-7beb5f863abd
# ╟─e01ec17e-c296-4cd9-b537-1682c1a702b8
# ╠═d5c2edd0-2319-414b-945c-4510e2026048
# ╟─0c91adcd-870e-4be6-96cc-01b097e35824
# ╟─f17f819e-bac3-4d1e-acfb-2114a2a81e13
# ╠═e99e96c7-9541-4559-a434-d388e52fb54b
# ╠═e84ee54c-21f5-423d-b6da-8ccf8d9b3d03
# ╟─b32100f6-cba9-4845-868e-256870c21248
# ╟─da130d6c-3171-461c-9466-7298f98e788d
# ╟─0ae37634-7677-473e-aba3-9e3c3aa0eb64
# ╟─a808f26c-7e50-4734-b22e-a201ea697409
# ╠═7c87c1a0-bf7d-450c-9387-4459d84fa7f4
# ╟─d5511c75-f07e-45c3-84ca-c529a49a0ff8
# ╠═188176cc-3096-4607-9451-d117f062e85a
# ╟─293f0857-be62-4d85-82d2-af684cbef39d
# ╠═42e76098-9f9e-4367-8c67-820797add3e8
# ╟─9dce55b0-469e-4443-a5e4-30ea45f2850c
# ╠═a0a70e81-7aa0-4cce-a578-1f16639f08a5
# ╠═a0da3f4e-3a7c-4284-a191-e9e22a8e3140
# ╟─053cd563-2b63-4361-aff3-1a9d5959419f
# ╟─a9ac73b3-e2ef-4b63-a14f-4979fa78a489
# ╟─1437c807-d9cd-4eeb-86ab-c92f6e98f5bd
# ╟─c81de3e4-8f1d-43ad-a8f5-f88cc3f83909
# ╟─0eb99dd9-5ab4-4c13-8155-904260257a32
# ╠═35abdaca-0aed-4e6e-905d-7fd2bd0ca730
# ╠═715f8b7d-440a-424a-bcba-73a53868c7e0
# ╠═673b8909-0822-4dd3-b155-ba702e3c88b6
# ╠═73d5b622-a9de-436d-8374-d272526dce0f
# ╠═6ad490c3-2924-48f2-83aa-7cd02b232126
# ╟─40e219fa-5840-4312-8d6a-e3a7bb5798d3
# ╟─2d233a53-1fd7-4e95-b664-06ec17c85225
# ╟─f4704dee-7d77-4ec2-be8e-190095d0ae9c
# ╟─4c4e7b5a-fa3d-45f3-a82c-2749cb701ba9
# ╟─9c3fc812-cd4b-48e3-bbaa-6bb723e74c69
# ╟─009809be-12b8-4aad-bb13-f654762e39e9
# ╟─d2ab2cd1-f659-47f0-87c9-f6ae2fd89912
# ╟─6725aa97-c688-4ed9-8f4d-38b3aa493019
# ╟─1d276df9-2c0f-4c91-b191-05d9d616e368
# ╠═ebcd2e41-5bbd-4f23-a84c-7f01fd1c8170
# ╟─159f5104-e267-4cd9-a350-1f565a25a1dd
# ╠═b15c961f-fc62-4a38-ad58-fc29f2947491
