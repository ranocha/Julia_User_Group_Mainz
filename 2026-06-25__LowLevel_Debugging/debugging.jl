### A Pluto.jl notebook ###
# v1.0.1

#> [frontmatter]
#> chapter = "2"
#> section = "3"
#> order = "7"
#> title = "Debugging"
#> date = "2025-05-28"
#> tags = ["module2", "track_principles"]
#> layout = "layout.jlhtml"
#> indepth_number = "2"
#> 
#>     [[frontmatter.author]]
#>     name = "Valentin Churavy"
#>     url = "https://vchuravy.dev"

using Markdown
using InteractiveUtils

# ╔═╡ 5f6deede-7e22-4ebf-ae1a-14d584595f17
begin
	using PlutoUI, PlutoTeachingTools
	using PlutoUI: Slider
	PlutoUI.TableOfContents(; depth=4)
end

# ╔═╡ bc354629-dd2a-4b5a-8e9c-378dcc357575
using Asciicast

# ╔═╡ 9182a8ea-7b06-4e4f-b0a3-35aa4cbfa456
using HypertextLiteral

# ╔═╡ aa46fef4-0c92-4947-ac64-f06ee31cb43f
ChooseDisplayMode()

# ╔═╡ d5799787-aa74-4966-b647-bab4fe1d0c3c
html"""
<h1> Low-Level debugging of Julia</h1>

<div style="text-align: center;">
Jun 25th 2025 <br>
Julia User Group
<br><br>
Valentin Churavy
<br><br>
High Performance Scientific Computing, University of Augsburg <br>
Numerical Mathematics, University of Mainz <br>
High-Performance Computing, Universität des Saarlandes
<br><br>
@vchuravy
</div>

 <table>
        <tr>
<td><img src="https://upload.wikimedia.org/wikipedia/commons/e/e9/Universit%C3%A4t_Augsburg_logo.svg" width=100px></td>

<td><img src="https://upload.wikimedia.org/wikipedia/commons/8/8a/Johannes_Gutenberg-Universit%C3%A4t_Mainz_logo.svg" width=200px></td>

<td><img src="https://upload.wikimedia.org/wikipedia/de/b/bc/Logo-Universit%C3%A4t_des_Saarlandes.svg" width=200px></td>
        </tr>
</table>
"""

# ╔═╡ 0be73b29-7780-4be0-bf07-9b62c99fc4b4
md"""
# Debugging
"""

# ╔═╡ 1109474f-1cf7-4293-b952-d80d15273de5
md"""
## Strategies / Ideas

- "Scientific method"
  - Form a hypothesis
  - Aim to invalidate it
  - One change at a time!
- Reduce, reduce, reduce!
  - The smaller the section of code, the easier it is to have a mental model of what is happeneing
- Create a git-repository \
  Goal: Keep track of changes we have so that we don't lose a bug.
- Is it "randomly" happening?

**What do we do after we fixed a bug?**
- Turn our MWE into a test! 
- The bug likely happened due to some confusion of assumptions! Does this pattern repeat elsewhere?
"""

# ╔═╡ d5e479ec-7234-415b-845d-07f29f68eeaf
md"""
### The importance of an MWE

Remember last weeks lecture on reproducibility! You might need to ask a friend or stranger for help!

- What version of Julia are you version
- Which packages
  - Project.toml / Manifest.toml

!!! note
    Make it as easy as possible for another person to reproduce your setup!
"""

# ╔═╡ c73a1c40-6426-4ced-a70b-909e54be4ced
md"""
### Common causes

- [Noteworthy Differences](https://docs.julialang.org/en/v1/manual/noteworthy-differences/)

##### Undefined memory

Julia allows for memory to be unitialized

```julia
x = Array{Float64}(undef, 10)
# or
similar(x)
```
"""

# ╔═╡ ef93ae34-fda8-488e-945b-6a3120de93c1
sum(Array{Float64}(undef, 10))

# ╔═╡ 8984c5c1-11a4-40c0-986d-2ed7db543098
let
	x = ones(10)
	similar(x)
end

# ╔═╡ 21185af9-6e23-4881-ba90-0fdd7af16471
md"""
this can "poison" your calculations
"""

# ╔═╡ a2562dc3-ace7-47f4-a121-f83b500e5af0
md"""
#### Aliasing of memory and variable scope

```julia
a = [0, 1, 0]
b = a
b[2] = 0
??? a
```

```julia
function f(a)
	a = 3
end

a = 2
f(a)
?a
```
"""

# ╔═╡ a8e2d9a5-750f-4924-b245-5a4e8a5dddb0
md"""
### Attaching a backtrace to a log message
"""

# ╔═╡ df4f823a-da8b-444d-b04c-507d5c1bce4f
try
	error("Ye who enters, abandons all hope!")
catch err
	bt = catch_backtrace()
	@info "Caught error" exception=(err, bt)
end

# ╔═╡ 951590c4-88df-446b-b98f-e38de4bf5f31
md"""
This can be useful if you want to figure out in which context a log message was fired!
"""

# ╔═╡ 4ec0864a-480c-4b5b-939c-57672bc1608b
function fib(x)
	if x == 0 || x == 1
		@info "Basecase" exception=(ErrorException(""), backtrace())
		return x
	end
	fib(x-1)+fib(x-2)
end

# ╔═╡ 75d871c4-2999-4b10-934f-b613410e9083
fib(1)

# ╔═╡ 3a78addc-74c9-4ed3-a78e-71a68534e3ce
fib(3)

# ╔═╡ 0866e2b8-5425-4cb6-94dc-b99421981a14
md"""
## Callstacks and backtraces!

Profilers and debuggers think in backtraces! When we execute a program, each function needs a little bit of space on the stack.

!!! note
    Heap and Stack are regions in memory, the `heap` contains all **allocated** data and the stack contains the local variables & co of the currently executing functions.

given the execution of `fib(3)`

```
fib(3)
	fib(2)
       fib(1)
       fib(0)
       +
    fib(1)
    +
```

Each function may execute more functions, but we only execute **one** function at a time. So we visit the call-graph in order and thus linearize the execution.

```
fib(3)
```

```
fib(3)
fib(2)
```

```
fib(3)
fib(2)
fib(1)
```

```
fib(3)
fib(2)
fib(0)
```

```
fib(3)
fib(1)
```
"""

# ╔═╡ 9237038c-0114-4d32-ab09-c02890e49b1d
md"""
When we collect a profile, or receive an exception, or are in a debugger at a breakpoint, we almost always will collect the current backtrace/stacktrace.

Letting us know "how" did we get here.

Information that is not in a backtrace is as an example loops!
"""

# ╔═╡ 8a99b12b-2d4f-426b-8f9d-4d0d79ac59ee
md"""
## Julia-based debugging tools
"""

# ╔═╡ 87f77c8a-bdab-4616-8518-04dd4ee59881
md"""
### Debugger.jl

Julia's native debugger. Built upon [JuliaInterpreter.jl](https://github.com/JuliaDebug/JuliaInterpreter.jl)

#### Pro's

- Close to "your code"
- Full step-through debugging with variable inspection

#### Con's

- Only as fast as JuliaInterpreter and thus can struggle with large codebases

##### Example

```julia
using Debugger

function foo(n)
    x = n+1
    ((BigInt[1 1; 1 0])^x)[2,1]
end

@enter foo(20)
```

"""

# ╔═╡ d9119daf-cd85-4000-ac82-05616efe75f5
md"""
## Infiltrator.jl

Since Debugger.jl needs to interpret all code to ensure breakpoints are triggered it can be slow on a very large code to reach the region of interest!

```julia
function f(x)
	out = []
	for i in x
		push!(out, 2i)
		@infiltrate
	end
	out
end
```

##### Tricks
- Combine with Revise.jl and use `Main.@infiltrate` to work on package code

"""

# ╔═╡ f486d4f3-b819-4154-98f6-c7790aaaf027
md"""
## VS-Code

Debugger.jl is integrated into [VS-Code]( 
https://www.julia-vscode.org/docs/dev/userguide/debugging/)

it has similar caveats as before, but you may prefer it over the pure text interface.

"""

# ╔═╡ 79fee4d7-3102-4398-8ba3-be11a015f089
md"""
## System-native debuggers
"""

# ╔═╡ aa5c57a6-ac91-4a6d-8839-16aa19d40ce3
md"""
### GDB/LLDB

GDB and LLDB are native debugger, they are closest to the true execution of Julia.
They are very powerful but they require are deeper understanding on how programs work.

##### Pro's
- Faster than `Debugger.jl` since it debug's native code
- Faithful: Debugging the code doesn't change it's behavior

##### Con's

- Low-level view
- Doesn't know much about Julia code

##### Tricks
- `ENABLE_GDBLISTENER=1 julia -g2` (and potentially a debug build of Julia)
- https://docs.julialang.org/en/v1/devdocs/debuggingtips/

For a debugging build of Julia built with:
```
> cat Make.user
LLVM_ASSERTIONS=1
FORCE_ASSERTIONS=1
override CC=gcc-14
override CXX=g++-14

# default to a debug build for better line number reporting
override JULIA_BUILD_MODE=debug
```
"""

# ╔═╡ 33bd3052-65c5-461a-a255-6b436ab3f1f4
md"""
### rr

!!! note
    A key challenge when debugging is that we often ask the question "how did we get here", but debuggers are really good at answering the question "where are we going from here". 

[RR](https://rr-project.org/) is a time-traveling debugger. Using `rr record` we record the execution of a program, and it's children. Using `rr replay` we can re-execute the program faith-fully. We can then use `gdb` to debug a program being replayed and do things like `reverse-continue`, e.g. execute a program backwards in time to answer questions like: When did this value in memory get set?

##### Caveat
- Linux only
- Works most reliably on Intel CPU (boooh!)
- Doesn't work with GPU codes.

"""

# ╔═╡ 76947d0a-681a-4805-b5f5-bac32751693d
md"""
## Examples
"""

# ╔═╡ c51dc4c5-697a-4708-b2a9-3eea5cb9120a
md"""
!!! note
    `juliaup` can get in the way of debugging. Use `Base.julia_cmd()` to get the  actual binary.
"""

# ╔═╡ 942181ed-8214-4097-8f01-5c1afad59717
Base.julia_cmd()

# ╔═╡ 556a40c2-53b4-4b9a-925f-61544bef76a6
md"""
### Simple
"""

# ╔═╡ 16830265-6aad-4d35-a726-7291584ab35a
md"""
- `@ccall jl_breakpoint(val::Any)::Cvoid`
- `jl_(...)`
"""

# ╔═╡ 029eceeb-77ee-4fae-8591-a0780f01cef6
md"""
```sh
ENABLE_GDBLISTENER=1 gdb --args julia -g2 01_simple.jl
```
"""

# ╔═╡ a14cdb36-b5b9-439f-a73f-21ac40f53262
md"""
### PETSc (sigsegv)

User report:
```quote
The issue appears if I am in the REPL and run `include("script.jl")`, in which PETSc gets initialized. The script compiles correctly, but then I get SIGSEGV as soon as I type a single character in the REPL.

...

It also works fine if I set `--threads=1`
```

```sh
> julia -t4,1 --project=.
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.12.6 (2026-04-09)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org release
|__/                   |

julia> include("script.jl")
PETSc initialized: true
Num Julia threads: 4

julia> abfish: Job 1, 'julia -t4,1 --project=.' terminated by signal SIGSEGV (Address boundary error)
```
"""

# ╔═╡ 2c02e9dc-12f0-426c-af27-ea5fcacd0d29
md"""
```sh
/home/vchuravy/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia -g2 -t4,1 --project=. -e 'include("script.jl")' -i
```
"""

# ╔═╡ b9276851-fc1b-47ee-8165-5c36ebfc954d
md"""
- `-e` Run this code upon starting Julia
- `-i` enter interactive mode
- `-t` use multi-threading
"""

# ╔═╡ 5f175532-d452-4035-b4ce-397d50941194
md"""
### Polyester
"""

# ╔═╡ 0186ae3f-41d9-41b4-b48b-c6973e8eca6f
md"""
!!! note
    Bugs with parallel systems often require a specific sequence of scheduling decisions. One of the biggest challenging it to actual reproduce the reliably.

"""

# ╔═╡ 0ccccedd-5bcd-40f7-bf79-57991d25f93e
md"""
```sh
watch -e julia +1.12 -g2 --project=.  --threads=3 --check-bounds=yes --code-coverage=none test_threaded.jl
```
"""

# ╔═╡ 42d659ef-9828-4181-8677-b3134bc17ee7
md"""
#### `rr  record` 
"""

# ╔═╡ 9c770ea7-63d2-4e18-bb8b-68c9b6341585
md"""
```sh
ENABLE_GDBLISTENER=1 watch -e rr record -h ~/src/julia-1.12/julia -g2 --project=.  --threads=3 --check-bounds=yes --code-coverage=none test_threaded.jl
```
"""

# ╔═╡ b017a171-3a85-4dd1-98f0-c1c35e30fbe4
md"""
- `watch -e` run until the program fails
- `rr record -h` "chaos mode" allow for random perturbation of scheduling
"""

# ╔═╡ d41ce914-23e6-4d5b-81d6-a47aa41b73ba
md"""
#### `rr replay`
"""

# ╔═╡ 493d00bd-170f-4418-9c5c-143fc8e0f1a2
md"""
```
Thread 1 received signal SIGSEGV, Segmentation fault.
julia_wake_thread!_2999 (_tid=2) at /home/vchuravy/.julia/packages/ThreadingUtilities/gpppI/src/threadtasks.jl:61
61	@noinline function wake_thread!(_tid::T) where {T<:Integer}
(rr) bt 9
#0  julia_wake_thread!_2999 (_tid=2) at /home/vchuravy/.julia/packages/ThreadingUtilities/gpppI/src/threadtasks.jl:61
#1  0x00007bdd34a16d6d in launch () at /home/vchuravy/.julia/packages/ThreadingUtilities/gpppI/src/threadtasks.jl:18
#2  launch_batched_thread! () at /home/vchuravy/.julia/packages/Polyester/almvr/src/batch.jl:107
#3  macro expansion () at /home/vchuravy/.julia/packages/Polyester/almvr/src/batch.jl:241
#4  _batch_no_reserve () at /home/vchuravy/.julia/packages/Polyester/almvr/src/batch.jl:168
#5  batch () at /home/vchuravy/.julia/packages/Polyester/almvr/src/batch.jl:334
#6  macro expansion () at /home/vchuravy/.julia/packages/Polyester/almvr/src/closure.jl:456
#7  macro expansion () at /home/vchuravy/.julia/packages/Trixi/fS5XP/src/auxiliary/auxiliary.jl:234
#8  julia_calc_interface_flux!_21198 (backend=..., surface_flux_values=<optimized out>, mesh=..., have_nonconservative_terms=..., equations=...,
    surface_integral=..., dg=..., cache=...) at /home/vchuravy/.julia/packages/Trixi/fS5XP/src/solvers/dgsem_tree/dg_2d.jl:606
(More stack frames follow...)
```
"""

# ╔═╡ f6dde184-9ef1-44f0-82ca-3e20a6b926d4
md"""
### Address-sanitizer
"""

# ╔═╡ 393391b3-a844-49b1-bbb2-f7565d390db7
md"""
Address sanitizer is compiler based instrumentation to catch out of bounds writes. It is quite costly, but in particular in C based application can save your bacon.
"""

# ╔═╡ 5b1b0048-038a-41fe-8a0e-188f3a6c4bd8
md"""
```sh
> julia 04_asan.jl
> ~/src/julia-1.11-asan/asan/julia -g2 04_asan.jl
```
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Asciicast = "2600d445-abca-43b9-92aa-ce144ac0b05b"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Asciicast = "~0.1.3"
HypertextLiteral = "~0.9.5"
PlutoTeachingTools = "~0.4.1"
PlutoUI = "~0.7.65"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.6"
manifest_format = "2.0"
project_hash = "8f54979265e90e4087c86070a9628114c6c441e7"

[[deps.ANSIColoredPrinters]]
git-tree-sha1 = "574baf8110975760d391c710b6341da1afa48d8c"
uuid = "a4c015fc-c6ff-483c-b24f-f7ea428134e9"
version = "0.0.1"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "6c3913f4e9bdf6ba3c08041a446fb1332716cbc2"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Asciicast]]
deps = ["Base64", "Dates", "Documenter", "JSON3", "Logging", "MarkdownAST", "REPL", "Random", "StructTypes", "UUIDs", "agg_jll", "pandoc_jll"]
git-tree-sha1 = "03bcba2dc0d2b1ce97b30d99cc054668e5dc9703"
uuid = "2600d445-abca-43b9-92aa-ce144ac0b05b"
version = "0.1.3"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Documenter]]
deps = ["ANSIColoredPrinters", "AbstractTrees", "Base64", "CodecZlib", "Dates", "DocStringExtensions", "Downloads", "Git", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "MarkdownAST", "Pkg", "PrecompileTools", "REPL", "RegistryInstances", "SHA", "TOML", "Test", "Unicode"]
git-tree-sha1 = "56e9c37b5e7c3b4f080ab1da18d72d5c290e184a"
uuid = "e30172f5-a6a5-5a46-863b-614d45cd2de4"
version = "1.17.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8f05e9a2e7c2e3eb524102bb2926c5743c07fbe1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.8.0+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Git]]
deps = ["Git_LFS_jll", "Git_jll", "JLLWrappers", "OpenSSH_jll"]
git-tree-sha1 = "824a1890086880696fc908fe12a17bcf61738bd8"
uuid = "d7ba0133-e1db-5d97-8f8c-041e4b3a1eb2"
version = "1.5.0"

[[deps.Git_LFS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8c66e385d631bb934ff05e76d4a566c640c8df69"
uuid = "020c3dae-16b3-5ae5-87b3-4cb189e250b2"
version = "3.7.1+0"

[[deps.Git_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "LibCURL_jll", "Libdl", "Libiconv_jll", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "0dd4cfb426924210c8f42742751cbde74b27bfa3"
uuid = "f8c6e375-362e-5223-8a59-34ff63f689eb"
version = "2.54.0+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "411eccfe8aba0814ffa0fdf4860913ed09c34975"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.14.3"

    [deps.JSON3.extensions]
    JSON3ArrowExt = ["ArrowTypes"]

    [deps.JSON3.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c0c9b76f3520863909825cbecdef58cd63de705a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.5+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazilyInitializedFields]]
git-tree-sha1 = "0f2da712350b020bc3957f269c9caad516383ee0"
uuid = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
version = "1.3.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MarkdownAST]]
deps = ["AbstractTrees", "Markdown"]
git-tree-sha1 = "93c718d892e73931841089cdc0e982d6dd9cc87b"
uuid = "d0879d2d-cac2-40c8-9cee-1863dc0c7391"
version = "0.1.3"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenSSH_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenSSL_jll", "Zlib_jll"]
git-tree-sha1 = "57baa4b81a24c2910afbb6d853aa0685e4312bf7"
uuid = "9bd350c2-7e96-507f-8002-3f2e150b4e1b"
version = "10.3.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "5d5e0a78e971354b1c7bff0655d11fdc1b0e12c8"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.4"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "537c439831c0f8d37265efe850ee5c0d9c7efbe4"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3151a0c8061cc3f887019beebf359e6c4b3daa08"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.65"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegistryInstances]]
deps = ["LazilyInitializedFields", "Pkg", "TOML", "Tar"]
git-tree-sha1 = "ffd19052caf598b8653b99404058fce14828be51"
uuid = "2792f1a3-b283-48e8-9a74-f99dce5104f3"
version = "0.1.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "159331b30e94d7b11379037feeb9b690950cace8"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.11.0"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.agg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5c0aa73c27ac8d9a416384cace91edee36cccad1"
uuid = "5233c530-8ff5-54e9-96dc-9768cc911b78"
version = "1.4.3+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.pandoc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "47fe6932c1d26391bb6b49b146b20f2e543d1433"
uuid = "c5432543-76ad-5c9d-82bf-db097047a5e2"
version = "3.6.4+0"
"""

# ╔═╡ Cell order:
# ╠═5f6deede-7e22-4ebf-ae1a-14d584595f17
# ╠═aa46fef4-0c92-4947-ac64-f06ee31cb43f
# ╠═bc354629-dd2a-4b5a-8e9c-378dcc357575
# ╠═9182a8ea-7b06-4e4f-b0a3-35aa4cbfa456
# ╟─d5799787-aa74-4966-b647-bab4fe1d0c3c
# ╟─0be73b29-7780-4be0-bf07-9b62c99fc4b4
# ╟─1109474f-1cf7-4293-b952-d80d15273de5
# ╟─d5e479ec-7234-415b-845d-07f29f68eeaf
# ╟─c73a1c40-6426-4ced-a70b-909e54be4ced
# ╠═ef93ae34-fda8-488e-945b-6a3120de93c1
# ╠═8984c5c1-11a4-40c0-986d-2ed7db543098
# ╟─21185af9-6e23-4881-ba90-0fdd7af16471
# ╟─a2562dc3-ace7-47f4-a121-f83b500e5af0
# ╟─a8e2d9a5-750f-4924-b245-5a4e8a5dddb0
# ╠═df4f823a-da8b-444d-b04c-507d5c1bce4f
# ╟─951590c4-88df-446b-b98f-e38de4bf5f31
# ╠═4ec0864a-480c-4b5b-939c-57672bc1608b
# ╠═75d871c4-2999-4b10-934f-b613410e9083
# ╠═3a78addc-74c9-4ed3-a78e-71a68534e3ce
# ╟─0866e2b8-5425-4cb6-94dc-b99421981a14
# ╟─9237038c-0114-4d32-ab09-c02890e49b1d
# ╟─8a99b12b-2d4f-426b-8f9d-4d0d79ac59ee
# ╟─87f77c8a-bdab-4616-8518-04dd4ee59881
# ╟─d9119daf-cd85-4000-ac82-05616efe75f5
# ╟─f486d4f3-b819-4154-98f6-c7790aaaf027
# ╟─79fee4d7-3102-4398-8ba3-be11a015f089
# ╟─aa5c57a6-ac91-4a6d-8839-16aa19d40ce3
# ╟─33bd3052-65c5-461a-a255-6b436ab3f1f4
# ╟─76947d0a-681a-4805-b5f5-bac32751693d
# ╟─c51dc4c5-697a-4708-b2a9-3eea5cb9120a
# ╠═942181ed-8214-4097-8f01-5c1afad59717
# ╟─556a40c2-53b4-4b9a-925f-61544bef76a6
# ╟─16830265-6aad-4d35-a726-7291584ab35a
# ╟─029eceeb-77ee-4fae-8591-a0780f01cef6
# ╟─a14cdb36-b5b9-439f-a73f-21ac40f53262
# ╟─2c02e9dc-12f0-426c-af27-ea5fcacd0d29
# ╟─b9276851-fc1b-47ee-8165-5c36ebfc954d
# ╟─5f175532-d452-4035-b4ce-397d50941194
# ╟─0186ae3f-41d9-41b4-b48b-c6973e8eca6f
# ╟─0ccccedd-5bcd-40f7-bf79-57991d25f93e
# ╟─42d659ef-9828-4181-8677-b3134bc17ee7
# ╟─9c770ea7-63d2-4e18-bb8b-68c9b6341585
# ╟─b017a171-3a85-4dd1-98f0-c1c35e30fbe4
# ╟─d41ce914-23e6-4d5b-81d6-a47aa41b73ba
# ╟─493d00bd-170f-4418-9c5c-143fc8e0f1a2
# ╟─f6dde184-9ef1-44f0-82ca-3e20a6b926d4
# ╟─393391b3-a844-49b1-bbb2-f7565d390db7
# ╟─5b1b0048-038a-41fe-8a0e-188f3a6c4bd8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
