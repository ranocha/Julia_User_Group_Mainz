# Introduction to HOHQMesh

This code was developed with Julia v1.10.5. To run the code,
you need to download and install Julia first, e.g., from the official
[download website](https://julialang.org/downloads/).
This presentation uses a
[Pluto.jl](https://github.com/fonsp/Pluto.jl) notebook.


## Working from a terminal (e.g., Linux or macOS)

Open a terminal and run

```bash
julia -e 'import Pkg; Pkg.activate(pwd()); Pkg.instantiate(); import Pluto; Pluto.run()'
```

in this directory. Then, open the Pluto notebook `intro_hohqmesh.jl`.


## Windows

Open the Julia REPL in this directory, e.g., by navigating to this directory
in the Windows File Explorer and typing "julia" the address bar. A console
window should open. There, execute

```julia
import Pkg; Pkg.activate(pwd()); Pkg.instantiate(); import Pluto; Pluto.run()
```

in the Julia REPL. Then, open the Pluto notebook `intro_hohqmesh.jl`.


## Run the examples

In addition to the notebook, the repository also includes the example explored in the Pluto notebook. It is recommended to run the `von_karman.jl` test case with multiple threads, starting the julia via the following command:

```bash
julia --project=. --threads=8
```

The mesh file for the Von Karman Vortex Sheets test case has to be generated first.

```julia
include("cylinder_mesh.jl")
```

The mesh will be then generated in the directory 'out' and will be automatically loaded in the `von_karman.jl` test case.