# Installation

To install the package, you can use the Julia package manager. Open the Julia REPL and run the following command:

```julia
using Pkg; Pkg.add("MyPkg")
```
or
```julia-repl
julia> ]

(@v1.11) pkg> add `PATH to your GitHub repo~
```

!!! info "Install from a specific branch"
    However, as the API is changing and not every new feature leads to a new release, one can also clone the main branch of the repository:
    ```julia
    add MyPkg#main
    ```

After installation, you can test the package by running the following commands:
```julia-repl
using MyPkg

julia> ]

(@v1.11) pkg> test MyPkg
```
The test will take a while, so grab a ☕️ or 🍵
