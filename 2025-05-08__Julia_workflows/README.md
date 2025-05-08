# Julia Workflows

These are some notes on possible Julia workflows. A nice ressource is
https://modernjuliaworkflows.org/


## Overview

I typically use a few different workflows based on the tasks

- Pluto.jl notebooks for teaching
- Code development while working on a paper
- Package development


For papers, I typically use a setup like

```
.
├── .gitignore
├── code
│   ├── code.jl
│   ├── Manifest.toml
│   ├── Project.toml
│   └── README.md
├── paper.tex
└── references.bib
```

I typically use [Visual Studio Code](https://code.visualstudio.com/)
to develop the code and write the paper (using the LaTeX and Julia
extensions). In addition, I often use a standard terminal next to
VSCode where I start a REPL session as

```julia
julia> using Revise; includet("code/code.jl")
```

Then, I can modify the source code and use the latest version in the
Julia REPL. The `code.jl` file typically begins with something like

```julia
# Install packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Load packages
using Trixi
```

When working on packages, I prefer a similar workflow - VSCode plus
terminal. Inside the package directory, I create a `run` directory.
While working on the package, I start Julia with the `run` directory
as project where I `Pkg.develop` the package I am working on. This
allows me to `Pkg.add` additional packages I use for development.
