# Julia Workflows and helpful VS Code Shortcuts

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

## Helpful VS Code shortcuts

- Open **Zen Mode** using (Ctrl + K) + Z
- **Move Code lines** with Alt + Up/Down
- You can execute code even when commented out by selecting it explicitly
- You can **rename variables** using `F2` or right-click and "Rename Symbol" (but I will like this is not 100% reliable) [more info](https://www.julia-vscode.org/docs/stable/userguide/editingcode/#Rename-symbol)  
- You can **display arrays interactively** using `vscodedisplay()` [more info](https://www.julia-vscode.org/docs/stable/userguide/grid/#Table-Viewer)
- You can **execute a full code block** (delimited by `#-`) using Alt + Enter [more info](https://www.julia-vscode.org/docs/stable/userguide/runningcode/#Julia:-Execute-Code-Cell-in-REPL)
- To **format your code**, you can press Ctrl + Shift + P and search for `format`.
- The **search function** (for searching files) Ctrl + P is **really powerful**.
- searching for functions and variables use Ctrl + T or Ctrl + P and start with #. [more info](https://www.julia-vscode.org/docs/stable/userguide/codenavigation/#Open-Symbol-by-Name)  
- By right clicking on VS Code in the task bar, you can **open a whole new instantiation of VS Code.**

## Mathpix

Is a free app which allows you to get Latex Code of pictures/screenshots.

Copy pasting the Latex Code of a equation as a command into your .jl file, will make **Copilot be extremely powerful when trying to implement the equation.**

## Set default values for Plots.jl:

```julia
default(
    grid=true,
    box=:on,
    size=(700, 500),
    dpi=300,
    titlefont=font(16),
    linewidth=3, gridlinewidth=2,
    markersize=8, markerstrokewidth=4,
    xtickfontsize=14, ytickfontsize=14,
    xguidefontsize=16, yguidefontsize=16,
    legendfontsize=14
)
```
