using Documenter

# If your package is not yet installed, you can use this line to load it
push!(LOAD_PATH,"../src/")
using MyPkg

@info "Making documentation..."

makedocs(
    sitename = "MyPkg",
    authors = "First Author etc",
    format = Documenter.HTML(;
    size_threshold_ignore = ["man/listfunctions.md"]), # easier local build
    modules = [MyPkg],

    pages = [
        "Home" => "index.md",
        "Installation" => "man/installation.md",
        "List of functions" => "man/listfunctions.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
