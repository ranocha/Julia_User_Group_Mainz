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