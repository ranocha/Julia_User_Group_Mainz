# START OF PROLOGUE
if isfile("libbinary_playground.dylib")
    const lib_BinaryPlayground = joinpath(pwd(),"libbinary_playground.dylib")
    println("Using locally compiled version of libbinary_playground.dylib")
else
    warning("libbinary_playground.dylib not found in current directory.")
end