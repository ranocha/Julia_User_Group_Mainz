using PETSc

petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)

println("PETSc initialized: ", PETSc.initialized(petsclib))
println("Num Julia threads: ", Threads.nthreads())

