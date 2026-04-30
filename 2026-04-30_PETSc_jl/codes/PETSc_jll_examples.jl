# This shows how to run some of the binaries of the PETSc_jll package
using PETSc_jll # download the PETSc binaries for your platform

# We can run example 19 with this
run(`$(PETSc_jll.ex19())`);

# If you want to see te convergence of the outer solver do this: 
run(`$(PETSc_jll.ex19()) -ksp_monitor`);

# help options 
run(`$(PETSc_jll.ex19()) -help`);

#Set up a multigrid solver:
run(`$(PETSc_jll.ex19()) -ksp_type gmres -pc_type mg -mg_levels 3 `);

#detailed timing of the full simulation:
run(`$(PETSc_jll.ex19()) -log_view`);

#Show convergence of nonlinear (outer=SNES) solver as well as that of the inner KSP solver:
run(`$(PETSc_jll.ex19()) -ksp_monitor -snes_monitor`);

