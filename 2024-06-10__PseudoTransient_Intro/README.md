This repository was created to host some examples for the Pseudotransient tutorial on Monday 10.June.2024.

The main idea is to start with the code "elliptic1d.jl"
This code has been extensively documented and shows all the main points.
The main steps are explained in detail in the power-point pdf (although with MATLAB examples).
The two following versions show accelerated versions (2nd order iterations of the main code).
The code "elliptic1d_acc.jl" shows the implementation of the acceleration (see residual definition).
The code "elliptic1d_acc_nlin.jl" shows the treatment of non-linearities  (see effective Keff).

In addition, a script comparing different function definitions is given (requires BenchmarkTools),
apart from this and basic plotting tools no special packages are needed.

Evangelos Moulas, Mainz 10-June 2024
