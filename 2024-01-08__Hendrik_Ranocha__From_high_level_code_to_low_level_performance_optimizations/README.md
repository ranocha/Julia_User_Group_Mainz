# From high-level code to low-level performance optimizations

I will introduce briefly how I use Julia in my daily work.
Afterwards, I will give a high-level overview of performance
aspects in Julia, ranging from basic concepts such as types
and type stability to low-level performance optimizations.

Run

```bash
julia --project=. -e 'import Pkg; Pkg.instantiate(); import Pluto; Pluto.run()'
```

in this directory. Then, open the Pluto notebook `talk.jl`.

This code was developed with Julia v1.10.0.
