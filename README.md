 # Installation guide

 This package requires the QP solver qpOASES in its sparse version. Currently, the Julia package for this qpOASES version is not registered. Therefore, both .jll packages for qpOASES and blockSQP have to be added by hand. The procedure of installing this package is hence:

 ```
git clone https://kosinus.math.uni-magdeburg.de/mathopt/software/blocksqp.jl

cd blocksqp.jl
julia --project
```
Then you press `]` to enter the REPL mode in Julia. Here, you can then add the two packages from Github
```julia
add https://github.com/chplate/qpoases_mumps https://github.com/chplate/blocksqp_mumps
using blockSQP
 ```