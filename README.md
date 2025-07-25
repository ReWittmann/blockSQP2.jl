 # Installation guide

 This package requires binaries for both the QP solver qpOASES in its sparse version, and blockSQP. Currently, the Julia packages containing these binaries are not registered. Therefore, both have to be added manually. The procedure for installing this package consists of two steps: 1) Cloning this repository, and 2) installing the packages with the binaries of qpOASES and blockSQP. 
 
 The commands are:

 ```
git clone https://kosinus.math.uni-magdeburg.de/mathopt/software/blocksqp.jl
```
or 
```
git clone git@kosinus.math.uni-magdeburg.de:mathopt/software/blocksqp.jl.git
```

And after that:
```
cd blocksqp.jl
julia install_deps.jl
```
This will install both packages into the blockSQP environment.