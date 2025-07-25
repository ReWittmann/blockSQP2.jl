using Pkg;
# Activate current blockSQP environment
Pkg.activate(@__DIR__)
# Add binary
Pkg.add(url="https://github.com/chplate/qpoases_mumps.git", rev="main")
Pkg.add(url="https://github.com/chplate/blocksqp_mumps.git", rev="main")
