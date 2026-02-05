using Pkg

Pkg.activate(@__DIR__)
Pkg.develop(path = joinpath(@__DIR__, ".."))
Pkg.develop(path = joinpath(@__DIR__, "../lib/NLPstructures"))
Pkg.add("OptimizationBase")
#Pkg.add("OptimizationMOI")
Pkg.add("Ipopt")

Pkg.add("QPALM")

Pkg.add("DifferentiationInterface")
Pkg.add("SparseConnectivityTracer")
Pkg.add("SparseMatrixColorings")
Pkg.add("SparseArrays")
Pkg.add("ForwardDiff")
Pkg.add("Zygote")

Pkg.add("CairoMakie")
