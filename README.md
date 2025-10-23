**blockSQP.jl** -- A Julia interface to blockSQP 2, a nonlinear  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; programming solver based on blockSQP by Dennis Janka.  
Copyright (c) 2025 Reinhold Wittmann <reinhold.wittmann@ovgu.de> and Christoph Plate <christoph.plate@ovgu.de>

## Installation

This package requires the compiled binary of the C interface to blockSQP_2 to be placed in the bin folder. The preferred way to obtain it is to fetch this package as a submodule for blockSQP 2 and have the build process generate and place the binary. In the future, a blockSQP_2_jll binary compatible with this version may be provided.  
This binary may link to, or include, compiled code subject to distinct licenses, see blockSQP_2/README.md for a list of dependencies and their licenses.

Afterwards, the package can be managed through the Julia package manager as usual, see <https://docs.julialang.org/en/v1/stdlib/Pkg/>.

### Licensing
blockSQP.jl is published under the very permissive zlib license, see LICENSE.txt.
