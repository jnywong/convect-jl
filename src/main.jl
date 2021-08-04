include("myroutines.jl")
using .myroutines
using BenchmarkTools

# Inputs
z0 = 0
z1 = 1
nz = 5 # no. of vertical gridpoints
n = 0 # wavenumber (for horizontal decomposition)
a = 1 # L/D aspect ratio

# Domain
z = range(z0,z1,length=nz)
dz = z[2] - z[1]
@btime x = myroutines.poisson(nz,dz,n,a)
