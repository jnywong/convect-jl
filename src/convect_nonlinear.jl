include("routines.jl")
include("data_utils.jl")
using .routines
using .data_utils
using BenchmarkTools
using Plots
using LaTeXStrings

# TODO: odd modes are not updated

# Inputs
nz = 101 # no. of vertical gridpoints
nn = 20 # no. of Fourier modes (excluding zeroth mode)
a = 3 # L/D aspect ratio
Ra = 1e6 # Rayleigh number
Pr = 0.5 # Prandtl number
dt = 3e-6 # timestep size
nt = 100 # no. of timesteps
nout = 10 # save output every nout timesteps
initOn = 1 # initialise run, otherwise load existing data
saveDir = string("/Users/wongj/Documents/convect-out/","2021-08-19") # save directory

# Vertical domain
z, dz = routines.zdomain(nz)

# Pre-allocate arrays for spectral coefficients
psi, tem, omg = routines.preallocate_spec(nz,nn)

# Time integration
dtemdt, domgdt, tem, omg, psi = routines.nonlinear_solver(z, dz, nz, nn, nt, nout, dt, a, Ra, Pr, psi, tem, omg, initOn, saveDir)