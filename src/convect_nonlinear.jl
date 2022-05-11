include("routines.jl")
include("data_utils.jl")
using .routines
using .data_utils
using Dates
using DelimitedFiles
# using BenchmarkTools
# using Plots
# using LaTeXStrings
t_start = now()

@time begin
    # Inputs
    nz = 101 # no. of vertical gridpoints
    nn = 50 # no. of Fourier modes (excluding zeroth mode)
    a = 3 # L/D aspect ratio
    Ra = 1e6 # Rayleigh number
    Pr = 0.5 # Prandtl number
    dt = 3e-6 # timestep size
    nt = 1e4 # no. of timesteps
    nout = 1e2 # save output every nout timesteps
    initOn = 1 # initialise run, otherwise load existing data
    saveDir = string("/rds/projects/2017/edmondac-rescomp-software/wongj/convect-jl/data/",today()) # save directory

    # Vertical domain
    z, dz = routines.zdomain(nz)

    # Pre-allocate arrays for spectral coefficients
    psi, tem, omg = routines.preallocate_spec(nz,nn)

    # Time integration
    dtemdt, domgdt, tem, omg, psi = routines.nonlinear_solver(z, dz, nz, nn, nt, nout, dt, a, Ra, Pr, psi, tem, omg, initOn, saveDir)
end

writedlm(string(saveDir,"/log.csv"), string((now()-t_start).value/1000))