include("routines.jl")
include("data_utils.jl")
using .routines
using .data_utils
using BenchmarkTools
using Plots
using LaTeXStrings

# Inputs
nz = 101 # no. of vertical gridpoints
nn = 50 # no. of Fourier modes (excluding zeroth mode)
a = 3 # L/D aspect ratio
Ra = 1e6 # Rayleigh number
Pr = 0.5 # Prandtl number
dt = 3e-6 # timestep size
nt = 1e6 # no. of timesteps
nout = 5e2 # save output every nout timesteps
initOn = 1 # initialise run, otherwise load existing data
saveDir = "/Users/wongj/Documents/convect-out/test" # save directory
zeroth = 1 # include zeroth order temperature in plot?

# Vertical domain
z, dz = routines.zdomain(nz)

# Pre-allocate arrays for spectral coefficients
psi, tem, omg = routines.preallocate_spec(nz,nn)

# Time integration
tem, omg, psi = routines.nonlinear_solver(z, dz, nz, nn, nt, nout, dt, a, Ra, Pr, psi, tem, omg, initOn, saveDir)

# Horizontal domain
# x, dx, nx = routines.xdomain(a,nz)

# # Spectral to spatial transforms in x-direction
# tem_full, omg_full, psi_full = routines.preallocate_spat(nz,nx)
# cosa = routines.cosines(a,x,nn,nx) # compute and store cosines
# sina = routines.sines(a,x,nn,nx) # compute and store sines
# tem_full = routines.ict(nn,nx,nz,cosa,tem,tem_full,zeroth) # inverse cosine transform
# omg_full = routines.ist(nn,nx,nz,sina,omg,omg_full) # inverse sine transform
# psi_full = routines.ist(nn,nx,nz,sina,psi,psi_full) 

# # Plot
# psiLim = maximum(abs.(psi_full))
# temLim = maximum(abs.(tem_full))
# if max(1,temLim)==1
#     cbar_range=(0,1)
# else
#     cbar_range=(-temLim, temLim)
# end
# psi_scale = psi_full/psiLim*cbar_range[2] # scale so that psi is within cbar range for plot

# p1 = heatmap(x, z, tem_full, c=:balance, levels = 100, clim=cbar_range, colorbar=true)
# contour!(x, z, psi_scale, c=:grays,levels=20,lw=2,colorbar = true)
# xlabel!(L"x")
# ylabel!(L"z")
# a_title = round(a,digits=2)
# title!(latexstring("\$Ra = $(Ra), Pr = $(Pr), a = $(a_title), N_z = $(nz), N_n = $(nn)\$"))
# plot(p1)

# psi, omg, tem = data_utils.load_data(saveDir)