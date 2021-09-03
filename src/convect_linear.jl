include("routines.jl")
using .routines
using BenchmarkTools
using Plots
using LaTeXStrings

# NOTE: too much power in each mode compared with nonlinear case 

# Inputs
nz = 101 # no. of vertical gridpoints
nn = 30 # no. of Fourier modes (excluding zeroth mode)
a = 5 # L/D aspect ratio
Ra = 2700 # Rayleigh number
Pr = 0.5 # Prandtl number
nt = 1e5 # no. of timesteps
nout = 1e3 # output every nout timesteps
zeroth = 0 # include zeroth order temperature in plot?
initOn = 1
saveDir = "/Users/wongj/Documents/convect-out/linear/2021-09-02"

# Vertical domain
z, dz = routines.zdomain(nz)

# Timestep size
dt = 0.9*dz^2/4 # (2.19)

# Pre-allocate arrays for spectral coefficients
psi, tem, omg = routines.preallocate_spec(nz,nn)

# Initial conditions
tem = routines.initial_linear_tem(nz,nn,z,tem)

# Time integration
dtemdt, domgdt, tem, omg, psi = routines.linear_solver(z, dz, nz, nn, nt, nout, dt, a, Ra, Pr, psi, tem, omg, initOn, saveDir)

# horizontal domain
x, dx, nx = routines.xdomain(a,nz)

# Spectral to spatial transforms in x-direction
tem_full, omg_full, psi_full = routines.preallocate_spat(nz,nx)
cosa = routines.cosines(a,x,nn,nx) # compute and store cosines
sina = routines.sines(a,x,nn,nx) # compute and store sines
tem_full = routines.ict(nn,nx,nz,cosa,tem,tem_full,zeroth) # inverse cosine transform
psi_full = routines.ist(nn,nx,nz,sina,psi,psi_full) # inverse sine transform

# Plot
psiLim = maximum(abs.(psi_full))
temLim = maximum(abs.(tem_full))
if max(1,temLim)==1
    cbar_range=(0,1)
else
    cbar_range=(-temLim, temLim)
end
psi_scale = psi_full/psiLim*cbar_range[2] # scale so that psi is within cbar range for plot

p1 = heatmap(x, z, tem_full, c=:balance, levels = 100, clim=cbar_range, colorbar=true)
contour!(x, z, psi_scale, c=:grays,levels=20,lw=2,colorbar = true)
xlabel!(L"x")
ylabel!(L"z")
a_title = round(a,digits=2)
title!(latexstring("\$Ra = $(Ra), Pr = $(Pr), a = $(a_title), N_z = $(nz), N_n = $(nn)\$"))
plot(p1)

# savefig("docs/figures/linear.png")