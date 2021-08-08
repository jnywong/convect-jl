include("routines.jl")
using .routines
using BenchmarkTools
using Plots
using LaTeXStrings

# Inputs
nz = 101 # no. of vertical gridpoints
nn = 3 # no. of Fourier modes (excluding 0)
a = sqrt(2) # L/D aspect ratio
Ra = 1500 # Rayleigh number
Pr = 0.3 # Prandtl number
nt = 10^5 # no. of timesteps
nout = 500 # output every nout timesteps

# Vertical domain
z, dz = routines.zdomain(nz)

# Timestep size
dt = 0.9*dz^2/4 # (2.19)

# Pre-allocation
psi, tem, omg, dtemdz2, domgdz2, dtemdt, domgdt = routines.preallocate_spec(nz,nn)

# Initial conditions
for n=1:1:nn
    tem[:,n,2] = routines.initial_tem(z)
end

# Time integration
tem, omg, psi = routines.sol(z, dz, nz, nn, nt, nout, dt, a, Ra, Pr, psi, tem, omg, dtemdz2, domgdz2, dtemdt, domgdt)

# horizontal domain
x, dx = routines.xdomain(a,nz)

# Spectral to spatial transforms in x-direction
tem_full, omg_full, psi_full = routines.preallocate_spat(nz,nx)
tem_full = routines.ref_tem(nx,nz,z,tem_full) # reference temperature
tem_full = routines.ict(a,x,nn,nx,nz,tem,tem_full) # inverse cosine transform
psi_full = routines.ist(a,x,nn,nx,nz,psi,psi_full) # inverse sine transform

# Plot
temLim = maximum(abs.(tem_full))
psiLim = maximum(abs.(psi_full))
p1 = heatmap(x, z, tem_full, c=:balance, levels = 100) #, clim=(-temLim,temLim))
contour!(x, z, psi_full, c=:grays,levels=20,lw=2,colorbar=true) #, clim=(-psiLim,psiLim))
xlabel!(L"x")
ylabel!(L"z")
a_title = round(a,digits=2)
title!(latexstring("\$Ra = $(Ra), Pr = $(Pr), a = $(a_title), N_z = $(nz), N_n = $(nn)\$"))
plot(p1)

# savefig("docs/figures/linear.png")