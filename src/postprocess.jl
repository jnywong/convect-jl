include("routines.jl")
include("data_utils.jl")
using .routines
using .data_utils
using BenchmarkTools
using Plots
using LaTeXStrings

# Inputs
saveDir = "/Users/wongj/Documents/convect-out/2021-08-19" # save directory
ndata = 0
zeroth=0

# Load
nz,nn,a,Ra,Pr,dt,nt,nout = data_utils.load_inputs(saveDir)
dtemdt, domgdt, tem, omg, psi = data_utils.load_data(saveDir, ndata)

# Domain
z, dz = routines.zdomain(nz)
x, dx, nx = routines.xdomain(a,nz)

# Spectral to spatial transforms in x-direction
tem_full, omg_full, psi_full = routines.preallocate_spat(nz,nx)
cosa = routines.cosines(a,x,nn,nx) # compute and store cosines
sina = routines.sines(a,x,nn,nx) # compute and store sines
tem_full = routines.ict(nn,nx,nz,cosa,tem,tem_full,zeroth) # inverse cosine transform
omg_full = routines.ist(nn,nx,nz,sina,omg,omg_full) # inverse sine transform
psi_full = routines.ist(nn,nx,nz,sina,psi,psi_full) 

# Plot
psiLim = maximum(abs.(psi_full))
temLim = maximum(abs.(tem_full))
if max(1,temLim)==1
    cbar_range=(0,1)
else
    cbar_range=(-temLim, temLim)
end
# cbar_range=(-temLim, temLim)
psi_scale = psi_full/psiLim*cbar_range[2] # scale so that psi is within cbar range for plot

p1 = heatmap(x, z, tem_full, c=:balance, levels = 100, clim=cbar_range, colorbar=true)
contour!(x, z, psi_scale, c=:grays,levels=20,lw=2,colorbar = true)
xlabel!(L"x")
ylabel!(L"z")
a_title = round(a,digits=2)
title!(latexstring("\$Ra = $(Ra), Pr = $(Pr), a = $(a_title), N_z = $(nz), N_n = $(nn)\$"))
plot(p1)
