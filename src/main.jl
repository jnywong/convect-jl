include("myroutines.jl")
using .myroutines
using BenchmarkTools
using Plots

# Inputs
nz = 5 # no. of vertical gridpoints
nn = 3 # no. of Fourier modes 
a = 1 # L/D aspect ratio
Ra = 10 # Rayleigh number
Pr = 1 # Prandtl number
nt = 10 # no. of timesteps

# Domain
z = range(0,1,length=nz)
dz = 1/(nz-1)
dt = 0.9*dz^2/4 # timestep size (2.19)

# Pre-allocation
psi = Array{Float64}(undef, (nz,nn))
global omg = Array{Float64}(undef, (nz,nn)) # FIX GLOBALS: avoid ambiguous scoping by wrapping everything in functions
global tem = Array{Float64}(undef, (nz,nn))
dtemdz2 = Array{Float64}(undef, (nz,nn))
domgdz2 = Array{Float64}(undef, (nz,nn))
dtemdt = Array{Float64}(undef, (nz,nn,2))
domgdt = Array{Float64}(undef, (nz,nn,2))

# Initial conditions
tem0 = sin.(pi*z) # satisfies (3.6) T(0)=T(1) at time t=0
global tem = repeat(tem0, outer = [1,nn]) # repeat in x-direction

# Time integration
global m = 0
while m<=nt
    for k=2:1:nz-1
        for n=1:1:nn
            # Second derivatives in tem and omg
            dtemdz2[k,n] = (tem[k+1,n] - 2*tem[k,n] + tem[k-1,n])/dz^2 # (2.16)
            domgdz2[k,n] = (omg[k+1,n] - 2*omg[k,n] + omg[k-1,n])/dz^2 # (2.16)
            # Update time dtemdt and domgdt
            dtemdt[k,n,2] = (n*pi/a)*psi[k,n] + (dtemdz2[k,n]-(n*pi/a)^2*tem[k,n]) # (3.3)
            domgdt[k,n,2] = Ra*Pr*(n*pi/a)*tem[k,n] + Pr*(domgdz2[k,n] - (n*pi/a)^2*omg[k,n]) # (3.4)
        end
    end

    # Update tem and omg using Adams Bashforth time integration
    global tem = tem + 0.5*dt*(3*dtemdt[:,:,2]-dtemdt[:,:,1]) # (2.18)
    global omg = omg + 0.5*dt*(3*domgdt[:,:,2]-domgdt[:,:,1]) # (2.18)
    # Update psi using poisson solver
    for n=1:1:nn
        psi[:,n] = myroutines.poisson(omg[:,n],nz,dz,n,a) # (3.5)
    end
    dtemdt[:,:,1] = dtemdt[:,:,2]
    domgdt[:,:,1] = domgdt[:,:,2]
    global m+=1
end

# Plot
p = Plots.plot(psi[3,:])

