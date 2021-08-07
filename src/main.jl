include("myroutines.jl")
using .myroutines
using FFTW
using BenchmarkTools
using Plots
using Printf

# Inputs
nz = 101 # no. of vertical gridpoints
nn = 3 # no. of Fourier modes (excluding 0)
a = sqrt(2) # L/D aspect ratio
Ra = 657.5 # Rayleigh number
Pr = 0.3 # Prandtl number
nt = 10^5 # no. of timesteps
nout = 500 # output every nout timesteps

# Domain
z = LinRange(0,1,nz)
dz = 1/(nz-1)
dt = 0.9*dz^2/4 # timestep size (2.19)
# dt = 3*10^(-6)

# Pre-allocation
psi = zeros(nz,nn,2)
tem = zeros(nz,nn,2)
omg = zeros(nz,nn,2)
dtemdz2 = zeros(nz,nn)
domgdz2 = zeros(nz,nn)
dtemdt = zeros(nz,nn,2)
domgdt = zeros(nz,nn,2)

# Initial conditions
for n=1:1:nn
    tem[:,n,2] = sin.(pi*z) # satisfies (3.6) T(0)=T(1) at time t=0
end

# Time integration
global m = 0
global time = 0
while m<=nt
    for k=2:1:nz-1
        for n=1:1:nn
            # Second derivatives in tem and omg
            dtemdz2[k,n] = (tem[k+1,n,2] - 2*tem[k,n,2] + tem[k-1,n,2])/dz^2 # (2.16)
            domgdz2[k,n] = (omg[k+1,n,2] - 2*omg[k,n,2] + omg[k-1,n,2])/dz^2 # (2.16)
            # Update time dtemdt and domgdt
            dtemdt[k,n,2] = (n*pi/a)*psi[k,n,2] + (dtemdz2[k,n]-(n*pi/a)^2*tem[k,n,2]) # (3.3)
            domgdt[k,n,2] = Ra*Pr*(n*pi/a)*tem[k,n,2] + Pr*(domgdz2[k,n] - (n*pi/a)^2*omg[k,n,2]) # (3.4)
        end
    end

    # Update tem and omg using Adams Bashforth time integration
    tem[:,:,2] = tem[:,:,1] + 0.5*dt*(3*dtemdt[:,:,2]-dtemdt[:,:,1]) # (2.18)
    omg[:,:,2] = omg[:,:,1] + 0.5*dt*(3*domgdt[:,:,2]-domgdt[:,:,1]) # (2.18)

    # Update psi using poisson solver
    for n=1:1:nn
        psi[:,n,2] = myroutines.poisson(omg[:,n,2],nz,dz,n,a) # (3.5)
    end

    # Diagnostics
    if mod(m,nout)==0
        @printf("time: %.2f    tem: %.5e    omg: %.5e    psi: %.5e \n ", time, log(abs(tem[round.(Int,nz/3),1,2]))-log(abs(tem[round.(Int,nz/3),1,1])),log(abs(omg[round.(Int,nz/3),1,2]))-log(abs(omg[round.(Int,nz/3),1,1])),log(abs(psi[round.(Int,nz/3),1,2]))-log(abs(psi[round.(Int,nz/3),1,1])))
    end

    # Store values for next timestep
    dtemdt[:,:,1] = dtemdt[:,:,2]
    domgdt[:,:,1] = domgdt[:,:,2]
    tem[:,:,1] = tem[:,:,2]
    omg[:,:,1] = omg[:,:,2]
    psi[:,:,1] = psi[:,:,2]

    global m+=1
    global time += dt

end

# test = FFTW.idct(tem,2)

# Spectral to spatial in x-direction
nx = Int(ceil(a*nz)) # no. of horizontal gridpoints
dx = a/(nx-1)
x = LinRange(0.0,a,nx)

# Store cosines and sines
cosa = zeros(nn,nx)
sina = zeros(nn,nx)
for n=1:1:nn
    cosa[n,:] = cos.(n*pi*x/a)
    sina[n,:] = sin.(n*pi*x/a)
end

tem_full = zeros(nz, nx)
psi_full = zeros(nz, nx)
# for i=1:1:nx
#     tem_full[:,i] = fill(1,nz) - z # reference temperature
# end
for i = 1:1:nx
    for k = 1:1:nz
        for n=1:1:nn
            tem_full[k,i] += tem[k,n,2]*cosa[n,i]
            psi_full[k,i] += psi[k,n,2]*sina[n,i]
        end
    end
end

# Plot
# p = Plots.plot(psi[3,:])
temLim = maximum(abs.(tem_full))
p1 = contourf(x, z, tem_full, c=:balance)
# contour!(x, z, psi_full, c=:grays,levels=50) #,clim=(-temLim,temLim))
p2 = contour(x, z, psi_full, c=:grays)
plot(p1,p2)