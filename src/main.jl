using LinearAlgebra

# Poisson solver

# Inputs
z0 = 0
z1 = 1
nz = 5 # no. of vertical gridpoints
n = 0 # wavenumber (for horizontal decomposition)
a = 1 # L/D aspect ratio

# Domain
z = range(z0,z1,length=nz)
dz = z[2] - z[1]

# Tridiagonal matrix
d = Array{Float64}(undef, nz) # diagonal
du = fill(-1/dz^2, nz-1) # upper diagonal
dl = fill(-1/dz^2, nz-1) # lower diagonal
for i=2:1:nz-1
    d[i] = (n*pi/a)^2 + 2/dz^2
end

# RHS
b = rand(nz)

# Boundary conditions
d[1] = 1
du[1] = 0
d[nz] = 1
dl[nz-1] = 0
b[1] = 0
b[nz] = 0

A = Tridiagonal(dl, d, du)
x = A\b
b_chk = A*x
err_chk = sum(b-b_chk)
println("Sum of errors is $err_chk")

