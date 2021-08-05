module myroutines
    using LinearAlgebra
    export poisson
    
    function poisson(b,nz,dz,n,a)
        d = Array{Float64}(undef, nz) # diagonal
        du = fill(-1/dz^2, nz-1) # upper diagonal
        dl = fill(-1/dz^2, nz-1) # lower diagonal
        for k=2:1:nz-1
            d[k] = (n*pi/a)^2 + 2/dz^2
        end

        # Boundary conditions
        d[1] = 1
        du[1] = 0
        d[nz] = 1
        dl[nz-1] = 0
        b[1] = 0
        b[nz] = 0

        A = Tridiagonal(dl, d, du)
        x = A\b

        return x
    end
end