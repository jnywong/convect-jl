module myroutines
    using LinearAlgebra
    export poisson
    
    function poisson(b,nz,dz,n,a)
        d = fill(0.0,nz) # diagonal
        du = fill(-1/dz^2, nz-1) # upper diagonal
        dl = fill(-1/dz^2, nz-1) # lower diagonal
        for k=2:1:nz-1
            d[k] = (n*pi/a)^2 + 2/dz^2
        end

        # Boundary conditions
        d[1] = 1 # (2.22a)
        du[1] = 0 #(2.22b)
        b[1] = 0 # (2.22c)
        dl[nz-1] = 0 # (2.22d)
        d[nz] = 1 # (2.22e)
        b[nz] = 0 # (2.22f)

        A = Tridiagonal(dl, d, du)
        x = A\b

        return x
    end
end