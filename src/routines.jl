module routines
    using Printf
    using LinearAlgebra
    export all
    
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

    function xdomain(a, nz)
        nx = Int(ceil(a*nz)) # no. of horizontal gridpoints
        dx = a/(nx-1)
        x = LinRange(0.0,a,nx)
        return x, dx
    end

    function zdomain(nz)
        z = LinRange(0,1,nz)
        dz = 1/(nz-1)
        return z, dz
    end

    function preallocate_spat(nz, nx)
        tem = zeros(nz, nx)
        psi = zeros(nz, nx)
        omg = zeros(nz, nx)
        return tem, omg, psi
    end

    function preallocate_spec(nz, nn)
        psi = zeros(nz,nn,2)
        tem = zeros(nz,nn,2)
        omg = zeros(nz,nn,2)
        dtemdz2 = zeros(nz,nn)
        domgdz2 = zeros(nz,nn)
        dtemdt = zeros(nz,nn,2)
        domgdt = zeros(nz,nn,2)
        return psi, tem, omg, dtemdz2, domgdz2, dtemdt, domgdt
    end

    function initial_tem(z)
        y = sin.(pi*z) # satisfies (3.6) T(0)=T(1) at time t=0
        return y
    end

    function ref_tem(nx,nz,z,tem_out)
        for i=1:1:nx
            tem_out[:,i] = fill(1,nz) - z # reference temperature
        end
        return tem_out
    end

    function first_deriv(k,n,dz,tem,omg,dtemdz2,domgdz2)
        # First z derivatives in tem and omg
        dtemdz1[k,n] = (tem[k+1,n,2] - tem[k-1,n,2])/(2*dz) # (2.15)
        domgdz1[k,n] = (omg[k+1,n,2] - omg[k-1,n,2])/(2*dz) # (2.15)
        return dtemdz1, domgdz1
    end

    function second_deriv(k,n,dz,tem,omg,dtemdz2,domgdz2)
        # Second z derivatives in tem and omg
        dtemdz2[k,n] = (tem[k+1,n,2] - 2*tem[k,n,2] + tem[k-1,n,2])/dz^2 # (2.16)
        domgdz2[k,n] = (omg[k+1,n,2] - 2*omg[k,n,2] + omg[k-1,n,2])/dz^2 # (2.16)
        return dtemdz2, domgdz2
    end

    function tem_eq(k,n,a,tem,psi,dtemdz2,dtemdt)
        # Update time deriv dtemdt
        dtemdt[k,n,2] = (n*pi/a)*psi[k,n,2] + (dtemdz2[k,n]-(n*pi/a)^2*tem[k,n,2]) # (3.3)
        return dtemdt
    end

    function omg_eq(k,n,a,Ra,Pr,tem,omg,domgdz2,domgdt)
        # Update time deriv domgt
        domgdt[k,n,2] = Ra*Pr*(n*pi/a)*tem[k,n,2] + Pr*(domgdz2[k,n] - (n*pi/a)^2*omg[k,n,2]) # (3.4)
        return domgdt 
    end

    function adamsbashforth(tem, omg, dt, dtemdt, domgdt)
        tem[:,:,2] = tem[:,:,1] + 0.5*dt*(3*dtemdt[:,:,2]-dtemdt[:,:,1]) # (2.18)
        omg[:,:,2] = omg[:,:,1] + 0.5*dt*(3*domgdt[:,:,2]-domgdt[:,:,1]) # (2.18)
        return tem, omg
    end

    function diagnostics(m,nz,nout,time,tem,omg,psi)
        if mod(m,nout)==0
            @printf("time: %.2f    tem: %.5e    omg: %.5e    psi: %.5e \n ", time, log(abs(tem[round.(Int,nz/3),1,2]))-log(abs(tem[round.(Int,nz/3),1,1])),log(abs(omg[round.(Int,nz/3),1,2]))-log(abs(omg[round.(Int,nz/3),1,1])),log(abs(psi[round.(Int,nz/3),1,2]))-log(abs(psi[round.(Int,nz/3),1,1])))
        end
    end

    function prepare(dtemdt, domgdt, tem, omg, psi)
        dtemdt[:,:,1] = dtemdt[:,:,2]
        domgdt[:,:,1] = domgdt[:,:,2]
        tem[:,:,1] = tem[:,:,2]
        omg[:,:,1] = omg[:,:,2]
        psi[:,:,1] = psi[:,:,2]
        return dtemdt, domgdt, tem, omg, psi
    end

    function linear_solver(z, dz, nz, nn, nt, nout, dt, a, Ra, Pr, psi, tem, omg, dtemdz2, domgdz2, dtemdt, domgdt)
        m = 0
        time = 0
        while m<=nt
            for k=2:1:nz-1
                for n=1:1:nn
                    dtemdz2, domgdz2 = second_deriv(k ,n, dz, tem, omg, dtemdz2, domgdz2)
                    dtemdt = tem_eq(k,n,a,tem,psi,dtemdz2,dtemdt)
                    domgdt = omg_eq(k,n,a,Ra,Pr,tem,omg,domgdz2,domgdt)
                end
            end

            # Update tem and omg using Adams Bashforth time integration
            tem, omg = adamsbashforth(tem, omg, dt, dtemdt, domgdt)

            # Update psi using poisson solver
            for n=1:1:nn
                psi[:,n,2] = poisson(omg[:,n,2],nz,dz,n,a) # (3.5)
            end

            # Diagnostics
            diagnostics(m, nz, nout, time, tem, omg, psi)

            # Prepare values for next timestep
            dtemdt, domgdt, tem, omg, psi = prepare(dtemdt, domgdt, tem, omg, psi)

            m+=1
            time += dt
        end
        return tem, omg, psi
    end

    function cosines(a, x, nn, nx)
        # Compute cosines and sines 
        cosa = zeros(nn,nx)
        for n=1:1:nn
            cosa[n,:] = cos.(n*pi*x/a)
        end
        return cosa
    end

    function sines(a, x, nn, nx)
        # Compute cosines and sines 
        sina = zeros(nn,nx)
        for n=1:1:nn
            sina[n,:] = sin.(n*pi*x/a)
        end
        return sina
    end    

    function ict(nn,nx,nz,cosa,coeffs,outfun)
        for i = 1:1:nx
            for k = 1:1:nz
                for n=1:1:nn
                    outfun[k,i] += coeffs[k,n,2]*cosa[n,i]
                end
            end
        end
        return outfun
    end

    function ist(nn,nx,nz,sina,coeffs,outfun)
        for i = 1:1:nx
            for k = 1:1:nz
                for n=1:1:nn
                    outfun[k,i] += coeffs[k,n,2]*sina[n,i]
                end
            end
        end
        return outfun
    end

end



