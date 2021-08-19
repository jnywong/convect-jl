module routines
    using Printf
    using LinearAlgebra
    using Distributions
    include("data_utils.jl")
    using .data_utils
    export all
    
    function poisson(b,nz,dz,n,a)
        d = fill(0.0,nz) # diagonal
        du = fill(-1/dz^2, nz-1) # upper diagonal
        dl = fill(-1/dz^2, nz-1) # lower diagonal
        for k=2:1:nz-1
            d[k] = ((n-1)*pi/a)^2 + 2/dz^2
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
        return x, dx, nx
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
        # Includes zeroth mode
        psi = zeros(nz,nn+1,2)
        tem = zeros(nz,nn+1,2)
        omg = zeros(nz,nn+1,2)
        return psi, tem, omg 
    end

    function initial_linear_tem(nz,nn,z,tem)
        tem[:,1,2] = fill(1,nz) - z # zeroth mode
        for n=2:1:nn+1
            tem[:,n,2] = sin.(pi*z) # satisfies (3.6) T(0)=T(1) at time t=0
        end
        return tem
    end

    function initial_nonlinear_tem(nz,nn,z,tem,fac=1e-3)
        tem[:,1,2] = fill(1,nz) - z # zeroth mode
        tem[:,2,2] = 0.01*sin.(pi*z) # §4.4.2 A nonlinear benchmark
        # for n=2:1:nn+1
        #     tem[:,n,2] = fac*rand(Uniform(-1, 1))*sin.(pi*z) # satisfies (3.6) T(0)=T(1) at time t=0
        # end
        return tem
    end

    function ref_tem(nx,nz,z,tem_out)
        for i=1:1:nx
            tem_out[:,i] = fill(1,nz) - z # reference temperature
        end
        return tem_out
    end

    function first_deriv(k,n,dz,y,dydz)
        # First derivative w.r.t z
        dydz[k,n] = (y[k+1,n,2] - y[k-1,n,2])/(2*dz) # (2.15)
        return dydz
    end

    function second_deriv(k,n,dz,y,dydz2)
        # Second derivative w.r.t z
        dydz2[k,n] = (y[k+1,n,2] - 2*y[k,n,2] + y[k-1,n,2])/dz^2 # (2.16)
        return dydz2
    end

    function lin_tem_eq(k,n,a,tem,psi,dtemdz2,dtemdt)
        # Update time deriv dtemdt
        dtemdt[k,n,2] = ((n-1)*pi/a)*psi[k,n,2] + (dtemdz2[k,n]-((n-1)*pi/a)^2*tem[k,n,2]) # (3.3)
        return dtemdt
    end

    function lin_omg_eq(k,n,a,Ra,Pr,tem,omg,domgdz2,domgdt)
        # Update time deriv domgt
        domgdt[k,n,2] = Ra*Pr*((n-1)*pi/a)*tem[k,n,2] + Pr*(domgdz2[k,n] - ((n-1)*pi/a)^2*omg[k,n,2]) # (3.4)
        return domgdt 
    end

    function nonlin_tem_eq(k,n,a,tem,psi,dtemdz2,dtemdt)
        # Update time deriv dtemdt
        dtemdt[k,n,2] = (dtemdz2[k,n]-((n-1)*pi/a)^2*tem[k,n,2]) # (3.3) with first term on RHS removed
        return dtemdt
    end

    function nonlin_omg_eq(k,n,a,Ra,Pr,tem,omg,domgdz2,domgdt)
        # Update time deriv domgt
        domgdt[k,n,2] = Ra*Pr*((n-1)*pi/a)*tem[k,n,2] + Pr*(domgdz2[k,n] - ((n-1)*pi/a)^2*omg[k,n,2]) # same as (3.4)
        return domgdt 
    end    

    function adamsbashforth(n, y, dydt, dt)
        y[:,n,2] = y[:,n,1] + 0.5*dt*(3*dydt[:,n,2]-dydt[:,n,1]) # (2.18)
        return y
    end

    function diagnostics(m,nz,nout,time,tem,omg,psi)
        # track n=1 mode over time at z = nz/3
        @printf("time: %.2f    tem: %.5e    omg: %.5e    psi: %.5e \n ", time, log(abs(tem[round.(Int,nz/3),2,2]))-log(abs(tem[round.(Int,nz/3),2,1])),log(abs(omg[round.(Int,nz/3),2,2]))-log(abs(omg[round.(Int,nz/3),2,1])),log(abs(psi[round.(Int,nz/3),2,2]))-log(abs(psi[round.(Int,nz/3),2,1])))
    end

    function prepare(dtemdt, domgdt, tem, omg, psi)
        dtemdt[:,:,1] = dtemdt[:,:,2]
        domgdt[:,:,1] = domgdt[:,:,2]
        tem[:,:,1] = tem[:,:,2]
        omg[:,:,1] = omg[:,:,2]
        psi[:,:,1] = psi[:,:,2]
        dtemdt[:,:,2] .=0
        domgdt[:,:,2] .=0
        return dtemdt, domgdt, tem, omg, psi
    end

    function linear_solver(z, dz, nz, nn, nt, nout, dt, a, Ra, Pr, psi, tem, omg)
        m = 0
        time = 0
        dtemdz2 = zeros(nz,nn+1)
        domgdz2 = zeros(nz,nn+1)
        dtemdt = zeros(size(tem))
        domgdt = zeros(size(omg))
        while m<=nt
            for k=2:1:nz-1
                for n=2:1:nn+1
                    dtemdz2 = second_deriv(k ,n, dz, tem, dtemdz2)
                    domgdz2 = second_deriv(k ,n, dz, omg, domgdz2)
                    dtemdt = lin_tem_eq(k,n,a,tem,psi,dtemdz2,dtemdt)
                    domgdt = lin_omg_eq(k,n,a,Ra,Pr,tem,omg,domgdz2,domgdt)
                end
            end

            # display(dtemdt[:,:,2])
            # Update tem and omg using Adams Bashforth time integration
            for n=2:1:nn+1
                tem = adamsbashforth(n, tem, dtemdt, dt)
                omg = adamsbashforth(n, omg, domgdt, dt)
            end

            # Update psi using poisson solver
            for n=2:1:nn+1
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

    function nonlinear_solver(z, dz, nz, nn, nt, nout, dt, a, Ra, Pr, psi, tem, omg, initOn, saveDir)
        if initOn==1
            rm(saveDir,recursive=true, force=true)
            mkdir(saveDir)
            data_utils.save_inputs(saveDir,nz,nn,a,Ra,Pr,dt,nt,nout)
            ndata = 0
            time = 0
            dtemdt = zeros(size(tem))
            domgdt = zeros(size(omg))
            tem = routines.initial_nonlinear_tem(nz,nn,z,tem)
        elseif initOn==0
            nz,nn,a,Ra,Pr,dt,nt,nout = data_utils.load_inputs(saveDir)
            time, ndata = data_utils.load_outputs(saveDir)
            dtemdt, domgdt, tem, omg, psi = data_utils.load_data(saveDir,ndata-1)
        end
        m = 0
        dtemdz1 = zeros(nz,nn+1)
        domgdz1 = zeros(nz,nn+1) 
        dpsidz1 = zeros(nz,nn+1)
        dtemdz2 = zeros(nz,nn+1)
        domgdz2 = zeros(nz,nn+1)
        while m<nt
            for k=2:1:nz-1
                # Linear terms
                for n=1:1:nn+1
                    dtemdz1 = first_deriv(k ,n, dz, tem, dtemdz1)
                    domgdz1 = first_deriv(k ,n, dz, omg, domgdz1)
                    dpsidz1 = first_deriv(k ,n, dz, psi, dpsidz1)
                    dtemdz2 = second_deriv(k ,n, dz, tem, dtemdz2)
                    domgdz2 = second_deriv(k ,n, dz, omg, domgdz2)
                    dtemdt = nonlin_tem_eq(k,n,a,tem,psi,dtemdz2,dtemdt)
                    domgdt = nonlin_omg_eq(k,n,a,Ra,Pr,tem,omg,domgdz2,domgdt)
                end
                # if k==33
                #     display(dtemdt[33,:,:])
                # end
                # Nonlinear terms
                for n1=2:1:nn+1
                    # Zeroth mode
                    dtemdt[k,1,2] += -pi/(2*a)*(n1-1)*(dpsidz1[k,n1]*tem[k,n1,2]+psi[k,n1,2]*dtemdz1[k,n1])
                end
                for n=2:1:nn+1
                    # n'= 0 mode
                    dtemdt[k,n,2] += -(n-1)*pi/a*psi[k,n,2]*dtemdz1[k,1]
                    # n'> 0
                    for n1=2:1:nn+1
                        n2 = zeros(Int8,3)
                        tem_term = zeros(Float64,3)
                        omg_term = zeros(Float64,3)
                        n2[1] = (n-1)-(n1-1)
                        n2[2] = (n-1)+(n1-1)
                        n2[3] = (n1-1)-(n-1)
                        for i=1:1:length(n2)
                            # Check if 2<=n<=nn+1, no contribution if not
                            if i==1 && n2[i]>=2 && n2[i]<=nn+1
                                tem_term[i] = -(n1-1)*dpsidz1[k,n2[i]]*tem[k,n1,2]+n2[i]*psi[k,n2[i],2]*dtemdz1[k,n1]
                                omg_term[i] = -(n1-1)*dpsidz1[k,n2[i]]*omg[k,n1,2]+n2[i]*psi[k,n2[i],2]*domgdz1[k,n1]
                            elseif i==2 && n2[i]>=2 && n2[i]<=nn+1
                                tem_term[i] = (n1-1)*dpsidz1[k,n2[i]]*tem[k,n1,2]+n2[i]*psi[k,n2[i],2]*dtemdz1[k,n1]
                                omg_term[i] = -(n1-1)*dpsidz1[k,n2[i]]*omg[k,n1,2]-n2[i]*psi[k,n2[i],2]*domgdz1[k,n1]
                            elseif i==3 && n2[i]>=2 && n2[i]<=nn+1
                                tem_term[i] = (n1-1)*dpsidz1[k,n2[i]]*tem[k,n1,2]+n2[i]*psi[k,n2[i],2]*dtemdz1[k,n1]
                                omg_term[i] = (n1-1)*dpsidz1[k,n2[i]]*omg[k,n1,2]-n2[i]*psi[k,n2[i],2]*domgdz1[k,n1]
                            end
                        end
                        dtemdt[k,n,2] += -pi/(2*a)*(sum(tem_term))
                        domgdt[k,n,2] += -pi/(2*a)*(sum(omg_term))
                    end
                end
            end
            
            # Update tem and omg using Adams Bashforth time integration
            for n=2:1:nn+1
                tem = adamsbashforth(n, tem, dtemdt, dt)
                omg = adamsbashforth(n, omg, domgdt, dt)
            end
            # Update psi using poisson solver
            for n=2:1:nn+1
                psi[:,n,2] = poisson(omg[:,n,2],nz,dz,n,a) # (3.5)
            end
            # Save and print diagnostics
            if mod(m,nout)==0 
                diagnostics(m, nz, nout, time, tem, omg, psi)
                # Prepare values for next timestep
                dtemdt, domgdt, tem, omg, psi = prepare(dtemdt, domgdt, tem, omg, psi)
                data_utils.save_data(saveDir,ndata,dtemdt, domgdt, tem, omg, psi)
                ndata+=1
            else
                dtemdt, domgdt, tem, omg, psi = prepare(dtemdt, domgdt, tem, omg, psi)
            end
            m+=1
            time += dt
            data_utils.save_outputs(saveDir, time, ndata)
        end
        return dtemdt, domgdt, tem, omg, psi
    end

    function cosines(a, x, nn, nx)
        # Compute cosines and sines 
        cosa = zeros(nn+1,nx)
        for n=1:1:nn+1
            cosa[n,:] = cos.((n-1)*pi*x/a)
        end
        return cosa
    end

    function sines(a, x, nn, nx)
        # Compute cosines and sines 
        sina = zeros(nn+1,nx)
        for n=1:1:nn+1
            sina[n,:] = sin.((n-1)*pi*x/a)
        end
        return sina
    end    

    function ict(nn,nx,nz,cosa,coeffs,outfun,zeroth=1)
        for i = 1:1:nx
            for k = 1:1:nz
                if zeroth==1 # include zeroth mode
                    for n=1:1:nn+1
                        outfun[k,i] += coeffs[k,n,2]*cosa[n,i]
                    end
                elseif zeroth==0 # exclude zeroth mode
                    for n=2:1:nn+1
                        outfun[k,i] += coeffs[k,n,2]*cosa[n,i]
                    end
                end
            end
        end
        return outfun
    end

    function ist(nn,nx,nz,sina,coeffs,outfun)
        for i = 1:1:nx
            for k = 1:1:nz
                for n=1:1:nn+1
                    outfun[k,i] += coeffs[k,n,2]*sina[n,i]
                end
            end
        end
        return outfun
    end

end



