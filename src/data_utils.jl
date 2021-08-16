module data_utils
    using Printf
    using HDF5
    export all

    function save_inputs(saveDir,nz,nn,a,Ra,Pr,dt,nt,nout)
        isdir(saveDir) || mkdir(saveDir)
        saveName=string(saveDir,"/inputs.h5")
        fid = HDF5.h5open(saveName,"w")
        fid["nz"] = nz
        fid["nn"] = nn
        fid["a"] = a
        fid["Ra"] = Ra
        fid["Pr"] = Pr
        fid["dt"] = dt
        fid["nt"] = nt
        fid["nout"] = nout
    end

    function save_data(saveDir, ndata, dtemdt, domgdt, tem, omg, psi)
        isdir(saveDir) || mkdir(saveDir)
        saveName=string(saveDir,"/data_",lpad(ndata,4,"0"),".h5")
        fid = HDF5.h5open(saveName,"w")
        fid["dtemdt"] = dtemdt
        fid["domgdt"] = domgdt
        fid["psi"] = psi
        fid["omg"] = omg
        fid["tem"] = tem
    end

    function save_outputs(saveDir,time,ndata)
        isdir(saveDir) || mkdir(saveDir)
        saveName=string(saveDir,"/outputs.h5")
        fid = HDF5.h5open(saveName,"w")
        fid["time"] = time
        fid["ndata"] = ndata
    end

    function load_inputs(saveDir)
        saveName=string(saveDir,"/inputs.h5")
        fid = HDF5.h5open(saveName,"r")
        nz = HDF5.read(fid["nz"])
        nn = HDF5.read(fid["nn"])
        a = HDF5.read(fid["a"])
        Ra = HDF5.read(fid["Ra"])
        Pr = HDF5.read(fid["Pr"])
        dt = HDF5.read(fid["dt"])
        nt = HDF5.read(fid["nt"])
        nout = HDF5.read(fid["nout"])
        return nz,nn,a,Ra,Pr,dt,nt,nout
    end

    function load_data(saveDir, ndata)
        saveName=string(saveDir,"/data_",lpad(ndata,4,"0"),".h5")
        fid = HDF5.h5open(saveName,"r")
        dtemdt = HDF5.read(fid["dtemdt"])
        domgdt = HDF5.read(fid["domgdt"])
        psi = HDF5.read(fid["psi"])
        omg = HDF5.read(fid["omg"])
        tem = HDF5.read(fid["tem"])
        return dtemdt, domgdt, tem, omg, psi
    end

    function load_outputs(saveDir)
        saveName=string(saveDir,"/outputs.h5")
        fid = HDF5.h5open(saveName,"r")
        time = HDF5.read(fid["time"])
        ndata = HDF5.read(fid["ndata"])
        return time, ndata
    end

end