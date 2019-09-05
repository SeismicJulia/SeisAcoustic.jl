function one_step_forward!(spt2::Snapshot{Tv}, spt1::Snapshot{Tv}, params::TdParams{Ti,Tv},
         fc::Vector{Tv}, tmp::Vector{Tv}, tmp_z1::Vector{Tv}, tmp_z2::Vector{Tv},
         tmp_x1::Vector{Tv}, tmp_x2::Vector{Tv}) where {Ti<:Int64, Tv<:AbstractFloat}

    # # local variables
    # fc = params.fc

    # number of byte for each element
    esize = sizeof(params.data_format)

    # save the pressure field
    tmp .= spt1.pz .+ spt1.px
    p    = pointer(tmp)

    # update vz=================================================================
    ilower = 1
    iupper = params.Nz
    for ix = 1 : params.Nx

        # copy one column
        copyto!(tmp_z1, 1, tmp, ilower, params.Nz)
        fill!(tmp_z2, 0.0)

        # dpdz
        for iz = 1 : params.order-1
            for idx = 1 : iz
                tmp_z2[iz] -= fc[idx] * tmp_z1[iz-idx+1]
            end
            for idx = 1 : params.order
                tmp_z2[iz] += fc[idx] * tmp_z1[iz+idx]
            end
        end

        for iz = params.order : params.Nz-params.order
            for idx = 1 : params.order
                tmp_z2[iz] += fc[idx] * (tmp_z1[iz+idx] - tmp_z1[iz-idx+1])
            end
        end

        for iz = params.Nz-params.order+1 : params.Nz
            for idx = 1 : params.order
                tmp_z2[iz] -= fc[idx] * tmp_z1[iz-idx+1]
            end
            for idx = 1 : params.Nz - iz
                tmp_z2[iz] += fc[idx] * tmp_z1[iz+idx]
            end
        end

        # update one column
        idx    = 0
        for iz = ilower : iupper
            idx = idx + 1
            spt2.vz[iz] = params.MvzBvz[idx]*spt1.vz[iz] + params.MvzBp[iz]*tmp_z2[idx]
        end

        # prepare for next column
        ilower = ilower + params.Nz
        iupper = iupper + params.Nz
    end

    # update vx=================================================================
    shift = 0
    for iz = 1 : params.Nz

        # copy one row
        BLAS.blascopy!(params.Nx, p+shift, params.Nz, tmp_x1, 1)
        fill!(tmp_x2, 0.0)

        # dpdx
        for ix = 1 : params.order-1
            for idx = 1 : ix
                tmp_x2[ix] -= fc[idx] * tmp_x1[ix-idx+1]
            end
            for idx = 1 : params.order
                tmp_x2[ix] += fc[idx] * tmp_x1[ix+idx]
            end
        end

        for ix = params.order : params.Nx-params.order
            for idx = 1 : params.order
                tmp_x2[ix] += fc[idx] * (tmp_x1[ix+idx] - tmp_x1[ix-idx+1])
            end
        end

        for ix = params.Nx-params.order+1 : params.Nx
            for idx = 1 : params.order
                tmp_x2[ix] -= fc[idx] * tmp_x1[ix-idx+1]
            end
            for idx = 1 : params.Nx - ix
                tmp_x2[ix] += fc[idx] * tmp_x1[ix+idx]
            end
        end

        idx = iz
        for ix = 1 : params.Nx
            spt2.vx[idx] = params.MvxBvx[ix]*spt1.vx[idx] + params.MvxBp[idx]*tmp_x2[ix]
            idx= idx + params.Nz
        end

        # prepare for next row
        shift = shift + esize
    end

    # update pz ================================================================
    ilower = 1
    iupper = params.Nz

    # loop over columns
    for ix = 1 : params.Nx

        # copy one column
        copyto!(tmp_z1, 1, spt2.vz, ilower, params.Nz)
        fill!(tmp_z2, 0.0)

        # dvdz
        for iz = 1 : params.order
            for idx = 1 : iz-1
                tmp_z2[iz] -= fc[idx] * tmp_z1[iz-idx]
            end
            for idx = 1 : params.order
                tmp_z2[iz] += fc[idx] * tmp_z1[iz+idx-1]
            end
        end

        for iz = params.order+1 : params.Nz-params.order+1
            for idx = 1 : params.order
                tmp_z2[iz] += fc[idx] * (tmp_z1[iz+idx-1]-tmp_z1[iz-idx])
            end
        end

        for iz = param.Nz-params.order+2 : params.Nz
            for idx = 1 : params.order
                tmp_z2[iz] -= fc[idx] * tmp_z1[iz-idx]
            end
            for idx = 1 : params.Nz-iz+1
                tmp_z2[iz] += fc[idx] * tmp_z1[iz+idx-1]
            end
        end

        idx    = 0
        for iz = ilower : iupper
            idx = idx + 1
            spt2.pz[iz] = params.MpzBpz[idx]*spt1.pz[iz] + params.MpzBvz[iz]*tmp_z2[idx]
        end

        # prepare for next column
        ilower = ilower + params.Nz
        iupper = iupper + params.Nz
    end

    # update px=================================================================
    shift = 0
    p = pointer(spt2.vx)
    for iz = 1 : params.Nz

        # copy one row
        BLAS.blascopy!(params.Nx, p+shift, params.Nz, tmp_x1, 1)
        fill!(tmp_x2, 0.0)

        # dvdx
        for ix = 1 : params.order
            for idx = 1 : ix
                tmp_x2[ix] -= fc[idx] * tmp_x1[ix-idx]
            end
            for idx = 1 : params.order
                tmp_x2[ix] += fc[idx] * tmp_x1[ix+idx-1]
            end
        end

        for ix = params.order+1 : params.Nx-params.order+1
            for idx = 1 : params.order
                tmp_x2[ix] += fc[idx] * (tmp_x1[ix+idx-1] - tmp_x1[ix-idx])
            end
        end

        for ix = params.Nx-params.order+2 : params.Nx
            for idx = 1 : params.order
                tmp_x2[ix] -= fc[idx] * tmp_x1[ix-idx]
            end
            for idx = 1 : params.Nx - ix + 1
                tmp_x2[ix] += fc[idx] * tmp_x1[ix+idx-1]
            end
        end

        idx = iz
        for ix = 1 : params.Nx
            spt2.vx[idx] = params.MvxBvx[ix]*spt1.vx[idx] + params.MvxBp[idx]*tmp_x2[ix]
            idx= idx + params.Nz
        end

        # prepare for next row
        shift = shift + esize
    end

end
