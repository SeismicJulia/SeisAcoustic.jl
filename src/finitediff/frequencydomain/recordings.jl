
struct MonochromaticRecordings{Ti<:Int64, Tv<:Float64}
     nr      :: Ti                  # number of receiver
     omega   :: Tv
     irz     :: Vector{Ti}          # vertical grid index of receivers
     irx     :: Vector{Ti}          # horizontal grid index of receivers
     spt2rec :: Vector{Ti}          # index mapping of receiver to snapshot
     p       :: Vector{Complex{Tv}} # recordings of pressure field
end

function MonochromaticRecordings(rz::Vector, rx::Vector, omega,
         params::FdParams; location_flag="index")

    # number of receivers
    nr  = length(rz)
    if length(rx) != nr
       error("length of receiver coordinate vector doesn't match")
    end

    irz = zeros(Int64, nr)
    irx = zeros(Int64, nr)

    # location are provided as index
    if location_flag == "index"
       for i = 1 : nr
           irz[i] = round(Int64, rz[i])
           irx[i] = round(Int64, rx[i])
       end

    # location are given as distance
    elseif location_flag == "distance"
       for i = 1 : nr
           irz[i] = round(Int64, rz[i]/params.h) + 1
           irx[i] = round(Int64, rx[i]/params.h) + 1
       end
    else
       error("wrong specification of receiver location")
    end

    # error checking
    for i = 1 : nr
        if irz[i] > params.nz || irx[i] > params.nx || irz[i] < 1 || irx[i] < 1
           error("receiver located outside of modeling area")
        end
    end

    # can't put receivers on the free surface
    if params.free_surface
       for i = 1 : nr
           if irz[i] == 1
              error("can't put receiver at free surface, move it deeper")
           end
       end
    end

    # the auxillary vector mapping snapshot to recordings
    spt2rec = zeros(Int64, nr)
    for i = 1 : nr
        spt2rec[i] = (irx[i]+params.npml-1) * params.Nz + irz[i] + params.ntop
    end

    # radian frequency
    omega = convert(Float64, omega)

    # initalize empty recordings
    rec= MonochromaticRecordings(nr, omega, irz, irx, spt2rec,
         zeros(Complex{Float64}, nr))

    return rec
end

"""
   sampling one snashot to fill the recordings
"""
function sample_spt2rec!(rec::MonochromaticRecordings, p::Vector{Complex{Tv}}) where {Tv<:Float64}

    # loop over receivers
    for i = 1 : rec.nr
        idx = rec.spt2rec[i]
        rec.p[i] = p[idx]
    end

    return nothing
end

"""
   the adjoint of sampling operator, inject the complex conjugate of recordings to snapshot
"""
function inject_rec2spt!(u::Vector{Complex{Tv}}, rec::MonochromaticRecordings) where {Tv<:Float64}

    # zerolize the vector u
    fill!(u, zero(Complex{Float64}))

    # loop over receivers
    for i = 1 : rec.nr
        idx    = rec.spt2rec[i]
        u[idx] = conj(rec.p[i])
    end

    return nothing
end

"""
   compute the residue between synthetic monochromatic slice and observed one
"""
function get_residue(dsyn::MonochromaticRecordings, dobs::MonochromaticRecordings)

    # check radian frequency
    if dsyn.omega != dobs.omega
       error("radian frequency doesn't match")
    end
    
    # check the number of receivers
    if dsyn.nr != dobs.nr
       error("number of receivers doesn't match")
    end

    # check the receiver's location
    for i = 1 : dsyn.nr
        if dsyn.irz[i] != dobs.irz[i] || dsyn.irx[i] != dobs.irx[i]
           error("the location of receivers doesn't match")
        end
    end

    p = zeros(Complex{Float64}, dsyn.nr)
    for i = 1 : dsyn.nr
        p[i] = dsyn.p[i] - dobs.p[i]
    end

    return MonochromaticRecordings(dsyn.nr, dsyn.omega,
           copy(dsyn.irz), copy(dsyn.irx), copy(dsyn.spt2rec), p)

end
