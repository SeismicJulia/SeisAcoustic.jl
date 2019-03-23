# convert IBM data to IEEE
primitive type IBMFloat32 32 end

ieeeOfPieces(fr::UInt32, exp::Int32, sgn::UInt32) = reinterpret(Float32, convert(UInt32,fr >>> 9) | convert(UInt32,exp << 23) | sgn) :: Float32

"""
   convert IBM Float32 to IEEE Float32
"""
function convert(::Type{Float32}, ibm::IBMFloat32)

    local fr::UInt32 = ntoh(reinterpret(UInt32, ibm))
    local sgn::UInt32 = fr & 0x80000000 # save sign
    fr <<= 1 # shift sign out
    local exp::Int32 = convert(Int32,fr >>> 25) # save exponent
    fr <<= 7 # shift exponent out

    if (fr == convert(UInt32,0))
       zero(Float32)
    else
       # normalize the signficand
       local norm::UInt32 = leading_zeros(fr)
       fr <<= norm
       exp = (exp << 2) - 130 - norm

       # exp <= 0 --> ieee(0,0,sgn)
       # exp >= 255 --> ieee(0,255,sgn)
       # else -> ieee(fr<<1, exp, sgn)
       local clexp::Int32 = exp & convert(Int32,0xFF)
       ieeeOfPieces(clexp == exp ? fr << 1 : convert(UInt32,0), clexp, sgn)
    end
end

"""
   the mapping between edbcb and ascii code, as the text header of segy file is
encode with edbcb format, so we need first translate it into ascii format,
then convert to character.
"""
const edbcb2ascii_table = [0, 1, 2, 3, 156, 9, 134,127,151,141,142, 11, 12, 13, 14, 15,
       16, 17, 18, 19,157,133,  8,135, 24, 25,146,143, 28, 29, 30, 31,
       128,129,130,131,132, 10, 23, 27,136,137,138,139,140,  5,  6,  7,
       144,145, 22,147,148,149,150,  4,152,153,154,155, 20, 21,158, 26,
       32,160,161,162,163,164,165,166,167,168, 91, 46, 60, 40, 43, 33,
       38,169,170,171,172,173,174,175,176,177, 93, 36, 42, 41, 59, 94,
       45, 47,178,179,180,181,182,183,184,185,124, 44, 37, 95, 62, 63,
       186,187,188,189,190,191,192,193,194, 96, 58, 35, 64, 39, 61, 34,
       195, 97, 98, 99,100,101,102,103,104,105,196,197,198,199,200,201,
       202,106,107,108,109,110,111,112,113,114,203,204,205,206,207,208,
       209,126,115,116,117,118,119,120,121,122,210,211,212,213,214,215,
       216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,
       123, 65, 66, 67, 68, 69, 70, 71, 72, 73,232,233,234,235,236,237,
       125, 74, 75, 76, 77, 78, 79, 80, 81, 82,238,239,240,241,242,243,
       92,159, 83, 84, 85, 86, 87, 88, 89, 90,244,245,246,247,248,249,
       48, 49, 50, 51, 52, 53, 54, 55, 56, 57,250,251,252,253,254,255]

"""
   create the mapping table which maps ascii code to edbcb code.
"""
function create_ascii2edbcb_table(e2a::Vector{Int64})

    a2e = zeros(Int64, 256)
    for i = 1 : length(e2a)
       	a2e[e2a[i]+1] = i-1
    end
    return a2e
end
const ascii2edbcb_table = create_ascii2edbcb_table(edbcb2ascii_table);

"""
   conversions between ascii and edbcb table
"""
function conversion_between_edbcb_ascii(a::Vector{UInt8}; flag="edbcb2ascii")

    r = zeros(UInt, length(a))
    if flag == "edbcb2ascii"

       idx = 1
       for c in a
           c1 = edbcb2ascii_table[c+1]
           r[idx] = c1
           idx = idx + 1
       end

    elseif flag == "ascii2edbcb"
       idx = 1
       for c in a
           c1 = ascii2edbcb_table[c+1]
           r[idx] = c1
           idx = idx + 1
       end

    end

    return convert(Vector{UInt8}, r)
end

"""
   read text header
"""
function print_text_header(path::String)

    # read the text header as UInt8
	  fid = open(path, "r")
		thdr = read(fid, 3200)
		close(fid)

    # convert edbcb coding to ascii coding
    thdr = conversion_between_edbcb_ascii(thdr; flag="edbcb2ascii")

    # convert each element to character
		thdr = convert(Vector{Char}, thdr)

    #print row-by-row and each row consist of 80 character
		for i = 1 : 40
			  il = (i-1)*80+1
				iu = il+80-1
			  println(join(thdr[il:iu]))
		end

		return thdr
end

"""
   put the text header to segy file
"""
function put_text_header(path::String, thdr::Vector{Char})

    #convert string of Char to a vector of UInt
    a = convert(Vector{UInt8}, thdr)

    # convert ascii coding to edbcb coding
    a = conversion_between_edbcb_ascii(a; flag="ascii2edbcb")

    # write to file
    fid = open(path, "w")
    write(fid, a)
    close(fid)

		return nothing
end

"""
   composite type of file header
*   ntrpe : number of traces per ensemble
*   dt    : sample interval in microseconds
*   dtfr  : sample interval in microseconds in original field recording
*   ns    : number of samples per data trace
*   nsfr  : number of samples per data tracf in original field recording
*   fmtc  : data sample format code(1=>IBM; 5=>IEEE)
*   expf  : expected number of trace per trace ensemble (such as CMP fold)
*   trsc  : trace sorting code(1=>as recorded, 2=>CDP, 3=>single fold profile, 4=>post stack
           5=>CSP, 6=>CRP, 7=>common offset, 8=>CMP, 9=>common conversion point)
*   sfs   : sweep frequency at start (Hz)
*   sfe   : sweep frequency at end (Hz)
*   slen  : sweep length (ms)
*   styp  : sweep type code (1=>linear, 2=>parabolic, 3=>exponential)
*   tnumsc: trace number sweep channel
*   stalens: sweep trace taper length in milliseconds at start
*   stalene: sweep trace taper length in milliseconds at end
*   taty   : taper type(1=>linear, 2=>cos^2, 3=>other)
*   corr   : correlated trace(1=>no, 2=>yes)
*   rgc    : gain recovered(1=>yes, 2=>no)
*   arm    : amplitude recover method(1=>none, 2=>spherical divergence, 3=>AGC, 4=>other)
*   unit   : measurement system(1=>meters, 2=>Feet)
*   pol    : impulse polarity(1=>increase pressure and upward movement is negative number,
                              2=>increase pressure and upward movement is positive number)
*   vpol   : Vibratory polarity code
*   fvn    : file reversion number
*   fltf   : fixed trace length flag
*   netfh  : number of extended text file header(0=> no extende, -1=>variable, positive=>exactly that many)
"""
mutable struct FileHeader
# corresponding to 3200-3259
	  jobid  :: Int32
	  linnum :: Int32
	  renum  :: Int32
	  ntrpe  :: Int16
	  natrpe :: Int16
	  dt     :: Int16
	  dtfr   :: Int16
	  ns     :: Int16
	  nsfr   :: Int16
	  fmtc   :: Int16
	  expf   :: Int16
	  trsc   :: Int16
	  vsumc  :: Int16
	  sfs    :: Int16
	  sfe    :: Int16
	  slen   :: Int16
	  styp   :: Int16
	  tnumsc :: Int16
	  stalens:: Int16
	  stalene:: Int16
	  taty   :: Int16
	  corr   :: Int16
	  rgc    :: Int16
	  arm    :: Int16
	  unit   :: Int16
	  pol    :: Int16
	  vpol   :: Int16
# 3260-3501 is unsigned, need to skip 242 byte
	  fltf   :: Int16 # 3502-3503
	  netfh  :: Int16 # 3504-3505
# then skip 94 byte to the start point of trace header
end

# record the start position of file header item
const file_header_location = Dict(
      :jobid   => 3200,
      :linnum  => 3204,
      :renum   => 3208,
      :ntrpe   => 3212,
      :natrpe  => 3214,
      :dt      => 3216,
      :dtfr    => 3218,
      :ns      => 3220,
      :nsfr    => 3222,
      :fmtc    => 3224,
      :expf    => 3226,
      :trsc    => 3228,
      :vsumc   => 3230,
      :sfs     => 3232,
      :sfe     => 3234,
      :slen    => 3236,
      :styp    => 3238,
      :tnumsc  => 3240,
      :stalens => 3242,
      :stalene => 3244,
      :taty    => 3246,
      :corr    => 3248,
      :rgc     => 3250,
      :arm     => 3252,
      :unit    => 3254,
      :pol     => 3256,
      :vpol    => 3258,
# 3260-3501 is unsigned, need to skip 242 byte
      :fltf    => 3502,  # 3502-3503
      :netfh   => 3504); # 3504-3505
  # then skip 94 byte to the start point of trace header

function init_FileHeader()
	  fhdr = FileHeader(0,0,0,0,0,0,0,0,0,0,
	                    0,0,0,0,0,0,0,0,0,0,
	                    0,0,0,0,0,0,0,0,0,)
	  return fhdr
end

function show(io::IO, fhdr::FileHeader)
    @printf("number of traces per ensemble                               : %6d\n", fhdr.ntrpe)
    @printf("sample interval in microseconds                             : %6d\n", fhdr.dt   )
    @printf("sample interval in microseconds in original field recording : %6d\n", fhdr.dtfr )
    @printf("number of samples per data trace                            : %6d\n", fhdr.ns   )
    @printf("number of samples per data tracf in original field recording: %6d\n", fhdr.nsfr )

    if fhdr.fmtc == 1
       @printf("data sample format                                          : IBM\n"            )
    elseif fhdr.fmtc == 5
       @printf("data sample format                                          : IEEE\n"            )
    end

    @printf("expected number of trace per trace ensemble                 : %6d\n", fhdr.expf )

    if fhdr.trsc == 1
       @printf("trace sorting code                                          : as recorded\n"      )
    elseif fhdr.trsc == 2
       @printf("trace sorting code                                          : CDP\n"      )
    elseif fhdr.trsc == 3
       @printf("trace sorting code                                          : single fold profile\n")
    elseif fhdr.trsc == 4
       @printf("trace sorting code                                          : post stack\n"      )
    elseif fhdr.trsc == 5
       @printf("trace sorting code                                          : common shot point\n")
    elseif fhdr.trsc == 6
       @printf("trace sorting code                                          : common receiver point\n")
    elseif fhdr.trsc == 7
       @printf("trace sorting code                                          : common offset profile\n")
    elseif fhdr.trsc == 8
       @printf("trace sorting code                                          : CMP\n"                  )
    elseif fhdr.trsc == 9
       @printf("trace sorting code                                          : common conversion point\n")
    end

    @printf("sweep frequency at start (Hz)                               : %6d\n", fhdr.sfs  )
    @printf("sweep frequency at end (Hz)                                 : %6d\n", fhdr.sfe  )
    @printf("sweep length (ms)                                           : %6d\n", fhdr.slen )

    if fhdr.styp == 1
       @printf("sweep type                                                  : linear\n"           )
    elseif fhdr.styp == 2
       @printf("sweep type                                                  : parabolic\n"        )
    elseif fhdr.styp == 3
       @printf("sweep type                                                  : exponential\n"      )
    end

    if fhdr.corr == 1
       @printf("correlated trace                                            : NO\n"      )
    elseif fhdr.corr == 2
       @printf("correlated trace                                            : YES\n"      )
    end

    if fhdr.rgc == 1
       @printf("gain recovered                                              : YES\n"      )
    elseif fhdr.rgc == 2
       @printf("gain recovered                                              : NO\n"      )
    end

    if fhdr.arm == 1
       @printf("amplitue recover method                                     : None\n"      )
    elseif fhdr.arm == 2
       @printf("amplitue recover method                                     : spherical divergence\n")
    elseif fhdr.arm == 3
       @printf("amplitue recover method                                     : AGC\n"              )
    elseif fhdr.arm == 4
       @printf("amplitue recover method                                     : other\n"            )
    end

    if fhdr.unit == 1
       @printf("measurement system                                          : Meters\n"          )
    elseif fhdr.unit == 2
       @printf("measurement system                                          : Feet\n"          )
    end

    if fhdr.pol == 1
       @printf("polarity : increase pressure and upward movement use negative number\n"        )
    elseif fhdr.pol == 2
       @printf("polarity : increase pressure and upward movement use positive number\n"        )
    end

    if fhdr.netfh == 0
       @printf("number of extended text file header                         : none\n"             )
    elseif fhdr.netfh == -1
       @printf("number of extended text file header                         : variable\n"         )
    elseif fhdr.netfh > 0
       @printf("number of extended text file header                         : %6d\n", fhdr.netfh)
    end

    return nothing
end

"""
   read 400 byte file header of SEGY file.
"""
function read_file_header(path::String; swap_bytes=true)

    stream = open(path, "r")

    # seek to the start point of file header
		seek(stream, file_header_location[:jobid])

    # initialize a empty file header
		fhdr = init_FileHeader()

    # read file header items
    # There are unsigned space in segy file, for loop can't be used
    fhdr.jobid   = read(stream, typeof(fhdr.jobid))
    fhdr.linnum  = read(stream, typeof(fhdr.linnum))
    fhdr.renum   = read(stream, typeof(fhdr.renum))
    fhdr.ntrpe   = read(stream, typeof(fhdr.ntrpe))
    fhdr.natrpe  = read(stream, typeof(fhdr.natrpe))
    fhdr.dt      = read(stream, typeof(fhdr.dt))
    fhdr.dtfr    = read(stream, typeof(fhdr.dtfr))
    fhdr.ns      = read(stream, typeof(fhdr.ns))
    fhdr.nsfr    = read(stream, typeof(fhdr.nsfr))
    fhdr.fmtc    = read(stream, typeof(fhdr.fmtc))
    fhdr.expf    = read(stream, typeof(fhdr.expf))
    fhdr.trsc    = read(stream, typeof(fhdr.trsc))
    fhdr.vsumc   = read(stream, typeof(fhdr.vsumc))
    fhdr.sfs     = read(stream, typeof(fhdr.sfs))
    fhdr.sfe     = read(stream, typeof(fhdr.sfe))
    fhdr.slen    = read(stream, typeof(fhdr.slen))
    fhdr.styp    = read(stream, typeof(fhdr.styp))
    fhdr.tnumsc  = read(stream, typeof(fhdr.tnumsc))
    fhdr.stalens = read(stream, typeof(fhdr.stalens))
    fhdr.stalene = read(stream, typeof(fhdr.stalene))
    fhdr.taty    = read(stream, typeof(fhdr.taty))
    fhdr.corr    = read(stream, typeof(fhdr.corr))
    fhdr.rgc     = read(stream, typeof(fhdr.rgc))
    fhdr.arm     = read(stream, typeof(fhdr.arm))
    fhdr.unit    = read(stream, typeof(fhdr.unit))
    fhdr.pol     = read(stream, typeof(fhdr.pol))
    fhdr.vpol    = read(stream, typeof(fhdr.vpol))

    skip(stream, 242) # skip the unsigned space
    fhdr.fltf    = read(stream, typeof(fhdr.fltf))
    fhdr.netfh   = read(stream, typeof(fhdr.netfh))

    # convert endian
    if swap_bytes
    	 for field in fieldnames(FileHeader)
    			 setfield!(fhdr, field, bswap(getfield(fhdr, field)))
    	 end
    end

    # skip to the start point of trace header
    # skip(stream, 94)
		close(stream)

	  return fhdr
end

"""
   write 400 byte file header of SEGY file.
"""
function put_file_header(path::String, fhdr::FileHeader; swap_bytes=true)

    # the start point for file header
    stream = open(path, "w+")
    seek(stream, file_header_location[:jobid])

    # initialize a temporary file header
    # keep the original file header untouch
    fh = init_FileHeader()
    if swap_bytes
       for field in fieldnames(FileHeader)
           setfield!(fh, field, bswap(getfield(fhdr, field)))
       end
    else
       for field in fieldnames(FileHeader)
           setfield!(fh, field, getfield(fhdr, field))
       end
    end

    # write file header
    write(stream, fh.jobid)
    write(stream, fh.linnum)
    write(stream, fh.renum)
    write(stream, fh.ntrpe)
    write(stream, fh.natrpe)
    write(stream, fh.dt)
    write(stream, fh.dtfr)
    write(stream, fh.ns)
    write(stream, fh.nsfr)
    write(stream, fh.fmtc)
    write(stream, fh.expf)
    write(stream, fh.trsc)
    write(stream, fh.vsumc)
    write(stream, fh.sfs)
    write(stream, fh.sfe)
    write(stream, fh.slen)
    write(stream, fh.styp)
    write(stream, fh.tnumsc)
    write(stream, fh.stalens)
    write(stream, fh.stalene)
    write(stream, fh.taty)
    write(stream, fh.corr)
    write(stream, fh.rgc)
    write(stream, fh.arm)
    write(stream, fh.unit)
    write(stream, fh.pol)
    write(stream, fh.vpol)

    # assign 0 to 242 unsigned bytes
    write(stream, zeros(Int8, 242))

    # write the last two field
    write(stream, fh.fltf)
    write(stream, fh.netfh)

    # assign 0 to 94 unsigned bytes
    write(stream, zeros(Int8, 94))
    close(stream)

    return nothing
end

"""
* tracl: trace number within line
* tracr: trace number within file
* fldr : original field record number
* tracf: trace number within the original field record
* ep   : energy source point number
* cdp  : Ensemble number (i.e. CDP, CMP, CRP,)
* cdpt : trace number within the ensemble
* trid : trace identification code (-1=other, 0=Unknown, 1=seismic, 2=dead, 3=Dummy, 4=Time break, 5=Uphole, 6=sweep
         7=timing, 8=waterbreak, 9=near field gun signature, 10=far field gun signature, 11=pressure sensor, 12=vertical component,
         13=cross-line component, 14=in-line component, 15=rotated vertical component, 16=rotated Transversal, 17=rotated radial component,
         18=vibrator reation mass, 19=vibrator baseplate, 19=vibrator estimated ground force, 21=vibrator reference
         22=time-velocity pairs, 23=optional use)
* nva    : number of vertically summed traces yield this trace
* nhs    : number of horizontally summed traces yield this trace
* duse   : data use (1=production, 2=Test)
* offset : distance between source and receiver
* gelev  : receiver elevation(above datum is positive, under is negative)
* selev  : surface elevation at source
* sdepth : source depth below surface (positive number)
* gdel   : datum elevation at receiver group
* sdel   : datum elevation at source
* swdep  : water depth at source
* gwdep  : water depth at receiver
* scalel : scalar applied to all evevations and depths (1,10,100,1000, positive used as multiplier, negative used as divisor)
* scalco : scalar applied to all coordinate (same as above case)
* sx     : source x
* sy     : source y
* gx     : receiver x
* gy     : receiver y
* counit : coordinate unit (1=length, 2=seconds of arc, 3=decimal degrees, 4=DMS)
* wevel  : weathering velocity (ft/s or meter/s)
* swevel : subweathing velocity
* sut    : uphole time at source in milliseconds
* gut    : uphole time at receiver in milliseconds
* sstat  : source static correction in milliseconds
* gstat  : group  static correction in milliseconds
* tstat  : total static applied in milliseconds
* laga   : lag time A
* lagb   : lag time B
* delrt  : delay recording time in millliseconds
* muts   : mute time - start time in milliseconds
* mute   : mute time - end time in milliseconds
* ns     : number of samples in this trace
* dt     : sample interval in microseconds
* gain   : gain type of field instruments (1=fixed, 2=binary, 3=floating point)
* igc    : instrument gain constant (dB)
* igi    : instrument early or initial gain (dB)
* corr   : correlated (1=no 2=yes)
* sfs    : sweep frequency at start
* sfe    : sweep frequency at end
* slen   : sweep length
* styp   : sweep type (1=linear 2=parabolic 3=exponential 4=other)
* stas   : sweep trace taper length at start in milliseconds
* stae   : sweep trace taper length at end in milliseconds
* tatyp  : taper type (1=linear, 2=cos^2, 3=other)
* afilf  : alias filter frequency (Hz), if applied
* afils  : alias filter slope
* nofilf : Notch filter frequency
* nofils : Notch filter slope
* lcf    : low-cut frequency
* hcf    : high-cut frequency
* lcs    : low-cut slope
* hcs    : high-cut slope
* year
* day
* hour
* minute
* sec
* timbas : timme basis code
* trwf   : trace weighting factor
* grnors : geophone group number of roll switch position one
* grnofr : geophone group number of roll switch position one within original field
* grnlof : geophone group number of last trace within original field
* gaps   : gap size
* otrav  : over travel
* CDPx   : CDP x coordinate, scala for coordinate need to apply
* CDPy   : CDP y
* inlineNum : for 3D post-stack data, inline number
* crosslineNum : for 3D post-stack data, crossline number
* shotNum: shot point number for poststack case
* scalarOnShotNum : applied to shot point number (=0, not scaled)
* didf   : device identifier
* scalart: scalar applied to acquisition time
* stype  : source type(1=vv, 2=vc, 3=vi, 4=iv, 5=ic, 6=ii, 7=div, 8=dic, 9=dii{v=vibrator v=vertical c=crossline i=inline i=impulse di=distributed impulse})
          source energy direction and source measurement, source measurement unit need to be handled.
"""
mutable struct TraceHeader
		tracl   :: Int32
		tracr   :: Int32
		fldr    :: Int32
		tracf   :: Int32
		ep      :: Int32
		cdp     :: Int32
		cdpt    :: Int32
    trid    :: Int16
		nva     :: Int16
		nhs     :: Int16
		duse    :: Int16
		offset  :: Int32
		gelev   :: Int32
		selev   :: Int32
		sdepth  :: Int32
		gdel    :: Int32
		sdel    :: Int32
		swdep   :: Int32
		gwdep   :: Int32
		scalel  :: Int16
		scalco  :: Int16
		sx      :: Int32
		sy      :: Int32
		gx      :: Int32
		gy      :: Int32
		counit  :: Int16
		wevel   :: Int16
		swevel  :: Int16
		sut     :: Int16
		gut     :: Int16
		sstat   :: Int16
		gstat   :: Int16
		tstat   :: Int16
		laga    :: Int16
		lagb    :: Int16
		delrt   :: Int16
		muts    :: Int16
		mute    :: Int16
		ns      :: Int16
		dt      :: Int16
		gain    :: Int16
		igc     :: Int16
		igi     :: Int16
		corr    :: Int16
		sfs     :: Int16
		sfe     :: Int16
		slen    :: Int16
		styp    :: Int16
		stas    :: Int16
		stae    :: Int16
		tatyp   :: Int16
		afilf   :: Int16
		afils   :: Int16
		nofilf  :: Int16
		nofils  :: Int16
		lcf     :: Int16
		hcf     :: Int16
		lcs     :: Int16
		hcs     :: Int16
		year    :: Int16
		day     :: Int16
		hour    :: Int16
		minute  :: Int16
		sec     :: Int16
		timbas  :: Int16
		trwf    :: Int16
		grnors  :: Int16
		grnofr  :: Int16
		grnlof  :: Int16
		gaps    :: Int16
		otrav   :: Int16
		CDPx    :: Int32
		CDPy    :: Int32
		inlineNum       :: Int32
		crosslineNum    :: Int32
		shotNum         :: Int32
		scalarOnShotNum :: Int16
		tvmu            :: Int16
    didf            :: Int16
		scalart         :: Int16
		stype           :: Int16
end

const trace_header_location = Dict(
    :tracl  => 0,
    :tracr  => 4,
    :fldr   => 8,
    :tracf  => 12,
    :ep     => 16,
    :cdp    => 20,
    :cdpt   => 24,
    :trid   => 28,
    :nva    => 30,
    :nhs    => 32,
    :duse   => 34,
    :offset => 36,
    :gelev  => 40,
    :selev  => 44,
    :sdepth => 48,
    :gdel   => 52,
    :sdel   => 56,
    :swdep  => 60,
    :gwdep  => 64,
    :scalel => 68,
    :scalco => 70,
    :sx     => 72,
    :sy     => 76,
    :gx     => 80,
    :gy     => 84,
    :counit => 88,
    :wevel  => 90,
    :swevel => 92,
    :sut    => 94,
    :gut    => 96,
    :sstat  => 98,
    :gstat  => 100,
    :tstat  => 102,
    :laga   => 104,
    :lagb   => 106,
    :delrt  => 108,
    :muts   => 110,
    :mute   => 112,
    :ns     => 114,
    :dt     => 116,
    :gain   => 118,
    :igc    => 120,
    :igi    => 122,
    :corr   => 124,
    :sfs    => 126,
    :sfe    => 128,
    :slen   => 130,
    :styp   => 132,
    :stas   => 134,
    :stae   => 136,
    :tatyp  => 138,
    :afilf  => 140,
    :afils  => 142,
    :nofilf => 144,
    :nofils => 146,
    :lcf    => 148,
    :hcf    => 150,
    :lcs    => 152,
    :hcs    => 154,
    :year   => 156,
    :day    => 158,
    :hour   => 160,
    :minute => 162,
    :sec    => 164,
    :timbas => 166,
    :trwf   => 168,
    :grnors => 170,
    :grnofr => 172,
    :grnlof => 174,
    :gaps   => 176,
    :otrav  => 178,
    :CDPx   => 180,
    :CDPy   => 184,
    :inlineNum       => 188,
    :crosslineNum    => 192,
    :shotNum         => 196,
    :scalarOnShotNum => 200,
    :tvmu            => 202,
# 10 bytes gap
    :didf            => 214,
    :scalart         => 216,
    :stype           => 218)
# 20 bytes gap

function init_TraceHeader()
	  h = TraceHeader(
	  0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,0,0,
	  0,0,0,0,0,0,0,0,0,0,
	  0)

	  return h
end

# show the content of trace header
function show(io::IO, h::TraceHeader)
    @printf("offset(distance between source and receiver)                : %6d\n", h.offset)
    @printf("receiver elevation (above datum is positive)                : %6d\n", h.gelev )
    @printf("source elevation (above datum is positive)                  : %6d\n", h.selev )
    @printf("source depth below surface(positve number)                  : %6d\n", h.sdepth)
    @printf("datum elevation at receiver group                           : %6d\n", h.gdel  )
    @printf("datum elevation at source location                          : %6d\n", h.sdel  )
    @printf("scalar applied to all evevations and depths(1,10,100...)    : %6d\n", h.scalel)
    @printf("scalar applied to all coordinate                            : %6d\n", h.scalco)
    @printf("source x                                                    : %6d\n", h.sx    )
    @printf("source y                                                    : %6d\n", h.sy    )
    @printf("receiver x                                                  : %6d\n", h.gx    )
    @printf("receiver y                                                  : %6d\n", h.gy    )

    if h.counit == 1
       @printf("coordinate unit                                             : length\n")
    elseif h.counit == 2
       @printf("coordinate unit                                             : seconds of arc\n")
    elseif h.counit == 3
       @printf("coordinate unit                                             : decimal degrees\n")
    elseif h.counit == 4
       @printf("coordinate unit                                             : DMS\n")
    end

    @printf("source static correction in milliseconds                    : %6d\n", h.sstat  )
    @printf("group  static correction in milliseconds                    : %6d\n", h.gstat  )
    @printf("total static applied in milliseconds                        : %6d\n", h.tstat  )
    @printf("number of samples in this trace                             : %6d\n", h.ns     )
    @printf("sample interval in microseconds                             : %6d\n", h.dt     )
    @printf("CDP x coordinate, scala for coordinate need to apply        : %6d\n", h.CDPx   )
    @printf("CDP y coordinate, scala for coordinate need to apply        : %6d\n", h.CDPy   )
    @printf("for 3D post-stack data, inline number                       : %6d\n", h.inlineNum   )
    @printf("for 3D post-stack data, crossline number                    : %6d\n", h.crosslineNum)

    return nothing
end

"""
   read the header of whole traces, return vector of TraceHeader
"""
function read_all_traces_header(path::String; swap_bytes=true)

    println("===================================================================\n")
		println("=======================Print text header===========================\n")
		text_header = print_text_header(path)

		file_header = read_file_header(path; swap_bytes=swap_bytes)
		println("===================================================================\n")
		println("=======================Print file header===========================\n")
    print(file_header)

    if file_header.netfh == -1
			 error("Does not support variable text file header")
    elseif file_header.netfh == 0
       file_hsize = 3600
    elseif file_header.netfh > 0
       file_hsize = 400 + 3200 * file_header.netfh
		end

    # compute the number of traces in this file
    file_size = filesize(path)
    num_samples = Int64(file_header.ns)
    trace_data_size = num_samples * 4
    num_traces  = round(Int64, (file_size - file_hsize) / (240 + trace_data_size))

    @printf("number of traces in this file                               : %6d\n", num_traces)

    fid = open(path, "r")
    seek(fid, file_hsize)

    vector_trace_header = Vector{TraceHeader}(num_traces)
		for i = 1 : num_traces
        h = init_TraceHeader()
        h.tracl  = read(fid, Int32)
        h.tracr  = read(fid, Int32)
        h.fldr   = read(fid, Int32)
        h.tracf  = read(fid, Int32)
        h.ep     = read(fid, Int32)
        h.cdp    = read(fid, Int32)
        h.cdpt   = read(fid, Int32)
        h.trid   = read(fid, Int16)
        h.nva    = read(fid, Int16)
        h.nhs    = read(fid, Int16)
        h.duse   = read(fid, Int16)
        h.offset = read(fid, Int32)
        h.gelev  = read(fid, Int32)
        h.selev  = read(fid, Int32)
        h.sdepth = read(fid, Int32)
        h.gdel   = read(fid, Int32)
        h.sdel   = read(fid, Int32)
        h.swdep  = read(fid, Int32)
        h.gwdep  = read(fid, Int32)
        h.scalel = read(fid, Int16)
        h.scalco = read(fid, Int16)
        h.sx     = read(fid, Int32)
        h.sy     = read(fid, Int32)
        h.gx     = read(fid, Int32)
        h.gy     = read(fid, Int32)
        h.counit = read(fid, Int16)
        h.wevel  = read(fid, Int16)
        h.swevel = read(fid, Int16)
        h.sut    = read(fid, Int16)
        h.gut    = read(fid, Int16)
        h.sstat  = read(fid, Int16)
        h.gstat  = read(fid, Int16)
        h.tstat  = read(fid, Int16)
        h.laga   = read(fid, Int16)
        h.lagb   = read(fid, Int16)
        h.delrt  = read(fid, Int16)
        h.muts   = read(fid, Int16)
        h.mute   = read(fid, Int16)
        h.ns     = read(fid, Int16)
        h.dt     = read(fid, Int16)
        h.gain   = read(fid, Int16)
        h.igc    = read(fid, Int16)
        h.igi    = read(fid, Int16)
        h.corr   = read(fid, Int16)
        h.sfs    = read(fid, Int16)
        h.sfe    = read(fid, Int16)
        h.slen   = read(fid, Int16)
        h.styp   = read(fid, Int16)
        h.stas   = read(fid, Int16)
        h.stae   = read(fid, Int16)
        h.tatyp  = read(fid, Int16)
        h.afilf  = read(fid, Int16)
        h.afils  = read(fid, Int16)
        h.nofilf = read(fid, Int16)
        h.nofils = read(fid, Int16)
        h.lcf    = read(fid, Int16)
        h.hcf    = read(fid, Int16)
        h.lcs    = read(fid, Int16)
        h.hcs    = read(fid, Int16)
        h.year   = read(fid, Int16)
        h.day    = read(fid, Int16)
        h.hour   = read(fid, Int16)
        h.minute = read(fid, Int16)
        h.sec    = read(fid, Int16)
        h.timbas = read(fid, Int16)
        h.trwf   = read(fid, Int16)
        h.grnors = read(fid, Int16)
        h.grnofr = read(fid, Int16)
        h.grnlof = read(fid, Int16)
        h.gaps   = read(fid, Int16)
        h.otrav  = read(fid, Int16)
        h.CDPx            = read(fid, Int32)
        h.CDPy            = read(fid, Int32)
        h.inlineNum       = read(fid, Int32)
        h.crosslineNum    = read(fid, Int32)
        h.shotNum         = read(fid, Int32)
        h.scalarOnShotNum = read(fid, Int16)
        h.tvmu            = read(fid, Int16)
        skip(fid, 10) # skip 10 bytes
        h.didf            = read(fid, Int16)
        h.scalart         = read(fid, Int16)
        h.stype           = read(fid, Int16)
        skip(fid, 20 + trace_data_size) # skip 20 bytes and trace data to the begining of next trace

        if swap_bytes
           for field in fieldnames(h)
               setfield!(h, field, bswap(getfield(h, field)))
           end
        end
        vector_trace_header[i] = h
	  end

    close(fid)
    return vector_trace_header

end

"""
   read the segy file, return FileHeader, vector of TraceHeader and matrix of data
"""
function read_segy_file(path::String; swap_bytes=true, data_format="NULL")

  println("===================================================================\n")
  println("=======================Print text header===========================\n")
  text_header = print_text_header(path)

  file_header = read_file_header(path; swap_bytes=swap_bytes)
  println("===================================================================\n")
  println("=======================Print file header===========================\n")
  print(file_header)

  if file_header.netfh == -1
     error("Does not support variable text file header")
  elseif file_header.netfh == 0
     file_hsize = 3600
  elseif file_header.netfh > 0
     file_hsize = 400 + 3200 * file_header.netfh
  end

  # compute the number of traces in this file
  file_size = filesize(path)
  num_samples = Int64(file_header.ns)
  trace_data_size = num_samples * 4
  num_traces  = round(Int64, (file_size - file_hsize) / (240 + trace_data_size))
  @printf("number of traces in this file                               : %6d\n", num_traces)

  # determine the data format
  if data_format == "NULL"
     if file_header.fmtc == 1
        data_format = "IBM"
        @printf("data format determined from file header                  : IBM\n")
     elseif file_header.fmtc == 5
        data_format = "IEEE"
        @printf("data format determined from file header                  : IEEE\n")
     end
  else
     @printf("data format determined from keyword input                : %s\n", data_format)
  end

  fid = open(path, "r")
  seek(fid, file_hsize)

  trace_header = Vector{TraceHeader}(num_traces)
  data = zeros(Float32, num_samples, num_traces)
  for i = 1 : num_traces

    h = init_TraceHeader()
    h.tracl  = read(fid, Int32)
    h.tracr  = read(fid, Int32)
    h.fldr   = read(fid, Int32)
    h.tracf  = read(fid, Int32)
    h.ep     = read(fid, Int32)
    h.cdp    = read(fid, Int32)
    h.cdpt   = read(fid, Int32)
    h.trid   = read(fid, Int16)
    h.nva    = read(fid, Int16)
    h.nhs    = read(fid, Int16)
    h.duse   = read(fid, Int16)
    h.offset = read(fid, Int32)
    h.gelev  = read(fid, Int32)
    h.selev  = read(fid, Int32)
    h.sdepth = read(fid, Int32)
    h.gdel   = read(fid, Int32)
    h.sdel   = read(fid, Int32)
    h.swdep  = read(fid, Int32)
    h.gwdep  = read(fid, Int32)
    h.scalel = read(fid, Int16)
    h.scalco = read(fid, Int16)
    h.sx     = read(fid, Int32)
    h.sy     = read(fid, Int32)
    h.gx     = read(fid, Int32)
    h.gy     = read(fid, Int32)
    h.counit = read(fid, Int16)
    h.wevel  = read(fid, Int16)
    h.swevel = read(fid, Int16)
    h.sut    = read(fid, Int16)
    h.gut    = read(fid, Int16)
    h.sstat  = read(fid, Int16)
    h.gstat  = read(fid, Int16)
    h.tstat  = read(fid, Int16)
    h.laga   = read(fid, Int16)
    h.lagb   = read(fid, Int16)
    h.delrt  = read(fid, Int16)
    h.muts   = read(fid, Int16)
    h.mute   = read(fid, Int16)
    h.ns     = read(fid, Int16)
    h.dt     = read(fid, Int16)
    h.gain   = read(fid, Int16)
    h.igc    = read(fid, Int16)
    h.igi    = read(fid, Int16)
    h.corr   = read(fid, Int16)
    h.sfs    = read(fid, Int16)
    h.sfe    = read(fid, Int16)
    h.slen   = read(fid, Int16)
    h.styp   = read(fid, Int16)
    h.stas   = read(fid, Int16)
    h.stae   = read(fid, Int16)
    h.tatyp  = read(fid, Int16)
    h.afilf  = read(fid, Int16)
    h.afils  = read(fid, Int16)
    h.nofilf = read(fid, Int16)
    h.nofils = read(fid, Int16)
    h.lcf    = read(fid, Int16)
    h.hcf    = read(fid, Int16)
    h.lcs    = read(fid, Int16)
    h.hcs    = read(fid, Int16)
    h.year   = read(fid, Int16)
    h.day    = read(fid, Int16)
    h.hour   = read(fid, Int16)
    h.minute = read(fid, Int16)
    h.sec    = read(fid, Int16)
    h.timbas = read(fid, Int16)
    h.trwf   = read(fid, Int16)
    h.grnors = read(fid, Int16)
    h.grnofr = read(fid, Int16)
    h.grnlof = read(fid, Int16)
    h.gaps   = read(fid, Int16)
    h.otrav  = read(fid, Int16)
    h.CDPx            = read(fid, Int32)
    h.CDPy            = read(fid, Int32)
    h.inlineNum       = read(fid, Int32)
    h.crosslineNum    = read(fid, Int32)
    h.shotNum         = read(fid, Int32)
    h.scalarOnShotNum = read(fid, Int16)
    h.tvmu            = read(fid, Int16)
    skip(fid, 10) # skip 10 bytes
    h.didf            = read(fid, Int16)
    h.scalart         = read(fid, Int16)
    h.stype           = read(fid, Int16)
    skip(fid, 20) # skip 20 bytes to the begining of trace data

    if swap_bytes
       for field in fieldnames(h)
           setfield!(h, field, bswap(getfield(h, field)))
       end
    end
    trace_header[i] = h

    if data_format == "IEEE"
       tmp = read(fid, Float32, num_samples)
		elseif data_format == "IBM"
	 		 tmp = read(fid, IBMFloat32, num_samples)
	 		 tmp = convert(Vector{Float32}, tmp)
		end
    data[:,i] = tmp
  end

  close(fid)
  return file_header, trace_header, data

end

"""
   only work when the data are sorted in shot gathers
"""
function create_acquisition_geometry(traces_header::Vector{TraceHeader})

    num_traces = length(traces_header)

    # determine the start and end trace for each shot
    trace_start = [1]
    trace_end   = [ ]
    for i  = 2 : num_traces
        if traces_header[i].sx != traces_header[i-1].sx || traces_header[i].sy != traces_header[i-1].sy
           push!(trace_end, i-1)
           push!(trace_start, i)
        end
    end
    push!(trace_end, num_traces)

    # determine the number of traces for each shot
    num_shots = length(trace_start)
    num_traces_per_shot = zeros(Int64, num_shots)
    for i = 1 : num_shots
        num_traces_per_shot[i] = trace_end[i] - trace_start[i] + 1
    end

    # get the source field_location
    source_x = zeros(Float32, num_shots)
    source_y = zeros(Float32, num_shots)
    receiver_x = Vector{Vector{Float32}}(num_shots)
    receiver_y = Vector{Vector{Float32}}(num_shots)
    for i = 1 : num_shots
        idx = trace_start[i]
        source_x[i] = traces_header[idx].sx
        source_y[i] = traces_header[idx].sy
        tmp_x = zeros(Float32, num_traces_per_shot[i])
        tmp_y = zeros(Float32, num_traces_per_shot[i])
        count = 0
        for j = trace_start[i] : trace_end[i]
            count = count + 1
            tmp_x[count] = traces_header[j].gx
            tmp_y[count] = traces_header[j].gy
        end
        receiver_x[i] = tmp_x
        receiver_y[i] = tmp_y
    end

    return source_x, source_y, receiver_x, receiver_y
end


"""
   create a lookup table to facilitate the reading of gathers(can be shot gather or CMP gathers)
"""
function create_lookup_table(path_txt::String, path::String; swap_bytes=true, data_format="NULL")

    # print file text header
    println("===================================================================\n")
    println("=======================Print text header===========================\n")
    text_header = print_text_header(path)

    # print file header
    file_header = read_file_header(path; swap_bytes=swap_bytes)
    println("===================================================================\n")
    println("=======================Print file header===========================\n")
    print(file_header)

    # compute the size of file header in terms of number of byte
    if file_header.netfh == -1
       error("Does not support variable text file header")
    elseif file_header.netfh == 0
       file_hsize = 3600
    elseif file_header.netfh > 0
       file_hsize = 400 + 3200 * file_header.netfh
    end
    @printf("the size of file header is                                  : %6d\n", file_hsize)

    # compute the number of traces in this file
    file_size = filesize(path)
    num_samples = Int64(file_header.ns)
    trace_data_size = num_samples * 4
    num_traces  = round(Int64, (file_size - file_hsize) / (240 + trace_data_size))
    @printf("number of traces in this file                               : %6d\n", num_traces)

    # determine the data format
    if data_format == "NULL"
       if file_header.fmtc == 1
          data_format = "IBM"
          @printf("data format determined from file header                     : IBM\n")
       elseif file_header.fmtc == 5
          data_format = "IEEE"
          @printf("data format determined from file header                     : IEEE\n")
       end
    else
       @printf("data format determined from keyword input                   : %s\n", data_format)
    end

    # skip text header and file header
    fid = open(path, "r")
    seek(fid, file_hsize)

    # read the header of all traces
    trace_header = Vector{TraceHeader}(num_traces)
    for i = 1 : num_traces
        h = init_TraceHeader()
        h.tracl  = read(fid, Int32)
        h.tracr  = read(fid, Int32)
        h.fldr   = read(fid, Int32)
        h.tracf  = read(fid, Int32)
        h.ep     = read(fid, Int32)
        h.cdp    = read(fid, Int32)
        h.cdpt   = read(fid, Int32)
        h.trid   = read(fid, Int16)
        h.nva    = read(fid, Int16)
        h.nhs    = read(fid, Int16)
        h.duse   = read(fid, Int16)
        h.offset = read(fid, Int32)
        h.gelev  = read(fid, Int32)
        h.selev  = read(fid, Int32)
        h.sdepth = read(fid, Int32)
        h.gdel   = read(fid, Int32)
        h.sdel   = read(fid, Int32)
        h.swdep  = read(fid, Int32)
        h.gwdep  = read(fid, Int32)
        h.scalel = read(fid, Int16)
        h.scalco = read(fid, Int16)
        h.sx     = read(fid, Int32)
        h.sy     = read(fid, Int32)
        h.gx     = read(fid, Int32)
        h.gy     = read(fid, Int32)
        h.counit = read(fid, Int16)
        h.wevel  = read(fid, Int16)
        h.swevel = read(fid, Int16)
        h.sut    = read(fid, Int16)
        h.gut    = read(fid, Int16)
        h.sstat  = read(fid, Int16)
        h.gstat  = read(fid, Int16)
        h.tstat  = read(fid, Int16)
        h.laga   = read(fid, Int16)
        h.lagb   = read(fid, Int16)
        h.delrt  = read(fid, Int16)
        h.muts   = read(fid, Int16)
        h.mute   = read(fid, Int16)
        h.ns     = read(fid, Int16)
        h.dt     = read(fid, Int16)
        h.gain   = read(fid, Int16)
        h.igc    = read(fid, Int16)
        h.igi    = read(fid, Int16)
        h.corr   = read(fid, Int16)
        h.sfs    = read(fid, Int16)
        h.sfe    = read(fid, Int16)
        h.slen   = read(fid, Int16)
        h.styp   = read(fid, Int16)
        h.stas   = read(fid, Int16)
        h.stae   = read(fid, Int16)
        h.tatyp  = read(fid, Int16)
        h.afilf  = read(fid, Int16)
        h.afils  = read(fid, Int16)
        h.nofilf = read(fid, Int16)
        h.nofils = read(fid, Int16)
        h.lcf    = read(fid, Int16)
        h.hcf    = read(fid, Int16)
        h.lcs    = read(fid, Int16)
        h.hcs    = read(fid, Int16)
        h.year   = read(fid, Int16)
        h.day    = read(fid, Int16)
        h.hour   = read(fid, Int16)
        h.minute = read(fid, Int16)
        h.sec    = read(fid, Int16)
        h.timbas = read(fid, Int16)
        h.trwf   = read(fid, Int16)
        h.grnors = read(fid, Int16)
        h.grnofr = read(fid, Int16)
        h.grnlof = read(fid, Int16)
        h.gaps   = read(fid, Int16)
        h.otrav  = read(fid, Int16)
        h.CDPx            = read(fid, Int32)
        h.CDPy            = read(fid, Int32)
        h.inlineNum       = read(fid, Int32)
        h.crosslineNum    = read(fid, Int32)
        h.shotNum         = read(fid, Int32)
        h.scalarOnShotNum = read(fid, Int16)
        h.tvmu            = read(fid, Int16)
        skip(fid, 10) # skip 10 bytes
        h.didf            = read(fid, Int16)
        h.scalart         = read(fid, Int16)
        h.stype           = read(fid, Int16)
        # skip to the begining of the header of next trace bytes
        skip(fid, 20 + trace_data_size)

        # swap byte
        if swap_bytes
           for field in fieldnames(h)
               setfield!(h, field, bswap(getfield(h, field)))
           end
        end

        # include in the vector
        trace_header[i] = h

    end
    close(fid)

    # count the number of shot gathers
    trace_start = [1]
    trace_end   = [ ]
    for i  = 2 : num_traces
        if trace_header[i].sx != trace_header[i-1].sx || trace_header[i].sy != trace_header[i-1].sy
           push!(trace_end, i-1)
           push!(trace_start, i)
        end
    end
    push!(trace_end, num_traces)
    num_shots = length(trace_start)

    # form a shot gather loopup table
    lookup = hcat(collect(1:num_shots), trace_start, trace_end)

    # create the lookup table into a text file
    fid = open(path_txt, "w")
    write(fid, join(["file_header_size = " "$file_hsize" "\n"]))
    write(fid, join(["num_samples      = " "$num_samples" "\n"]))
    write(fid, join(["num_traces       = " "$num_traces" "\n"]))
    write(fid, join(["num_shots        = " "$num_shots" "\n"]))
    write(fid, join(["swap_bytes       = " "$swap_bytes" "\n"]))
    write(fid, join(["data_format      = " data_format "\n"]))

    # write the trace index table
    writedlm(fid,  lookup)
    close(fid)
end


"""
   read the lookup table for shot gathers
   example: (file_header_size, num_samples, num_traces, swap_bytes, data_format, shot_idx) = read_lookup_table(path);
"""
function read_lookup_table(path::String)

    fid = open(path, "r")

    # read the file_header_size
    tmp = readline(fid); loc = rsearchindex(tmp, " = ");
    file_header_size = parse(Int32, tmp[loc+3:end])

    # read the number of samples per trace
    tmp = readline(fid); loc = rsearchindex(tmp, " = ");
    num_samples = parse(Int32, tmp[loc+3:end])

    # read the number of traces
    tmp = readline(fid); loc = rsearchindex(tmp, " = ");
    num_traces = parse(Int32, tmp[loc+3:end])

    # read the number of shots
    tmp = readline(fid); loc = rsearchindex(tmp, " = ");
    num_shots = parse(Int32, tmp[loc+3:end])

    # read the flag of byte swap
    tmp = readline(fid); loc = rsearchindex(tmp, " = ");
    swap_bytes = parse(Bool, tmp[loc+3:end])

    # read the data format of the data sampale
    tmp = readline(fid); loc = rsearchindex(tmp, " = ");
    data_format = tmp[loc+3:end]

    # close the file
    close(fid)

    # read shot start and end index
    shot_idx = readdlm(path, skipstart=6)
    shot_idx = convert(Matrix{Int32}, shot_idx)

    if num_shots != size(shot_idx, 1)
       error("the information on number of shots doesn't match")
    end

    return file_header_size, num_samples, num_traces, swap_bytes, data_format, shot_idx

end

"""
   read one shot gather based on look up table
"""
function read_one_shot(path_sgy::String, path_txt::String, idx::Integer)

    # read the look up table
    (file_header_size, num_samples, num_traces, swap_bytes, data_format, shot_idx) = read_lookup_table(path_txt)

    # the size of a trace(header plus data)
    trace_size = 240 + num_samples * 4

    # start and end trace of this shot
    trace_start   = shot_idx[idx, 2]
    trace_end     = shot_idx[idx, 3]
    num_receivers = trace_end - trace_start + 1

    # the start position
    location = file_header_size + (trace_start-1) * trace_size
    fid = open(path_sgy, "r")
    seek(fid, location)

    # read the shot gather
    trace_header = Vector{TraceHeader}(num_receivers)
    data         = zeros(Float32, num_samples, num_receivers)
    for i = 1 : num_receivers
        h = init_TraceHeader()
        h.tracl  = read(fid, Int32)
        h.tracr  = read(fid, Int32)
        h.fldr   = read(fid, Int32)
        h.tracf  = read(fid, Int32)
        h.ep     = read(fid, Int32)
        h.cdp    = read(fid, Int32)
        h.cdpt   = read(fid, Int32)
        h.trid   = read(fid, Int16)
        h.nva    = read(fid, Int16)
        h.nhs    = read(fid, Int16)
        h.duse   = read(fid, Int16)
        h.offset = read(fid, Int32)
        h.gelev  = read(fid, Int32)
        h.selev  = read(fid, Int32)
        h.sdepth = read(fid, Int32)
        h.gdel   = read(fid, Int32)
        h.sdel   = read(fid, Int32)
        h.swdep  = read(fid, Int32)
        h.gwdep  = read(fid, Int32)
        h.scalel = read(fid, Int16)
        h.scalco = read(fid, Int16)
        h.sx     = read(fid, Int32)
        h.sy     = read(fid, Int32)
        h.gx     = read(fid, Int32)
        h.gy     = read(fid, Int32)
        h.counit = read(fid, Int16)
        h.wevel  = read(fid, Int16)
        h.swevel = read(fid, Int16)
        h.sut    = read(fid, Int16)
        h.gut    = read(fid, Int16)
        h.sstat  = read(fid, Int16)
        h.gstat  = read(fid, Int16)
        h.tstat  = read(fid, Int16)
        h.laga   = read(fid, Int16)
        h.lagb   = read(fid, Int16)
        h.delrt  = read(fid, Int16)
        h.muts   = read(fid, Int16)
        h.mute   = read(fid, Int16)
        h.ns     = read(fid, Int16)
        h.dt     = read(fid, Int16)
        h.gain   = read(fid, Int16)
        h.igc    = read(fid, Int16)
        h.igi    = read(fid, Int16)
        h.corr   = read(fid, Int16)
        h.sfs    = read(fid, Int16)
        h.sfe    = read(fid, Int16)
        h.slen   = read(fid, Int16)
        h.styp   = read(fid, Int16)
        h.stas   = read(fid, Int16)
        h.stae   = read(fid, Int16)
        h.tatyp  = read(fid, Int16)
        h.afilf  = read(fid, Int16)
        h.afils  = read(fid, Int16)
        h.nofilf = read(fid, Int16)
        h.nofils = read(fid, Int16)
        h.lcf    = read(fid, Int16)
        h.hcf    = read(fid, Int16)
        h.lcs    = read(fid, Int16)
        h.hcs    = read(fid, Int16)
        h.year   = read(fid, Int16)
        h.day    = read(fid, Int16)
        h.hour   = read(fid, Int16)
        h.minute = read(fid, Int16)
        h.sec    = read(fid, Int16)
        h.timbas = read(fid, Int16)
        h.trwf   = read(fid, Int16)
        h.grnors = read(fid, Int16)
        h.grnofr = read(fid, Int16)
        h.grnlof = read(fid, Int16)
        h.gaps   = read(fid, Int16)
        h.otrav  = read(fid, Int16)
        h.CDPx            = read(fid, Int32)
        h.CDPy            = read(fid, Int32)
        h.inlineNum       = read(fid, Int32)
        h.crosslineNum    = read(fid, Int32)
        h.shotNum         = read(fid, Int32)
        h.scalarOnShotNum = read(fid, Int16)
        h.tvmu            = read(fid, Int16)
        skip(fid, 10) # skip 10 bytes
        h.didf            = read(fid, Int16)
        h.scalart         = read(fid, Int16)
        h.stype           = read(fid, Int16)
        skip(fid, 20) # skip 20 bytes to the begining of trace data

        if swap_bytes
           for field in fieldnames(h)
               setfield!(h, field, bswap(getfield(h, field)))
           end
        end
        trace_header[i] = h

        if data_format == "IEEE"
           tmp = read(fid, Float32, num_samples)
  		  elseif data_format == "IBM"
  	 		   tmp = read(fid, IBMFloat32, num_samples)
  	 		   tmp = convert(Vector{Float32}, tmp)
  		  end
        data[:,i] = tmp
    end

    return trace_header, data
end


# i = 50;
# shot = data[:, trace_start[i]:trace_end[i]]
# SeisPlot(shot, pclip=90);
#
# path = "/Users/wenlei/Desktop/for_Wenlei/data_Cynthia.sgy"
# (file_header, trace_header, data) = read_segy_file(path);
# (f, a) = amplitude_spectra(d[:,300], 0.004, 2048; fhigh_end=100);
#
# path = "/Users/wenlei/Desktop/for_Wenlei/source_wavelet.sgy"
# (file_header1, trace_header1, wavelet) = read_segy_file(path);
#
# path = "/Users/wenlei/Desktop/for_Wenlei/velocity.sgy"
# (file_header2, trace_header2, velocity) = read_segy_file(path);
#
# SeisPlot(wavelet, pclip=99);
# SeisPlot(velocity, vmin=minimum(velocity), vmax=maximum(velocity), cmap="rainbow", wbox=10, hbox=4); colorbar();
#
# (trace_start, trace_end, source_x, source_y, receiver_x, receiver_y) = create_acquisition_geometry(trace_header)
