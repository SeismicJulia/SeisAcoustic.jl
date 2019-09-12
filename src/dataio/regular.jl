"""
   the header for regularly sampled data, the length of each string field
is 80 and support maximum of 9 dimensions.
"""
struct RegularSampleHeader{Ti<:Int64, Tv<:Float64}
    #
    n1 :: Ti
    n2 :: Ti
    n3 :: Ti
    n4 :: Ti
    n5 :: Ti
    n6 :: Ti
    n7 :: Ti
    n8 :: Ti
    n9 :: Ti
    #
    o1 :: Tv
    o2 :: Tv
    o3 :: Tv
    o4 :: Tv
    o5 :: Tv
    o6 :: Tv
    o7 :: Tv
    o8 :: Tv
    o9 :: Tv
    #
    d1 :: Tv
    d2 :: Tv
    d3 :: Tv
    d4 :: Tv
    d5 :: Tv
    d6 :: Tv
    d7 :: Tv
    d8 :: Tv
    d9 :: Tv
    #
    label1 :: String
    label2 :: String
    label3 :: String
    label4 :: String
    label5 :: String
    label6 :: String
    label7 :: String
    label8 :: String
    label9 :: String
    #
    unit1  :: String
    unit2  :: String
    unit3  :: String
    unit4  :: String
    unit5  :: String
    unit6  :: String
    unit7  :: String
    unit8  :: String
    unit9  :: String
    #
    title  :: String
    data_format :: DataType
end

"""
   the byte location of each field of the header of Uniform sampled data
"""
const field_location = Dict(
    # domension
    :n1 => 0,
    :n2 => 8,
    :n3 => 16,
    :n4 => 24,
    :n5 => 32,
    :n6 => 40,
    :n7 => 48,
    :n8 => 56,
    :n9 => 64,
    # origin
    :o1 => 72,
    :o2 => 80,
    :o3 => 88,
    :o4 => 96,
    :o5 => 104,
    :o6 => 112,
    :o7 => 120,
    :o8 => 128,
    :o9 => 136,
    # interval
    :d1 => 144,
    :d2 => 152,
    :d3 => 160,
    :d4 => 168,
    :d5 => 176,
    :d6 => 184,
    :d7 => 192,
    :d8 => 200,
    :d9 => 208,
    # label
    :label1 => 216,
    :label2 => 296,
    :label3 => 376,
    :label4 => 456,
    :label5 => 536,
    :label6 => 616,
    :label7 => 696,
    :label8 => 776,
    :label9 => 856,
    # unit
    :unit1 => 936,
    :unit2 => 1016,
    :unit3 => 1096,
    :unit4 => 1176,
    :unit5 => 1256,
    :unit6 => 1336,
    :unit7 => 1416,
    :unit8 => 1496,
    :unit9 => 1576,
    # title
    :title => 1656,
    :data_format => 1736,
    # data_format is saved as UInt8 (1 byte)
    :data  => 1737
);

"""
   constructor for the header of regularly sampled data
"""
function RegularSampleHeader(;n1=0, n2=1, n3=1, n4=1, n5=1, n6=1, n7=1, n8=1, n9=1,
         o1=0., o2=0., o3=0., o4=0., o5=0., o6=0., o7=0., o8=0., o9=0.,
         d1=0., d2=0., d3=0., d4=0., d5=0., d6=0., d7=0., d8=0., d9=0.,
         label1="", label2="", label3="", label4="", label5="", label6="", label7="", label8="", label9="",
         unit1="" , unit2="" , unit3="" , unit4="" , unit5="" , unit6="" , unit7="" , unit8="" , unit9="" ,
         title="" , data_format=Float32)

    if !(data_format <: Number)
       error("data_format must be a legal subtype of number")
    end

    hdr = RegularSampleHeader(Int64(n1)  , Int64(n2)  , Int64(n3)  , Int64(n4)  , Int64(n5)  , Int64(n6)  , Int64(n7)  , Int64(n8)  , Int64(n9)  ,
                              Float64(o1), Float64(o2), Float64(o3), Float64(o4), Float64(o5), Float64(o6), Float64(o7), Float64(o8), Float64(o9),
                              Float64(d1), Float64(d2), Float64(d3), Float64(d4), Float64(d5), Float64(d6), Float64(d7), Float64(d8), Float64(d9),
                              label1, label2, label3, label4, label5, label6, label7, label8, label9,
                              unit1 , unit2 , unit3 , unit4 , unit5 , unit6 , unit7 , unit8 , unit9 ,
                              title , data_format)

    return hdr

end

"""
   constructor for the header of regularly sampled data
"""
function RegularSampleHeader(d::Array{Tv}; n1=0, n2=1, n3=1, n4=1, n5=1, n6=1, n7=1, n8=1, n9=1,
         o1=0., o2=0., o3=0., o4=0., o5=0., o6=0., o7=0., o8=0., o9=0.,
         d1=0., d2=0., d3=0., d4=0., d5=0., d6=0., d7=0., d8=0., d9=0.,
         label1="", label2="", label3="", label4="", label5="", label6="", label7="", label8="", label9="",
         unit1="" , unit2="" , unit3="" , unit4="" , unit5="" , unit6="" , unit7="" , unit8="" , unit9="" ,
         title="" ) where {Tv<:Number}

    # get the dimensions of data
    dim  = size(d)
    ndim = length(dim)
    if ndim == 1
       n1 = dim[1]
    elseif ndim == 2
       n1 = dim[1]; n2 = dim[2];
    elseif ndim == 3
       n1 = dim[1]; n2 = dim[2]; n3 = dim[3];
    elseif ndim == 4
       n1 = dim[1]; n2 = dim[2]; n3 = dim[3]; n4 = dim[4];
    elseif ndim == 5
       n1 = dim[1]; n2 = dim[2]; n3 = dim[3]; n4 = dim[4]; n5 = dim[5];
    elseif ndim == 6
       n1 = dim[1]; n2 = dim[2]; n3 = dim[3]; n4 = dim[4]; n5 = dim[5]; n6 = dim[6];
    elseif ndim == 7
       n1 = dim[1]; n2 = dim[2]; n3 = dim[3]; n4 = dim[4]; n5 = dim[5]; n6 = dim[6]; n7 = dim[7];
    elseif ndim == 8
       n1 = dim[1]; n2 = dim[2]; n3 = dim[3]; n4 = dim[4]; n5 = dim[5]; n6 = dim[6]; n7 = dim[7]; n8 = dim[8];
    elseif ndim == 9
       n1 = dim[1]; n2 = dim[2]; n3 = dim[3]; n4 = dim[4]; n5 = dim[5]; n6 = dim[6]; n7 = dim[7]; n8 = dim[8]; n9 = dim[9];
    end

    # get data format
    data_format = eltype(d)

    hdr = RegularSampleHeader(n1, n2, n3, n4, n5, n6, n7, n8, n9,
                              o1, o2, o3, o4, o5, o6, o7, o8, o9,
                              d1, d2, d3, d4, d5, d6, d7, d8, d9,
                              label1, label2, label3, label4, label5, label6, label7, label8, label9,
                              unit1 , unit2 , unit3 , unit4 , unit5 , unit6 , unit7 , unit8 , unit9 ,
                              title , data_format)
    return hdr

end

"""
   construct data based on the header information
"""
function zeros(hdr::RegularSampleHeader)

    # get the dimensions of data
    if hdr.n9 > 1
       return zeros(hdr.data_format, hdr.n1, hdr.n2, hdr.n3, hdr.n4, hdr.n5, hdr.n6, hdr.n7, hdr.n8, hdr.n9)
    elseif hdr.n8 > 1
       return zeros(hdr.data_format, hdr.n1, hdr.n2, hdr.n3, hdr.n4, hdr.n5, hdr.n6, hdr.n7, hdr.n8)
    elseif hdr.n7 > 1
       return zeros(hdr.data_format, hdr.n1, hdr.n2, hdr.n3, hdr.n4, hdr.n5, hdr.n6, hdr.n7)
    elseif hdr.n6 > 1
       return zeros(hdr.data_format, hdr.n1, hdr.n2, hdr.n3, hdr.n4, hdr.n5, hdr.n6)
    elseif hdr.n5 > 1
       return zeros(hdr.data_format, hdr.n1, hdr.n2, hdr.n3, hdr.n4, hdr.n5)
    elseif hdr.n4 > 1
       return zeros(hdr.data_format, hdr.n1, hdr.n2, hdr.n3, hdr.n4)
    elseif hdr.n3 > 1
       return zeros(hdr.data_format, hdr.n1, hdr.n2, hdr.n3)
    elseif hdr.n2 > 1
       return zeros(hdr.data_format, hdr.n1, hdr.n2)
    elseif hdr.n1 > 1
       return zeros(hdr.data_format, hdr.n1)
    else
       error("empty header")
    end
end

"""
   overloading the show function for the header of regularly sampled data
"""
function show(io::IO, hdr::RegularSampleHeader)

    @printf("n1= %6d, o1=%6.2f, d1=%6.2f, label1=%10s, unit1=%10s\n", hdr.n1, hdr.o1, hdr.d1, hdr.label1, hdr.unit1)
    @printf("n2= %6d, o2=%6.2f, d2=%6.2f, label2=%10s, unit2=%10s\n", hdr.n2, hdr.o2, hdr.d2, hdr.label2, hdr.unit2)
    @printf("n3= %6d, o3=%6.2f, d3=%6.2f, label3=%10s, unit3=%10s\n", hdr.n3, hdr.o3, hdr.d3, hdr.label3, hdr.unit3)
    @printf("n4= %6d, o4=%6.2f, d4=%6.2f, label4=%10s, unit4=%10s\n", hdr.n4, hdr.o4, hdr.d4, hdr.label4, hdr.unit4)
    @printf("n5= %6d, o5=%6.2f, d5=%6.2f, label5=%10s, unit5=%10s\n", hdr.n5, hdr.o5, hdr.d5, hdr.label5, hdr.unit5)
    @printf("n6= %6d, o6=%6.2f, d6=%6.2f, label6=%10s, unit6=%10s\n", hdr.n6, hdr.o6, hdr.d6, hdr.label6, hdr.unit6)
    @printf("n7= %6d, o7=%6.2f, d7=%6.2f, label7=%10s, unit7=%10s\n", hdr.n7, hdr.o7, hdr.d7, hdr.label7, hdr.unit7)
    @printf("n8= %6d, o8=%6.2f, d8=%6.2f, label8=%10s, unit8=%10s\n", hdr.n8, hdr.o8, hdr.d8, hdr.label8, hdr.unit8)
    @printf("n9= %6d, o9=%6.2f, d9=%6.2f, label9=%10s, unit9=%10s\n", hdr.n9, hdr.o9, hdr.d9, hdr.label9, hdr.unit9)

    @printf("title=%20s\n", hdr.title)
    @printf("data_format=%16s", hdr.data_format)

end

"""
   write regularly sampled data to binary file
"""
function write_RSheader(path::String, hdr::RegularSampleHeader)

    #for each field, Int64, Float64,
    fid = open(path, "w")
    for n = 1 : 27
        a = getfield(hdr, n)
        write(fid, a)
    end

    # 80 byte string
    for n = 28 : 46
        str = getfield(hdr, n)
        tmp = convert(UInt8, '\n') * ones(UInt8, 80)
        for i = 1 : length(str)
            tmp[i] = convert(UInt8, str[i])
        end
        write(fid, tmp)
    end

    #data_format: 1->Float32
                # 2->Float64
                # 3->Complex64
                # 4->Complex128
                # 5->Int32
                # 6->Int64
    if hdr.data_format == Float32
       write(fid, UInt8(1))
    elseif hdr.data_format == Float64
       write(fid, UInt8(2))
    elseif hdr.data_format == Complex64
       write(fid, UInt8(3))
    elseif hdr.data_format == Complex128
       write(fid, UInt8(4))
    elseif hdr.data_format == Int32
       write(fid, UInt8(5))
    elseif hdr.data_format == Int64
       write(fid, UInt8(6))
    end
    # the size of header is 1001 byte

    flush(fid)
    return fid

end

function write_RSdata(path::String, hdr::RegularSampleHeader, d::Array{Tv}) where {Tv<:Number}

    # check the data_format
    if hdr.data_format != eltype(d)
       error("the data format is wrong")
    end

    # check the length of data
    N = hdr.n1 * hdr.n2 * hdr.n3 * hdr.n4 * hdr.n5 * hdr.n6 * hdr.n7 * hdr.n8 * hdr.n9
    if N != length(d)
       error("the length of data is wrong")
    end

    #for each field, Int64, Float64,
    fid = open(path, "w")
    for n = 1 : 27
        a = getfield(hdr, n)
        write(fid, a)
    end

    # 80 byte string
    for n = 28 : 46
        str = getfield(hdr, n)
        tmp = convert(UInt8, '\n') * ones(UInt8, 80)
        for i = 1 : length(str)
            tmp[i] = convert(UInt8, str[i])
        end
        write(fid, tmp)
    end

    #data_format: 1->Float32
                # 2->Float64
                # 3->Complex64
                # 4->Complex128
                # 5->Int32
                # 6->Int64
    if hdr.data_format == Float32
       write(fid, UInt8(1))
    elseif hdr.data_format == Float64
       write(fid, UInt8(2))
    elseif hdr.data_format == Complex{Float32}
       write(fid, UInt8(3))
    elseif hdr.data_format == Complex{Float64}
       write(fid, UInt8(4))
    elseif hdr.data_format == Int32
       write(fid, UInt8(5))
    elseif hdr.data_format == Int64
       write(fid, UInt8(6))
    end
    # the size of header is 1001 byte

    # write the data
    write(fid, vec(d))
    close(fid)
end

"""
   read the header of regularly sampled data file
"""
function read_RSheader(path::String)

    fid = open(path, "r")

    # the size of each dimension
    n1  = read(fid, Int64)
    n2  = read(fid, Int64)
    n3  = read(fid, Int64)
    n4  = read(fid, Int64)
    n5  = read(fid, Int64)
    n6  = read(fid, Int64)
    n7  = read(fid, Int64)
    n8  = read(fid, Int64)
    n9  = read(fid, Int64)

    # start index
    o1  = read(fid, Float64)
    o2  = read(fid, Float64)
    o3  = read(fid, Float64)
    o4  = read(fid, Float64)
    o5  = read(fid, Float64)
    o6  = read(fid, Float64)
    o7  = read(fid, Float64)
    o8  = read(fid, Float64)
    o9  = read(fid, Float64)

    # interval aling each dimension
    d1  = read(fid, Float64)
    d2  = read(fid, Float64)
    d3  = read(fid, Float64)
    d4  = read(fid, Float64)
    d5  = read(fid, Float64)
    d6  = read(fid, Float64)
    d7  = read(fid, Float64)
    d8  = read(fid, Float64)
    d9  = read(fid, Float64)

    # read 11 string
    indicator = convert(UInt8, '\n')
    str = Vector{String}(undef, 19)
    a = zeros(UInt8, 80 * 19)
    read!(fid, a)
    a = reshape(a, 80, 19)
    for i = 1 : 19
        k = 1
        while a[k,i] != indicator
              k = k + 1
        end
        k = k-1
        str[i] = String(a[1:k,i])
    end

    # data format
    indicator = read(fid, UInt8)
    if indicator == 1
       data_format = Float32
    elseif indicator == 2
       data_format = Float64
    elseif indicator == 3
       data_format = Complex{Float32}
    elseif indicator == 4
       data_format = Complex{Float64}
    elseif indicator == 5
       data_format = Int32
    elseif indicator == 6
       data_format = Int64
    end

    # construct a header
    return RegularSampleHeader(n1, n2, n3, n4, n5, n6, n7, n8, n9,
                               o1, o2, o3, o4, o5, o6, o7, o8, o9,
                               d1, d2, d3, d4, d5, d6, d7, d8, d9,
                               str[1] , str[2] , str[3] , str[4] , str[5] , str[6] , str[7] , str[8] , str[9] ,
                               str[10], str[11], str[12], str[13], str[14], str[15], str[16], str[17], str[18],
                               str[19], data_format)
end

"""
   read the regularly sampled data, produce two output, one is header and the other
is the data
"""
function read_RSdata(path::String)

    fid = open(path, "r")

    # the size of each dimension
    n1  = read(fid, Int64)
    n2  = read(fid, Int64)
    n3  = read(fid, Int64)
    n4  = read(fid, Int64)
    n5  = read(fid, Int64)
    n6  = read(fid, Int64)
    n7  = read(fid, Int64)
    n8  = read(fid, Int64)
    n9  = read(fid, Int64)

    # start index
    o1  = read(fid, Float64)
    o2  = read(fid, Float64)
    o3  = read(fid, Float64)
    o4  = read(fid, Float64)
    o5  = read(fid, Float64)
    o6  = read(fid, Float64)
    o7  = read(fid, Float64)
    o8  = read(fid, Float64)
    o9  = read(fid, Float64)

    # interval aling each dimension
    d1  = read(fid, Float64)
    d2  = read(fid, Float64)
    d3  = read(fid, Float64)
    d4  = read(fid, Float64)
    d5  = read(fid, Float64)
    d6  = read(fid, Float64)
    d7  = read(fid, Float64)
    d8  = read(fid, Float64)
    d9  = read(fid, Float64)

    # read 11 string
    indicator = convert(UInt8, '\n')
    str = Vector{String}(undef, 19)
    a = zeros(UInt8, 80 * 19)
    read!(fid, a)
    a = reshape(a, 80, 19)
    for i = 1 : 19
        k = 1
        while a[k,i] != indicator
              k = k + 1
        end
        k = k-1
        str[i] = String(a[1:k,i])
    end

    # data format
    indicator = read(fid, UInt8)
    if indicator == 1
       data_format = Float32
    elseif indicator == 2
       data_format = Float64
    elseif indicator == 3
       data_format = Complex{Float32}
    elseif indicator == 4
       data_format = Complex{Float64}
    elseif indicator == 5
       data_format = Int32
    elseif indicator == 6
       data_format = Int64
    end

    # construct a header
    hdr =  RegularSampleHeader(n1, n2, n3, n4, n5, n6, n7, n8, n9,
                               o1, o2, o3, o4, o5, o6, o7, o8, o9,
                               d1, d2, d3, d4, d5, d6, d7, d8, d9,
                               str[1] , str[2] , str[3] , str[4] , str[5] , str[6] , str[7] , str[8] , str[9] ,
                               str[10], str[11], str[12], str[13], str[14], str[15], str[16], str[17], str[18],
                               str[19], data_format)

    N = n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8 * n9
    data = Vector{data_format}(undef, N)

    # read data
    read!(fid, data)
    if n9 != 1
       data = reshape(data, n1, n2, n3, n4, n5, n6, n7, n8, n9)
    elseif n8 != 1
       data = reshape(data, n1, n2, n3, n4, n5, n6, n7, n8)
    elseif n7 != 1
       data = reshape(data, n1, n2, n3, n4, n5, n6, n7)
    elseif n6 != 1
       data = reshape(data, n1, n2, n3, n4, n5, n6)
    elseif n5 != 1
       data = reshape(data, n1, n2, n3, n4, n5)
    elseif n4 != 1
       data = reshape(data, n1, n2, n3, n4)
    elseif n3 != 1
       data = reshape(data, n1, n2, n3)
    elseif n2 != 1
       data = reshape(data, n1, n2)
    end
    close(fid)

    return hdr, data
end

function l2norm_rsf(path::String)
    (hdr, d) = read_RSdata(path)
    return norm(d)
end
