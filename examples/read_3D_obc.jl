using SeisAcoustic

# the path to segy file
work_dir = "/home/wgao1/australia_OBN";

# work_dir = "/Users/wenlei/Desktop/OBC_2D";
path_sgy = joinpath(work_dir, "hydro.segy");

# read text file header
text_header = read_text_header(path_sgy);
print_text_header(text_header);

# read file header
file_header = read_file_header(path_sgy);

# number of traces
tmp = filesize(path_sgy);
num_traces = convert(Int64, (tmp - 3600) / (240+sizeof(Float32)*file_header.ns))

# extract all the traces' header from segy file
traces_header = extract_traces_header(path_sgy; print_interval=100000);

# print the value of the fields of trace header
fields = [:sx, :sy, :gx, :gy]
coor = extract_field_trace_header(traces_header, fields; istart=1, iend=length(traces_header));

fields = [:ep]
ep = extract_field_trace_header(traces_header, fields; istart=1, iend=length(traces_header));

# extract the coordinate
function sort_CSG(coor::Matrix{Real})

    while

end

for i = 1 : 500
    sx = traces_header[i].sx
    sy = traces_header[i].sy
    gx = traces_header[i].gx
    gy = traces_header[i].gy
    println("sx: $sx, sy: $sy, gx: $gx, gy:$gy")
    # x = traces_header[i].cdp
    # println("$x")
end

(sx, sy, gx, gy) = create_acquisition_geometry(traces_header);

path_tab = joinpath(work_dir, "lookup_table.txt");
create_lookup_table(path_tab, path_sgy);

work_dir = "/Users/wenlei/Desktop/OBC_2D";
path_sgy = joinpath(work_dir, "hydro.segy");
path_tab = joinpath(work_dir, "lookup_table.txt");
(fh_size, ns, ntrace, swap_bytes, data_format, shot_idx) = read_lookup_table(path_tab);

i = 16
(trace_header, d) = read_one_shot(path_sgy, path_tab, i);
figure(); imshow(d, vmax=10, vmin=-10, cmap="seismic", aspect=0.2);



i = 1000
(shot_x, shot_y) = extract_shot_gather(sx[i], sy[i], traces_header);
figure(figsize=(8,8)); scatter(gx, gy, s=5, c="b"); scatter(sx, sy, s=5, c="r");
scatter(shot_x, shot_y, s=10, c="g"); scatter(sx[i], sy[i], s=80, c="r");
xlim(351689, 366488); ylim(7780276, 7794926);


function extract_shot_gather(sx::T, sy::T, traces_header::Vector{TraceHeader}; print_interval=100000) where {T<:Int32}

    num = length(traces_header)
    rx  = Vector{Int32}(undef, 0)
    ry  = Vector{Int32}(undef, 0)

    for i = 1 : num
        if sx == traces_header[i].sx && sy == traces_header[i].sy
           push!(rx, traces_header[i].gx)
           push!(ry, traces_header[i].gy)
        end

        if (i % print_interval) == 0
           println("$i")
        end
    end
    return rx, ry
end

# # extract unique source and receiver's location
# (sx, sy, gx, gy) = extract_sensors_location(traces_header; print_interval=100000);

# temporarily save the unique position of sources and receivers
# path_sco = joinpath(work_dir, "unique_source_location.bin");
# fp = open(path_sco, "w")
# write(fp, convert(Int32, length(sx)))
# write(fp, sx, sy); close(fp);
path_sco = "/Users/wenlei/Desktop/unique_source_location.bin"
fp = open(path_sco, "r")
num_s = read(fp, Int32)
sx = Vector{Int32}(undef, num_s); read!(fp, sx);
sy = Vector{Int32}(undef, num_s); read!(fp, sy);

path_rco = "/Users/wenlei/Desktop/unique_receiver_location.bin"
fp = open(path_rco, "r")
num_r = read(fp, Int32)
gx = Vector{Int32}(undef, num_r); read!(fp, gx);
gy = Vector{Int32}(undef, num_r); read!(fp, gy);

# create the acquisition geometry for each shot
(src_x, src_y, rec_x, rec_y) = create_acquisition_geometry(traces_header);

figure(figsize=(8,8)); scatter(gx, gy, s=5, c="b"); scatter(sx, sy, s=5, c="r");
i = 20000; scatter(rec_x[i], rec_y[i], s=10, c="g"); scatter(src_x[i], src_y[i], s=80, c="r");
xlim(351689, 366488); ylim(7780276, 7794926);


(usx1, usy1) = remove_repeat(sx, sy);
(usx2, usy2) = remove_repeat(src_x, src_y);


figure(figsize=(8,8)); scatter(usx1[1:100], usy1[1:100], s=5, c="b");
xlim(351689, 366488); ylim(7780276, 7794926);


function min_distance(sx::Vector{T}, sy::Vector{T}) where {T<:Real}

    r = Inf32
    num = length(sx)

    for i = 1 : length(sx)

    end

end

function remove_repeat(sx::Vector{T}, sy::Vector{T}; tol=0.1) where {T<:Real}

    num = length(sx)
    usx = [sx[1]]
    usy = [sy[1]]

    res = zeros(eltype(sx), 2)
    dis = Inf32

    for i = 2 : num

        # source location
        dis = Inf32
        for j = 1 : length(usx)
            res[1] = usx[j] - sx[i]
            res[2] = usy[j] - sx[i]
            tmp = norm(res)
            if tmp < dis
               dis = tmp
            end
        end
        if dis > tol
           push!(usx, sx[i])
           push!(usy, sy[i])
        end

        println("$i")
    end

    return usx, usy
end





path_rco = joinpath(work_dir, "unique_receiver_location.bin");
fp = open(path_rco, "w")
write(fp, convert(Int32, length(gx)))
write(fp, gx, gy); close(fp);


# create the unique source location
sh = RegularSampleHeader(n1=length(sx), n2=2, label1="coordinate", label2="dimension",
     title="unique source coordinate", data_format=Int32);
src_loc = hcat(sx, sy)

path_sloc = joinpath(work_dir, "unique_source_location.rsd");
write_RSdata(path_sloc, sh, hcat(sx, sy))

#create the unique receiver location
gh = RegularSampleHeader(n1=length(gx), n2=2, label1="coordinate", label2="dimension",
     title="unique receiver coordinate", data_format=Float32);
rec_loc = hcat(gx, gy)













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
