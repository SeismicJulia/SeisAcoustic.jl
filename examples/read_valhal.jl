using SeisAcoustic

# the path to segy file
work_dir = "/home/wgao1/BP_valhall";
path_sgy = joinpath(work_dir, "obc.sgy");

# read text file header
thdr = read_text_header(path_sgy);
print_text_header(thdr)

# read file header
fhdr = read_file_header(path_sgy);

# total number of traces
tmp = filesize(path_sgy);
num_traces = convert(Int64, (tmp - 3600) / (240+sizeof(Float32)*fhdr.ns))

# read all the trace header
traces_header = extract_traces_header(path_sgy; print_interval=10000);
# for i = 1 : 500
#     sx = traces_header[i].sx
#     sy = traces_header[i].sy
#     gx = traces_header[i].gx
#     gy = traces_header[i].gy
#     println("sx: $sx, sy: $sy, gx: $gx, gy:$gy")
#     # x = traces_header[i].cdp
#     # println("$x")
# end

# write_traces_header

# read_traces_header





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
