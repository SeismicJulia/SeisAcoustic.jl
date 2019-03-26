using SeisAcoustic, PyPlot

# the path to segy file
work_dir = "/home/wgao1/BP_valhall/line2";
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
for i = 1 : 512
    sx = traces_header[i].sx
    sy = traces_header[i].sy
    gx = traces_header[i].gx
    gy = traces_header[i].gy
    println("sx: $sx, sy: $sy, gx: $gx, gy:$gy")
    # x = traces_header[i].cdp
    # println("$x")
end

# write_traces_header
path_thdr = joinpath(work_dir, "traces_header.bin");
write_traces_header(path_thdr, traces_header);

# read_traces_header
thdr = read_traces_header(path_thdr);

# creat acquisition geometry from shot gather
(sx, sy, gx, gy) = create_acquisition_geometry(thdr; num_component=4);
i = 1; figure(figsize(8,8));
scatter(gx[i], gy[i], s=1 , c="g");
scatter(sx[i], sy[i], s=40, c="r");

# create quick lookup table
path_tab = joinpath(work_dir, "lookup_table.txt");
create_lookup_table(path_tab, path_sgy);

# read one shot gather from segy file
idx = 1;
(th, data) = read_one_shot(path_sgy, path_tab, idx);
figure(); imshow(d, vmax=10, vmin=-10, cmap="seismic", aspect=0.2);
