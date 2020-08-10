path = joinpath(homedir(), "Desktop/traces_header_train.bin");
traces_header1 = read_traces_header(path);

path = joinpath(homedir(), "Desktop/traces_header_test.bin");
traces_header2 = read_traces_header(path);

# creat acquisition geometry from shot gather
(sx1, sy1, gx1, gy1) = create_acquisition_geometry(traces_header1; num_component=4);
(sx2, sy2, gx2, gy2) = create_acquisition_geometry(traces_header2; num_component=4);

i = 1; figure(figsize=(6,6));
scatter(gx2[i]/1000, gy2[i]/1000, s=1, c="g", label="receivers");
scatter(sx2[1:630]/1000, sy2[1:630]/1000, s=1, c="b", label="training");
# testing data set

# scatter(sx2[1:605]/1000, sy2[1:605]/1000, s=1, c="r", label="testing");
# legend();
xlabel("UTM easting (km)")
ylabel("UTM northing (km)")

gx = gx2[1]; gy = gy2[1];

il=1; iu=70;
scatter(gx[il:iu], gy[il:iu], s=1);
num_line = zeros(1)
il = []; push!(il,1)
iu = [];
for i = 2 : length(gy)
    if gy[i] - gy[i-1] > 3000
       num_line[1] = num_line[1] + 1
       push!(iu, i-1)
       push!(il, i)
    end
end
push!(iu, length(gy))

for i = 1 : 13
scatter(gx[il[i]:iu[i]], gy[il[i]:iu[i]],s=1)
end

il=1; iu=50;
scatter(gx[il:iu]/1000, gy[il:iu]/1000, s=1);
path = "/Users/wenlei/Desktop/gx.bin"
fid = open(path, "w")
write(fid, )


i = 1; figure(figsize=(8,8));
scatter(gx[i]/1000, gy[i]/1000, s=1, c="g", label="receivers");
scatter(sx[1:630]/1000, sy[1:630]/1000, s=1, c="b", label="training");
# testing data set

# scatter(sx2[1:630]/1000, sy2[1:630]/1000, s=1, c="r", label="testing");
legend();
xlabel("UTM easting (km)")
ylabel("UTM northing (km)")


# read the trace header
work_dir = "/home/wgao1/BP_valhall/line1";
path_sgy = joinpath(work_dir, "line_666_668.segy");

# read text file header
text_header = read_text_header(path_sgy);
print_text_header(text_header)

# read file header
file_header = read_file_header(path_sgy);

# total number of traces
tmp = filesize(path_sgy);
num_traces = convert(Int64, (tmp - 3600) / (240+sizeof(Float32) * file_header.ns))

# read all the trace header
traces_header = extract_traces_header(path_sgy; print_interval=100000);

# create quick lookup table for shot gathers
path_header = joinpath(work_dir, "traces_header_test.bin");
write_traces_header(path_header, traces_header);
