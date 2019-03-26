using SeisAcoustic, PyPlot

# the path to segy file
work_dir = "/home/wgao1/BP_valhall/line2";
path_sgy = joinpath(work_dir, "obc.sgy");

# read text file header
text_header = read_text_header(path_sgy);
print_text_header(text_header)

# read file header
file_header = read_file_header(path_sgy);

# total number of traces
tmp = filesize(path_sgy);
num_traces = convert(Int64, (tmp - 3600) / (240+sizeof(Float32)*fhdr.ns))

# read all the trace header
traces_header = extract_traces_header(path_sgy; print_interval=10000);

# write_traces_header
path_thdr = joinpath(work_dir, "traces_header.bin");
write_traces_header(path_thdr, traces_header);

# read_traces_header
thdr = read_traces_header(path_thdr);

# create quick lookup table
path_tab = joinpath(work_dir, "lookup_table.txt");
create_lookup_table(path_tab, path_sgy);

# creat acquisition geometry from shot gather
(sx, sy, gx, gy) = create_acquisition_geometry(thdr; num_component=4);
i = 1; figure(figsize=(8,8));
scatter(gx[i], gy[i], s=1 , c="g");
scatter(sx[i], sy[i], s=1, c="r");

# read one shot gather from segy file
idx = 317;
(th1, data1) = read_one_shot(path_sgy, path_tab, idx);
idx = 946;
(th2, data2) = read_one_shot(path_sgy, path_tab, idx);

# test write segy file
path_out = joinpath(work_dir, "test.segy")
write_segy_file(path_out, text_header, fhdr, th1, data1);

(text3, fh3, th3, data3) = read_segy_file(path_out);



# double check the traces are organized in the same other for any shot gather
function foo(th1, th2)
  tmp = 0.0
  for i = 1 : 500
      gx1 = th1[i].gx; gx2 = th2[i].gx;
      gy1 = th1[i].gy; gy2 = th2[i].gy;
      println("gx: $gx1 vs $gx2, gy: $gy1 vs $gy2")
      # tmp = tmp + abs(gx1-gx2) + abs(gy1-gy2)
  end
  return tmp
end

# create two cubes for the two source line, assume the boat is driving in the opposite direction
# first source line 317-631



# second source line 946-632



# seperated into component
p1 = data[:,1:4:end]
p2 = data[:,2:4:end]
p3 = data[:,3:4:end]
p4 = data[:,4:4:end]

figure(figsize=(10,10)); imshow(p4, vmax=0.01, vmin=-0.01, cmap="seismic");

# plotting with su
tmp = "/home/wgao1/BP_valhall/line2/figure/p1.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(p1)); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/p2.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(p2)); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/p3.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(p3)); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/p4.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(p4)); close(tmpfid);

psimage < /home/wgao1/BP_valhall/line2/figure/p1.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=2414 d2num=200 f2num=200 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/p1.eps
psimage < /home/wgao1/BP_valhall/line2/figure/p2.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=2414 d2num=200 f2num=200 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/p2.eps
psimage < /home/wgao1/BP_valhall/line2/figure/p3.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=2414 d2num=200 f2num=200 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/p3.eps
psimage < /home/wgao1/BP_valhall/line2/figure/p4.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=2414 d2num=200 f2num=200 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/p4.eps

# separate into two lines
for i = 2 : length(sy)
    if abs(sy[i] - sy[i-1]) > 4000
       println("$i")
    end
end
# organize the data into two cube


# 1. gaining control
# 2. time shift in time domain

# 3. sorting in common receiver gather

# 4. testing write segy file
