using SeisAcoustic

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
scatter(gx[i], gy[i], s=1, c="g");
scatter(sx[i], sy[i], s=1, c="r");

# read one shot gather from segy file
idx = 317;
(th1, data1) = read_one_shot(path_sgy, path_tab, idx);

# dimensions of data cube
num_samples   = Int64(file_header.ns)
num_receivers = 2414

# prepare one common source line shot 317 - 631
idx_l=317; idx_u=631;
(traces_header, data) = forward(path_sgy, path_tab, idx_l, idx_u, num_samples, num_receivers);
path_out = joinpath(work_dir, "shot_line1.segy");
write_segy_file(path_out, text_header, file_header, traces_header, data);

# prepare one common source line shot 946-632
idx_l=632; idx_u=946;
(traces_header, data) = backward(path_sgy, path_tab, idx_l, idx_u, num_samples, num_receivers);
path_out = joinpath(work_dir, "shot_line2.segy");
write_segy_file(path_out, text_header, file_header, traces_header, data);

# read the first generated segy file
path_out1 = joinpath(work_dir, "shot_line1.segy");
(text_header1, file_header1, traces_header1, data1) = read_segy_file(path_out1);
data1 = reshape(data1, num_samples, num_receivers, num_sources);

# read the second generated segy file
path_out2 = joinpath(work_dir, "shot_line2.segy");
(text_header2, file_header2, traces_header2, data2) = read_segy_file(path_out2);
data2 = reshape(data2, num_samples, num_receivers, num_sources);


idx = 500
imshow(data1[:,1:idx], vmax=0.001, vmin=-0.001, cmap="gray", aspect=0.3)
tb  = zeros(Int64, idx)
for i = 1 : idx
    tb[i] = round(Int64, traces_header1[i].offset / 1500 / 0.004)
end
plot(1:idx, tb)





function test(trace_header)

    distance  = sqrt( (trace_header.gx - trace_header.sx)^2 + (trace_header.gy - trace_header.sy)^2)

    return distance

end


function automatic_gian_control(path_in::Ts, path_out::Ts) where {Ts<:String}

    (text_header, file_header, traces_header, data) = read_segy_file(path_in)

end

num = 300
offset = zeros(Int32, num)
for i = 1 : num
    offset[i] = traces_header1[i].offset
end

din = data1[:,1:num];
function foo(offset, din)

    (nt, n1) = size(din)
    dout     = similar(din)
    for i = 1 : n1
        d = agc(offset[i], din[:,i])
        dout[:,i] = d[:]
    end
    return dout
end

dout = foo(offset, din);

"""
   automatic gaining control applied to one trace, hwl is the half window length in time,
"""
function agc(offset, trace::Vector{Tv}; half_window_length=0.5, vel=1500, dt=0.004, drms=1.0) where {Tv<:Real}

    nt   = length(trace)
    data = copy(trace)

    dtime = offset / vel
    istart= floor(Int64, dtime / dt)
    hwl   = round(Int64, half_window_length/dt)
    wl    = 2 * hwl + 1
    ratio = 1.0
    crms  = 0.0

    if istart - hwl < 1
       istart = hwl + 1
    end

    iend = nt - hwl
    if iend <= istart
       println("window size is too big")
       return data
    end

    # compute desired rms amplitude from seismic data
    if drms == 0.0
       idx  = istart
       il   = istart - hwl
       iu   = istart + hwl
       for i = il : iu
           drms = drms + trace[i]^2
       end
       drms = sqrt(drms / wl)
       crms = drms
    end

    # boost the middle part of the data
    for idx = istart+1 : iend
        crms= sqrt((crms^2 * wl - trace[idx-hwl-1]^2 + trace[idx+hwl]^2)/wl)
        ratio = drms / crms
        data[idx] = data[idx] * ratio
    end

    # boost the last hwl samples
    for idx = iend+1 : nt
        data[idx] = data[idx] * ratio
    end

    return data

end

# start from the first sample
function agc(trace::Vector{Tv}; half_window_length=0.5, dt=0.004, drms=1.0, delta=1e-16) where {Tv<:Real}

    nt   = length(trace)
    data = copy(trace)

    hwl = round(Int64, half_window_length/dt)
    wl  = 2 * hwl + 1

    if wl > nt
       println("window size is bigger than the length of input data")
       return data
    end
    # the rest of the code can handle the case when nt = wl

    # the start and end index of window center
    istart= hwl + 1
    iend  = nt  - hwl

    # set crms as the rms amplitude of the first window
    crms  = 0.0
    il = istart - hwl
    iu = istart + hwl
    for i = il : iu
        crms = crms + trace[i]^2
    end
    crms = sqrt(crms / wl)

    # the ratio between drms and crms
    ratio = drms / (crms+delta)

    # apply to the samples in the first half of window
    for i = 1 : istart
        data[i] = data[i] * ratio
    end

    # loop over from the second to the last window
    for idx = istart+1 : iend
        crms= sqrt((crms^2 * wl - trace[idx-hwl-1]^2 + trace[idx+hwl]^2)/wl)
        ratio = drms / (crms+delta)
        data[idx] = data[idx] * ratio
    end

    # boost the last hwl samples
    for idx = iend+1 : nt
        data[idx] = data[idx] * ratio
    end

    return data

end

function tsquare_gain(d; dt=0.004)
    nt  = length(d)
    dout= copy(d)
    for i = 1 : nt
        t = (i-1) * dt
        dout[i] = dout[i] * t^2
    end
    return dout
end

gx1 = zeros(Int32, num_receivers)
gy1 = zeros(Int32, num_receivers)
gx2 = zeros(Int32, num_receivers)
gy2 = zeros(Int32, num_receivers)
for i = 1 : num_receivers
    gx1[i] = traces_header1[i].gx
    gy1[i] = traces_header1[i].gy
    gx2[i] = traces_header2[i].gx
    gy2[i] = traces_header2[i].gy
end

num_sources = 315
sx1 = zeros(Int32, num_sources)
sy1 = zeros(Int32, num_sources)
sx2 = zeros(Int32, num_sources)
sy2 = zeros(Int32, num_sources)
for i = 1 : num_sources
    idx = (i-1)*num_receivers + 1
    sx1[i] = traces_header1[idx].sx
    sy1[i] = traces_header1[idx].sy
    sx2[i] = traces_header2[idx].sx
    sy2[i] = traces_header2[idx].sy
end

function forward(path_sgy, path_tab, idx_l, idx_u, num_samples, num_receivers)

    # number of shot in this source line
    num_shots = idx_u - idx_l + 1

    # allocate memory for all the traces
    data = zeros(Float32, num_samples, num_receivers, num_shots);

    # trace header
    traces_header = Vector{TraceHeader}(undef, 0);

    ncount = 1
    for idx = idx_l : idx_u

        # read one shot gather
        (th, d) = read_one_shot(path_sgy, path_tab, idx)

        # take out the fourth component
        th = th[4:4:end]
        d  = d[:,4:4:end]

        # put to the main one
        traces_header = vcat(traces_header, th)
        data[:,:,ncount] = d

        println("$ncount")
        #jump to next shot
        ncount = ncount + 1
    end

    return traces_header, data
end

function backward(path_sgy, path_tab, idx_l, idx_u, num_samples, num_receivers)

    # number of shot in this source line
    num_shots = idx_u - idx_l + 1

    # allocate memory for all the traces
    data = zeros(Float32, num_samples, num_receivers, num_shots);

    # trace header
    traces_header = Vector{TraceHeader}(undef, 0);

    ncount = 1
    for idx = idx_u : -1 : idx_l

        # read one shot gather
        (th, d) = read_one_shot(path_sgy, path_tab, idx)

        # take out the fourth component
        th = th[4:4:end]
        d  = d[:,4:4:end]

        # put to the main one
        traces_header = vcat(traces_header, th)
        data[:,:,ncount] = d

        println("$ncount")
        #jump to next shot
        ncount = ncount + 1
    end

    return traces_header, data
end


# test write segy file: checked
path_out = joinpath(work_dir, "test.segy")
write_segy_file(path_out, text_header, file_header, th1, data1);
(text3, fh3, th3, data3) = read_segy_file(path_out);







# double check the traces are organized in the same other for any shot gather
function foo(th1, th2)
  tmp = 0.0
  # for i = 1 : length(th1)
  #     # gx1 = th1[i].gx; gx2 = th2[i].gx;
  #     # gy1 = th1[i].gy; gy2 = th2[i].gy;
  #     # println("gx: $gx1 vs $gx2, gy: $gy1 vs $gy2")
  #     # tmp = tmp + abs(gx1-gx2) + abs(gy1-gy2)
  #     for field in fieldnames(TraceHeader)
  #         tmp = tmp + abs(getfield(th1[i], field) - getfield(th2[i], field))
  #     end
  # end
  for field in fieldnames(FileHeader)
      tmp = tmp + abs(getfield(th1, field) - getfield(th2, field))
  end
  return tmp
end

# create two cubes for the two source line, assume the boat is driving in the opposite direction
# first source line 317-631


# plot common shot gather
tmp = "/home/wgao1/BP_valhall/line2/figure/p1.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,:,10])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/p2.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,:,110])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/p3.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,:,210])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/p4.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,:,310])); close(tmpfid);

psimage < /home/wgao1/BP_valhall/line2/figure/p1.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=2414 d2num=200 f2num=200 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/p1.eps
psimage < /home/wgao1/BP_valhall/line2/figure/p2.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=2414 d2num=200 f2num=200 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/p2.eps
psimage < /home/wgao1/BP_valhall/line2/figure/p3.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=2414 d2num=200 f2num=200 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/p3.eps
psimage < /home/wgao1/BP_valhall/line2/figure/p4.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=2414 d2num=200 f2num=200 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/p4.eps

# plot common receiver gather
tmp = "/home/wgao1/BP_valhall/line2/figure/r1.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,200,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r2.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,400,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r3.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,600,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r4.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,800,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r5.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,1000,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r6.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,1200,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r7.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,1400,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r8.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,1600,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r9.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,1800,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r10.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,2000,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r11.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,2200,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r12.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data1[:,2400,:])); close(tmpfid);

psimage < /home/wgao1/BP_valhall/line2/figure/r1.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r1.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r2.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r2.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r3.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r3.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r4.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r4.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r5.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r5.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r6.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r6.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r7.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r7.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r8.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r8.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r9.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r9.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r10.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r10.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r11.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r11.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r12.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r12.eps

# plot the second shot line
tmp = "/home/wgao1/BP_valhall/line2/figure/r21.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data2[:,200,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r22.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data2[:,400,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r23.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data2[:,600,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r24.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data2[:,800,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r25.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data2[:,1000,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r26.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data2[:,1200,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r27.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data2[:,1400,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r28.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data2[:,1600,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r29.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data2[:,1800,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r210.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data2[:,2000,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r211.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data2[:,2200,:])); close(tmpfid);
tmp = "/home/wgao1/BP_valhall/line2/figure/r212.bin"; tmpfid = open(tmp, "w"); write(tmpfid, vec(data2[:,2400,:])); close(tmpfid);

psimage < /home/wgao1/BP_valhall/line2/figure/r21.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r21.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r22.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r22.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r23.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r23.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r24.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r24.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r25.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r25.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r26.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r26.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r27.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r27.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r28.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r28.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r29.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r29.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r210.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r210.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r211.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r211.eps
psimage < /home/wgao1/BP_valhall/line2/figure/r212.bin xbox=0.5 ybox=0.5 height=8 width=12 labelsize=15 n1=2000 label1="Time" n2=315 d2num=50 f2num=50 label2="Trace" perc=90 > /home/wgao1/BP_valhall/line2/figure/r212.eps

scp -r wgao1@saig-ml.physics.ualberta.ca:/home/wgao1/BP_valhall/line2/figure /Users/wenlei/Desktop

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
# 2. time shift in frequency domain
"""
Shifting the shots according to the fire time, double the time samples padding at head
d: 3D Array, first dimension is time, second dimension is channel, third dimension is shot
"""
function DataShiftPadHead(d::Array{Tv,3}, shift::Vector{Tv}, dt::Tv) where Tv<:Float64
    (nt, nr, ns) = size(d)
    # double the time samples
    df = vcat(zeros(Float64, nt, nr, ns), d)
    df = fft(df, 1)
    wn = collect(0:nt)*2*pi/(2*nt*dt)
    ps = zeros(Complex{Float64}, nt+1)
    for i2 = 1 : ns   # loop over shots
        for iw = 1 : nt+1
            ps[iw] = exp(-im * wn[iw] * shift[i2])
        end
        for i1 = 1 : nr
            for iw = 1 : nt+1
                df[iw,i1,i2] = df[iw,i1,i2] * ps[iw]
            end
            for iw = 2 : nt
                idx = 2*nt - iw + 2
                df[idx,i1,i2] = conj(df[iw,i1,i2])
            end
        end
    end
    ds = real(ifft(df,1))[nt+1:end,:,:]
    return ds
end

"""
Shifting the shots according to the fire time, padding zeros at tail to make the number of
time samples equal to nextpow2.
d: 3D Array, first dimension is time, second dimension is channel, third dimension is shot
"""
function DataShiftPadTail(d::Array{Tv,3}, shift::Vector{Tv}, dt::Tv) where Tv<:Float64
    (nt, nr, ns) = size(d)
    nf = nextpow2(nt); hf = floor(Int64, nf/2);
    df = vcat(d, zeros(Float64, nf-nt, nr, ns)); df = fft(df, 1);
    wn = collect(0:hf)*2*pi/(nf*dt)
    ps = zeros(Complex{Float64}, nf+1)
    for i2 = 1 : ns   # loop over shots
        for iw = 1 : hf+1
            ps[iw] = exp(-im * wn[iw] * shift[i2])
        end
        for i1 = 1 : nr
            for iw = 1 : hf+1
                df[iw,i1,i2] = df[iw,i1,i2] * ps[iw]
            end
            for iw = 2 : hf
                idx = nf - iw + 2
                df[idx,i1,i2] = conj(df[iw,i1,i2])
            end
        end
    end
    ds = real(ifft(df,1))[1:nt,:,:]
    return ds
end

"""
Shifting the shots according to the fire time, no padding is applied
d: 3D Array, first dimension is time, second dimension is channel, third dimension is shot
"""
function DataShiftNoPad(d::Array{Tv,3}, shift::Vector{Tv}, dt::Tv) where Tv<:Float64
    (nt, nr, ns) = size(d)
    nf = floor(Int64, nt/2)  # assume the number of samples is even
    df = fft(d, 1)
    wn = collect(0:nf)*2*pi/(nt*dt)
    ps = zeros(Complex{Float64}, nf+1)
    for i2 = 1 : ns   # loop over shots
        for iw = 1 : nf+1
            ps[iw] = exp(-im * wn[iw] * shift[i2])
        end
        for i1 = 1 : nr
            for iw = 1 : nf+1
                df[iw,i1,i2] = df[iw,i1,i2] * ps[iw]
            end
            for iw = 2 : nf
                idx = nt - iw + 2
                df[idx,i1,i2] = conj(df[iw,i1,i2])
            end
        end
    end
    ds = real(ifft(df,1))
    return ds
end


# 3. sorting in common receiver gather

# 4. testing write segy file
