using SeisAcoustic

# the path to segy file which include 3 sourec lines 4 component OBC data
work_dir = "/home/wgao1/BP_valhall/line2";
path_sgy = joinpath(work_dir, "shot_line2.sgy");

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

# create quick lookup table for shot gathers
path_tab = joinpath(work_dir, "lookup_table.txt");
create_lookup_table(path_tab, path_sgy);

# creat acquisition geometry from shot gather
(sx, sy, gx, gy) = create_acquisition_geometry(traces_header; num_component=4);
i = 1; figure(figsize=(8,8));
scatter(gx[i], gy[i], s=1, c="g");
scatter(sx[i], sy[i], s=1, c="r");

# read one shot gather from segy file
idx = 317; (th, shot) = read_one_shot(path_sgy, path_tab, idx);

# dimensions of data cube
dt = 0.004
num_samples   = 2000
num_receivers = 2414
num_sources   = 315

# prepare one common source line shot 317 - 631
idx_l=317; idx_u=631;
(traces_header, data) = forward(path_sgy, path_tab, idx_l, idx_u, num_samples, num_receivers);

# write this source line as segy file
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

# clipping the common receiver gathers simultaneously
dc1 = copy(data1);
dc2 = copy(data2);
for i = 1 : num_receivers
    tmp1 = dc1[:,i,:]
    tmp2 = dc2[:,i,:]
    val  = quantile(abs.(vcat(vec(tmp1), vec(tmp2))), 0.90)
    dc1[:,i,:] = clamp.(tmp1, -val, val)
    dc2[:,i,:] = clamp.(tmp2, -val, val)
    println("$i")
end

# save the clipped data
path_out1 = joinpath(work_dir, "clip_line1.segy");
write_segy_file(path_out1, text_header1, file_header1, traces_header1, dc1);

path_out2 = joinpath(work_dir, "clip_line2.segy");
write_segy_file(path_out2, text_header2, file_header2, traces_header2, dc2);

# plotting
idx = 2414
figure(); imshow(dc1[:,idx,:], cmap="gray", aspect=0.2)
figure(); imshow(dc2[:,idx,:], cmap="gray", aspect=0.2)

# begin for pseudo-blending
# random time shift the second shot line
tau1 = 4.0 * rand(num_sources);
# make the time delay is the integer times of time sample interval
for i = 1 : num_sources
    tmp = round(Int64, tau1[i]/dt)
    tau1[i] = tmp * dt
end
path_rts = joinpath(work_dir, "random_time_shift1.bin")
fp = open(path_rts, "w");
write(fp, num_sources);
write(fp, tau1); close(fp);

# read the random time shift
tau2 = 4.0 * rand(num_sources);
# make the time delay is the integer times of time sample interval
for i = 1 : num_sources
    tmp = round(Int64, tau2[i]/dt)
    tau2[i] = tmp * dt
end
path_rts = joinpath(work_dir, "random_time_shift2.bin")
fp = open(path_rts, "w");
write(fp, num_sources);
write(fp, tau2); close(fp);

# read the clipped seismic data
path_clip1 = joinpath(work_dir, "clip_line1.segy");
(text_header1, file_header1, traces_header1, d_clip1) = read_segy_file(path_clip1);
d_clip1 = reshape(d_clip1, num_samples, num_receivers, num_sources);

path_clip2 = joinpath(work_dir, "clip_line2.segy");
(text_header2, file_header2, traces_header2, d_clip2) = read_segy_file(path_clip2);
d_clip2 = reshape(d_clip2, num_samples, num_receivers, num_sources);

# show the clipped seismic data
idx = 50
figure(); imshow(d_clip1[:,:,idx], cmap="gray")
figure(); imshow(d_clip1[:,:,idx], cmap="gray")

# pseudo blending, align the shots from first source line
path_tau1 = joinpath(work_dir, "random_time_shift1.bin")
fp = open(path, "r")
num_tau = rand(fp, Int64)
tau1    = zeros(Float64, num_tau)
read!(fp, tau1); close(fp)
d_blend1 = pseudo_blending(d_clip1, d_clip2, tau1);

# pseudo blending, align the shots from second source line
path_tau2 = joinpath(work_dir, "random_time_shift2.bin")
fp = open(path_tau2, "r")
num_tau = read(fp, Int64)
tau2    = zeros(Float64, num_tau)
read!(fp, tau2); close(fp)
d_blend2 = pseudo_blending(d_clip2, d_clip1, tau2);

# write the blended shot gathers into segy file
path_out = joinpath(work_dir, "blend_line1.segy");
write_segy_file(path_out, text_header1, file_header1, traces_header1, d_blend1);

path_out = joinpath(work_dir, "blend_line2.segy");
write_segy_file(path_out, text_header2, file_header2, traces_header2, d_blend2);

# write the data into binary file can be accessed for training
path_in = joinpath(work_dir, "input_train1.bin")
fp = open(path_in, "w")
write(fp, num_samples, num_receivers, num_sources);
write(fp, vec(d_blend1)); close(fp);

path_in = joinpath(work_dir, "output_train1.bin")
fp = open(path_in, "w")
write(fp, num_samples, num_receivers, num_sources);
write(fp, vec(d_clip1)); close(fp);

# for the second part of training data
path_in = joinpath(work_dir, "input_train2.bin")
fp = open(path_in, "w")
write(fp, num_samples, num_receivers, num_sources);
write(fp, vec(d_blend2)); close(fp);

path_in = joinpath(work_dir, "output_train2.bin")
fp = open(path_in, "w")
write(fp, num_samples, num_receivers, num_sources);
write(fp, vec(d_clip2)); close(fp);

# examine the figure
idx = 1111
figure(); imshow(d_clip1[:,idx,:], cmap="gray", aspect=0.2)
figure(); imshow(d_clip2[:,idx,:], cmap="gray", aspect=0.2)

# show blended common receiver gathers
idx = 100
figure(); imshow(d_clip1[:,idx,:], cmap="gray", aspect=0.2)
figure(); imshow(d_blend1[:,idx,:], cmap="gray", aspect=0.2)

figure(); imshow(d_clip2[:,idx,:], cmap="gray", aspect=0.2)
figure(); imshow(d_blend2[:,idx,:], cmap="gray", aspect=0.2)


"""
   pseudo blending function
"""
function pseudo_blending(d1::Array{Tv,3}, d2::Array{Tv,3}, tau) where {Tv<:AbstractFloat}

    # size of data cube
    (num_samples, num_receivers, num_sources) = size(d1)

    # check the size of another cucbe
    if num_samples != size(d2,1) || num_receivers != size(d2,2) || num_sources != size(d2,3) || length(tau) != num_sources
       error("size mismatch")
    end

    # copy one cube
    db = copy(d1)

    # pseudo belending of two sources
    for i = 1 : num_sources
        istart = round(Int64, tau[i]/dt) + 1

        # loop over receivers
        for j = 1 : num_receivers

            #loop over samples
            for k = istart : num_samples
                idx = k - istart + 1
                db[k,j,i] = db[k,j,i] + d2[idx,j,i]
            end
        end

        println("finished $i sources")
    end

    return db
end

"""
   agc based on rms amplitude
"""
function rms_agc(data; hwl=0.25)

    n = size(data, 2)
    dout = copy(data)
    for i = 1 : n
        dout[:,i] = agc(dout[:,i]; half_window_length=hwl)
    end
    return dout
end

"""
   t-square amplitude gaining to compensate geometrical spreading
"""
function tsquare(d; dt=0.004)

    (nt, n1) = size(d)
    dout= copy(d)

    for i = 1 : nt
        t = (i-1) * dt

        for j = 1 : n1
            dout[i,j] = dout[i,j] * t^2
        end
    end
    return dout
end

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

# plot commands via seismic unix
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

# download the figure to local machine
scp -r wgao1@saig-ml.physics.ualberta.ca:/home/wgao1/BP_valhall/line2/figure /Users/wenlei/Desktop

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
