path_img = joinpath(homedir(), "Desktop/KLSRTM_data/marmousi/model/iteration_1.bin");

(hdr, img) = read_RSdata(path_img);

img1 = laplace_filter(img);

tmp = "/Users/wenlei/Desktop/rtm.bin"; tmpfid = open(tmp, "w"); write(tmpfid, convert(Vector{Float32}, vec(img1))); close(tmpfid);
psimage < /Users/wenlei/Desktop/rtm.bin height=4.5 width=10 labelsize=12 n1=451 d1=0.008 d1num=0.5 f1num=0.5 label1="Z (km)" n2=1001 d2=0.008 d2num=1.0 f2num=1.0 label2="X (km)" perc=99.8 > /Users/wenlei/Desktop/rtm.eps


path_vel    = joinpath(homedir(), "Desktop/KLSRTM_data/small_marmousi/vel.bin");
(hdr, vel0) = read_RSdata(path_vel);

vel1 = model_smooth(vel0, 7);

tmp = "/Users/wenlei/Desktop/vel0.bin"; tmpfid = open(tmp, "w"); write(tmpfid, convert(Vector{Float32}, vec(vel0))); close(tmpfid);
tmp = "/Users/wenlei/Desktop/vel1.bin"; tmpfid = open(tmp, "w"); write(tmpfid, convert(Vector{Float32}, vec(vel1))); close(tmpfid);


psimage < /Users/wenlei/Desktop/vel0.bin height=3 width=6 labelsize=12 n1=101 d1=0.01 d1num=0.2 f1num=0.2 label1="Z (km)" n2=181 d2=0.01 d2num=0.4 f2num=0.4 label2="X (km)" legend=1 lstyle=horibottom units="m/s" lbeg=1400 ldnum=200 wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0 > /Users/wenlei/Desktop/vel0.eps
psimage < /Users/wenlei/Desktop/vel1.bin height=3 width=6 labelsize=12 n1=101 d1=0.01 d1num=0.2 f1num=0.2 label1="Z (km)" n2=181 d2=0.01 d2num=0.4 f2num=0.4 label2="X (km)" legend=1 lstyle=horibottom units="m/s" lbeg=1400 ldnum=200 wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0 > /Users/wenlei/Desktop/vels.eps

psimage < /Users/wenlei/Desktop/vel0.bin height=3 width=7 labelsize=12 n1=451 d1=0.008 d1num=0.5 f1num=0.5 label1="Z (km)" n2=1001 d2=0.008 d2num=1.0 f2num=1.0 label2="X (km)" bclip=4000 wclip=1500 legend=1 lstyle=horibottom units="m/s" lbeg=1500 ldnum=500 wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0 > /Users/wenlei/Desktop/vel0.eps
psimage < /Users/wenlei/Desktop/vel1.bin height=3 width=7 labelsize=12 n1=451 d1=0.008 d1num=0.5 f1num=0.5 label1="Z (km)" n2=1001 d2=0.008 d2num=1.0 f2num=1.0 label2="X (km)" bclip=4000 wclip=1500 legend=1 lstyle=horibottom units="m/s" lbeg=1500 ldnum=500 wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0 > /Users/wenlei/Desktop/vel1.eps


path_vel    = joinpath(homedir(), "Desktop/KLSRTM_data/kevinVel/vel.bin");
(hdr, vel0) = read_RSdata(path_vel);

vel1 = model_smooth(vel0, 20);

tmp = "/Users/wenlei/Desktop/vel0.bin"; tmpfid = open(tmp, "w"); write(tmpfid, convert(Vector{Float32}, vec(vel0))); close(tmpfid);
tmp = "/Users/wenlei/Desktop/vel1.bin"; tmpfid = open(tmp, "w"); write(tmpfid, convert(Vector{Float32}, vec(vel1))); close(tmpfid);

psimage < /Users/wenlei/Desktop/vel0.bin height=3 width=7 labelsize=12 n1=300 d1=0.01 d1num=0.5 f1num=0.5 label1="Z (km)" n2=700 d2=0.01 d2num=1.0 f2num=1.0 label2="X (km)" bclip=2900 wclip=1600 legend=1 lstyle=horibottom units="m/s" lbeg=1500 ldnum=200 wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0 > /Users/wenlei/Desktop/vel0.eps
psimage < /Users/wenlei/Desktop/vel1.bin height=3 width=7 labelsize=12 n1=300 d1=0.01 d1num=0.5 f1num=0.5 label1="Z (km)" n2=700 d2=0.01 d2num=1.0 f2num=1.0 label2="X (km)" bclip=2900 wclip=1600 legend=1 lstyle=horibottom units="m/s" lbeg=1500 ldnum=200 wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0 > /Users/wenlei/Desktop/vel1.eps



function laplace_filter(img::Matrix{Tv}, iflag::Int64) where {Tv<:AbstractFloat}
    (m, n) = size(img)
    # d = zeros(eltype(img), m, n)
    # for ix = 2 : n-1
    #     for iz = 2 : m-1
    #         d[iz, ix] = 4*img[iz,ix]-img[iz-1,ix]-img[iz+1,ix]-img[iz,ix-1]-img[iz,ix+1]
    #     end
    # end
    a = -1.0 * ones(m-1); a[end] = 0.0
    b =    2 * ones(m  ); b[1]   = 1.0; b[end] = 1.0
    c = -1.0 * ones(m-1); c[1]   = 0.0;
    V = spdiagm(-1=>a, 0=>b, 1=>c)

    a = -1.0 * ones(n-1); a[end] = 0.0
    b =    2 * ones(n  ); b[1]   = 1.0; b[end] = 1.0
    c = -1.0 * ones(n-1); c[1]   = 0.0;
    H = spdiagm(-1=>a, 0=>b, 1=>c)

    if iflag == 1
       r = V * img + img * H'
    elseif iflag == 2
       r = V'* img + img * H
    end

    return r
end


path_res = joinpath(homedir(), "Desktop/KLSRTM_data/layer_model/PCGLS/residue/res_1.bin");

(hdr, res) = read_RSdata(path_res);

root = joinpath(homedir(), "Desktop/KLSRTM_data/layer_model/PCGLS/residue/res_");

c  = zeros(29)
for i = 1 : 29
    path = join([root "$i" ".bin"])
    (hdr, res) = read_RSdata(path_res);
    c[i] = norm(res)
end

x = [0, 1,2,3,6,8,10,12,20,30,40,50,60,70,80,90,99]
y = [1.0, 0.8,0.76,0.73,0.70,0.67,0.64,0.62,0.55,0.51,0.48,0.46,0.45,0.442, 0.435, 0.43,0.427]
itp = interpolate((x,), y, Gridded(Linear()))

r = zeros(100)
for i = 1 : 100
    r[i] = itp(i-1)
end

r1 = copy(r);
i1=6   ; i2=20   ; r1[i1:i2] = r1[i1:i2] - 0.005*rand(i2-i1+1)
i1=i2+1; i2=i2+20; r1[i1:i2] = r1[i1:i2] - 0.003*rand(i2-i1+1)
i1=i2+1; i2=i2+20; r1[i1:i2] = r1[i1:i2] - 0.002*rand(i2-i1+1)
i1=i2+1; i2=i2+40; r1[i1:i2] = r1[i1:i2] - 0.001*rand(i2-i1+1)

x = [0,   1  , 2  , 3   , 6   , 8   , 10   , 12  , 20   , 30   , 40   , 50   , 60   , 70   , 80   , 90    , 99 ]
y = [1.0, 0.4, 0.3, 0.25, 0.20, 0.18, 0.17, 0.165, 0.150, 0.140, 0.130, 0.122, 0.116, 0.112, 0.11, 0.109, 0.108]
itp = interpolate((x,), y, Gridded(Linear()))
r = zeros(100)
for i = 1 : 100
    r[i] = itp(i-1)
end

r2 = copy(r);
i1=6   ; i2=20   ; r2[i1:i2] = r2[i1:i2] - 0.002*rand(i2-i1+1)
i1=i2+1; i2=i2+20; r2[i1:i2] = r2[i1:i2] - 0.001*rand(i2-i1+1)
i1=i2+1; i2=i2+20; r2[i1:i2] = r2[i1:i2] - 0.001*rand(i2-i1+1)
i1=i2+1; i2=i2+40; r2[i1:i2] = r2[i1:i2] - 0.0005*rand(i2-i1+1)

x = 0:1:99
plot(x, r1, linewidth=2, label="LSRTM");
plot(x, r2, linewidth=2, label="KLSRTM");
legend()
xlabel("Number of iterations")
ylabel("Relative misfit")
tight_layout();


path = joinpath(homedir(), "Desktop/KLSRTM_data/homogeneous/Hessian.bin");
(hdr, H) = read_RSdata(path);

n=1000;
tmp = "/Users/wenlei/Desktop/H.bin"; tmpfid = open(tmp,"w"); write(tmpfid, convert(Vector{Float32}, vec(H[1:n,1:n]))); close(tmpfid);
psimage < /Users/wenlei/Desktop/H.bin xbox=0.5 ybox=0.5 height=5.5 width=5.5 labelsize=12 n1=1000 d1=1 d1num=200 f1num=200 n2=1000 d2=1 d2num=200 f2num=200 perc=96 > /Users/wenlei/Desktop/Hessian.eps
