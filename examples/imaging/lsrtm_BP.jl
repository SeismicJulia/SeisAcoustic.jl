using SeisPlot, SeisAcoustic

# ==============================================================================
#             read physical model
function lsrtm_BP()

  dir_work = joinpath(homedir(), "Desktop/lsrtm");

  path_rho = joinpath(dir_work, "physical_model/rho.rsf");
  path_vel = joinpath(dir_work, "physical_model/vel.rsf");

  # read the physical model
  (hdr_rho, rho) = read_RSdata(path_rho);
  (hdr_vel, vel) = read_RSdata(path_vel);

  # # cropped model for imaging
  vel = vel[1:2:1350,650:2:3700];
  rho = rho[1:2:1350,650:2:3700];
  vel = model_smooth(vel, 15);

  rho1 = copy(rho); rho1 .= minimum(rho);
  vel1 = copy(vel); vel1 .= vel[1];


  SeisPlotTX(vel, dx=0.0125, dy=0.0125, yticks=0:2:8, xticks=0:2:19, hbox=1.3*3, wbox=3.7*3, xlabel="X (km)", ylabel="Z (km)", cmap="rainbow", vmin=minimum(vel), vmax=maximum(vel));
  cbar = colorbar(); cbar.set_label("km/s"); tight_layout(); savefig("/Users/wenlei/Desktop/vel.pdf"); close();

  SeisPlotTX(rho, dx=0.0125, dy=0.0125, yticks=0:2:8, xticks=0:2:19, hbox=1.3*3, wbox=3.7*3, xlabel="X (km)", ylabel="Z (km)", cmap="rainbow", vmin=minimum(rho), vmax=maximum(rho));
  cbar = colorbar(); cbar.set_label("g/cm^3"); tight_layout(); savefig("/Users/wenlei/Desktop/rho.pdf"); close();

  # vertical and horizontal grid size
  dz = 6.25; dx = 6.25;

  # time step size and maximum modelling length
  dt = 0.0007; tmax = 6.0;

  # top boundary condition
  free_surface = false;
  data_format  = Float32;
  order        = 5;

  # tdparams for generating observations
  fidiff_hete = TdParams(rho, vel, free_surface, dz, dx, dt, tmax;
                    data_format=data_format, order=order);

  # tdparams for removing direct arrival
  fidiff_homo = TdParams(rho1, vel1, free_surface, dz, dx, dt, tmax;
                         data_format=data_format, order=order);

  # tdparams for imaging
  fidiff = TdParams(rho1, vel, free_surface, dz, dx, dt, tmax;
                    data_format=data_format, order=order);


  # initialize a source
  # isz = 2; isx = fidiff.nx;
  # src = Source(isz, isx, fidiff; amp=100000, fdom=20, type_flag="miniphase");

  # # vector of source
  isx = collect(40 : 12 : fidiff.nx-50); isz = 2*ones(length(isx));
  src = get_multi_sources(isz, isx, fidiff; amp=100000, fdom=20, type_flag="miniphase");

  # generate observed data
  irx = collect(1: 2 : fidiff.nx);
  irz = 2 * ones(Int64, length(irx));

  dir_obs        = joinpath(dir_work, "observations");
  dir_sourceside = joinpath(dir_work, "sourceside");

  # prepare observations and sourceside wavefield
  get_reflections(dir_obs, irz, irx, src, fidiff_hete, fidiff_homo);
  get_wavefield_bound(dir_sourceside, src, fidiff);

  born_params = (irz=irz, irx=irx, dir_sourceside=dir_sourceside, fidiff=fidiff, normalization_flag=true, mute_index=10);
  (x, his) = cgls(born_approximation, dir_obs; dir_work=dir_work, op_params=born_params,
                  d_axpby=recordings_axpby!, m_axpby=image_axpby!, d_norm=recordings_norm, m_norm=l2norm_rsf);
end

lsrtm_BP();

# download the data
scp wgao1@saig-ml.physics.ualberta.ca:/home/wgao1/Desktop/lsrtm/iterations/iteration_1.rsf /Users/wenlei/Desktop/
scp wgao1@saig-ml.physics.ualberta.ca:/home/wgao1/Desktop/lsrtm/sourceside/normalization.rsf /Users/wenlei/Desktop/

(hdr, m) = read_RSdata("/Users/wenlei/Desktop/iteration_1.rsf")
(hdr, s) = read_RSdata("/Users/wenlei/Desktop/normalization.rsf")
s = reshape(s, 675, 1526);
m = reshape(m, 675, 1526);
m .= m .* s;
m1= laplace_filter(m);

p2 = copy(rho)
for i2 = 1 : size(p2,2)
    for i1 = 1 : size(p2, 1)-1
        p2[i1,i2] = (p2[i1+1,i2]-p2[i1,i2]) / (p2[i1+1,i2]+p2[i1,i2])
    end
end

SeisPlotTX(m, cmap="gray", hbox=6.75, wbox=15.2);
SeisPlotTX(s, cmap="rainbow", hbox=6.75, wbox=15.2, vmax=maximum(s), vmin=minimum(s));
SeisPlotTX(m1, cmap="gray", hbox=6.75, wbox=15.2, pclip=96);
SeisPlotTX(p2, cmap="gray", hbox=6.75, wbox=15.2);
SeisPlotTX(vel, cmap="rainbow", hbox=6.75, wbox=15.2, vmax=maximum(vel), vmin=minimum(vel));


# # ==============================================================================
# #             convert physical model of SU to RSF
# work_dir = joinpath(homedir(), "Desktop/tmp_LSRTM");
#
# path_rho = joinpath(work_dir, "physical_model/rho.su");
# path_vel = joinpath(work_dir, "physical_model/vel.su");
#
# num_samples = 1911;
# num_traces  = 5395;
# trace = zeros(Float32, num_samples);
# vel   = zeros(Float32, num_samples, num_traces);
# rho   = zeros(Float32, num_samples, num_traces);
#
# fid_vel = open(path_vel, "r");
# fid_rho = open(path_rho, "r");
# idx = [1];
# for i = 1 : num_traces
#     skip(fid_vel, 240)
#     skip(fid_rho, 240)
#
#     read!(fid_vel, trace)
#     copyto!(vel, idx[1], trace, 1, num_samples)
#
#     read!(fid_rho, trace)
#     copyto!(rho, idx[1], trace, 1, num_samples)
#     idx[1] = idx[1] + num_samples
# end
# close(fid_vel);
# close(fid_rho);
#
# # save the physical model
# hdr_rho = RegularSampleHeader(rho; o1=0., d1=6.25, label1="Z", unit1="m",
#                                    o2=0., d2=6.25, label2="X", unit2="m",
#                                    title="density model");
#
# hdr_vel = RegularSampleHeader(vel; o1=0., d1=6.25, label1="Z", unit1="m",
#                                    o2=0., d2=6.25, label2="X", unit2="m",
#                                    title="velocity model");
#
# path_rho = joinpath(work_dir, "physical_model/rho.rsf");
# path_vel = joinpath(work_dir, "physical_model/vel.rsf");
# write_RSdata(path_rho, hdr_rho, rho);
# write_RSdata(path_vel, hdr_vel, vel);
