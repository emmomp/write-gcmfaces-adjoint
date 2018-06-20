% Example calling of write_adj2netcdf and write_OF2netcdf
addpath ~/Matlab/gcmfaces  % Path to gcmfaces toolbox, available from https://github.com/gaelforget/gcmfaces
addpath ~/Matlab/m_map     % Path to m_map toolbox (if you want plots), available from facebook

% Load gcmfaces and grid
gcmfaces_global;
gloc = '~/data/orchestra/grid/'; % Edit this to your grid location
grid_load(gloc,5,'compact');

% Define variables you want to write to netcdf
variables = {'adxx_empmr','adxx_qnet','adxx_tauu','adxx_tauv','ADJtheta','ADJsalt'};
% Location of experiment outputs adxx_*12.data and ADJ*.data
expt = '/data/smurphs/emmomp/orchestra/experiments/run_ad.8yr.SOpv3.00.atl/';
% Mask name for objective function files
maskname='pvmask3_atl_mld00julnov_max_mask';

% Call write_adj2netcdf with options as defined in 'adj_netcdf_options.m'
% Type "help write_adj2netcdf" for more info
write_adj2netcdf(expt,variables,mygrid,'both')

% Summarise what has been written in example netcdfs
ncdisp([expt variables{1} '.0001.nc'])
ncdisp([expt variables{1} '.05deg.nc'])

write_OF2netcdf(expt,maskname,mygrid,'both')

% Summarise what has been written in example netcdfs
ncdisp([expt maskname '.0001.nc'])
ncdisp([expt maskname '.05deg.nc'])

% Example plots

% Read first variable first timestep data on original llc grid
var1_llc = read_nctiles([expt variables{1}],variables{1},1);
% Read first variable first timestep data on interpolated 0.5 deg grid
vinfo = ncinfo([expt variables{1} '.05deg.nc'],variables{1});
varSize = vinfo.Size; % Get variable size
stride = varSize;
stride(end) = 1; % Create stride vector to read first time step only
start = ones(size(varSize));
var1_05deg = ncread([expt variables{1} '.05deg.nc'],variables{1},start,stride);
lat_05deg = ncread([expt variables{1} '.05deg.nc'],'lat');
lon_05deg = ncread([expt variables{1} '.05deg.nc'],'lon');

%Create plot
figure
subplot(2,1,1)
m_map_gcmfaces(var1_llc,1)
title([variables{1} ' at first time step on original grid'])
subplot(2,1,2)
imagesc(lat_05deg,lon_05deg,var1_05deg')
axis xy
colorbar
title([variables{1} ' at first time step on interpolated grid'])

% Read spatial mask on original llc grid
maskc_llc = read_nctiles([expt maskname],'maskC',1);
% Read spatial mask on interpolated 0.5 deg grid
maskc_05deg = ncread([expt maskname '.05deg.nc'],'maskC');

%Create plot
figure
subplot(2,1,1)
m_map_gcmfaces(sum(maskc_llc,3),1)
title('Vertical sum of Objective Function maskC on original grid')
subplot(2,1,2)
imagesc(lat_05deg,lon_05deg,sum(maskc_05deg,3)')
axis xy
colorbar
title('Vertical sum of Objective Function maskC on interpolated grid')




