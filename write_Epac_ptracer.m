% Example calling of write_adj2netcdf and write_OF2netcdf
addpath ~/Matlab/gcmfaces  % Path to gcmfaces toolbox, available from https://github.com/gaelforget/gcmfaces
addpath ~/Matlab/m_map     % Path to m_map toolbox (if you want plots), available from https://www.eoas.ubc.ca/~rich/map.html

% Load gcmfaces and grid
gcmfaces_global;
gloc = '~/data/orchestra/grid/'; % Edit this to your grid location
grid_load(gloc,5,'compact');

% Define variables you want to write to netcdf
variables = {'ADJptracer01'};
% Location of experiment outputs adxx_*12.data and ADJ*.data
expt = '/data/smurphs/emmomp/orchestra/experiments/run_ad.8yr.SOpv3.00.pac.ptracer/';
% Mask name for objective function files
maskname='pvmask3_Epac_mld00julnov_max_mask';

% Call write_adj2netcdf with options as defined in 'adj_netcdf_options.m'
% Type "help write_adj2netcdf" for more info
write_adj2netcdf(expt,variables,mygrid,'nctiles')

% Summarise what has been written in example netcdfs
ncdisp([expt variables{1} '.0001.nc'])

write_OF2netcdf(expt,maskname,mygrid,'nctiles')

% Summarise what has been written in example netcdfs
ncdisp([expt maskname '.0001.nc'])





