# write-gcmfaces-adjoint

A set of matlab files designed to take the output of MITgcm adjoint runs on the LLC90 grid and write to
netcdf files in either the nctiles format or 1/2 degree interpolated.

See write2netcdf_example.m for example calling of both write_adj2netcdf (which takes sensitivities and 
writes to netcdf) and write_OF2netcdf.m (which writes maskfiles that define experiment objective functions to netcdf).

adj_netcdf_options.h and OF_netcdf_options.h contain inputs that are experiment dependent and should both 
be edited for different experiments.

my*.m are slightly edited versions of gcmfaces files that are called by the write_*.m files

Requirements:
gcmfaces toolbox from https://github.com/gaelforget/gcmfaces

