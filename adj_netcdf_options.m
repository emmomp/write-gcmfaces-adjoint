% Options for function write_adj2netcdf 
% EDIT to provide details of relevent experiment

descr = 'Results of SO Sensitivity Experiment';
tileno = gcmfaces_loc_tile(90,90); % For standard 13 tile output
readme = {'8 year adjoint sensitivities to ECCOv4_r2 Southern Ocean 2000 Winter Mode Water Heat Content',...
    %'Objective Function Defined in file pvmask3_ind_mld00julnov_max_mask.[ind].nc'};
    'Objective Function Defined in file pvmask3_Epac_mld00julnov_max_mask.[ind].nc'};

missval = nan;
fillval = nan;
% deltat from mitgcm run
timestep = 3600.; 
% Start of Experiment
date0_num = datenum('1993-01-01 12:00:00');
% Date that is considered "lag zero"
date_lag0 = datenum('2000-07-01 12:00:00');

% Define longnames and units for various fields here
key = {'adxx_empmr','adxx_qnet','adxx_tauu','adxx_tauv','ADJtheta','ADJsalt','ADJkin','ADJdyn','ADJptracer01'};

vals = {'surface fresh water flux E-P-R','surface heat flux Qnet','surface zonal wind stress','surface meridional wind stress'...
    '3-D potential temperature','3-D salinity','3-D kinematic temperature changes at constant density',...
    '3-D dynamic temperature changes at constant salinity','3-D passive tracer'};    
longn = containers.Map(key,vals);

%unit_of = 'Degree C'; %Objective Function units
unit_of = 'Tracer concentration';
vals = {'m/s','W/m^2','N/m^2','N/m^2','Degree C','psu','Degree C','Degree C','Tracer conc.'};
unit_v = containers.Map(key,vals);