
descr = 'Objective Function for Southern Ocean Adjoint Sensitivity Experiments';
tileno = gcmfaces_loc_tile(90,90); % For standard 13 tile output
readme = {'Objective Function as defined by space and time masks, maskC and maskT',...
    'Space mask defined by horizontal extent of winter mode water pools and mld',...
    'Vertical extent defined by max Jul-Nov mld','Time mask defined July to Nov'};
units = '1';
nt = 8*12; %Length of time mask in months
y1 = 1993; % Start year
m1 = 1; % Start month