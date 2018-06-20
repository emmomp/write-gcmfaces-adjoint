function write_OF2netcdf(expt,maskname,mygrid,mode)
% Writes mask files that define adjoint experiment objective functions to netcdf format
%INPUTS:
%        expt- Location of mask files
%        maskname - maskname to convert to netcdf format, there should be
%        two masks of the form [maskname]C and [maskname]T
%        mygrid - the mygrid object created by grid_load
%        mode - 'nctiles','interp' or 'both' for output in nctiles format,
%        interpolated to 1/2 degree grid, or both.
%        ### Define variable attributes etc in *OF_netcdf_options.m* ###

% Load options
OF_netcdf_options;

display('Reading Mask')
%Open spatial mask
maskC = read_bin([expt maskname 'C']);
%open T mask
fid = fopen([expt maskname 'T']);
maskT =fread(fid,'float32','b');
fclose(fid);
maskT = maskT(1:nt);

display(['Masks ' maskname 'C and ' maskname 'T found'])

fout = [expt maskname];

time = zeros(nt,3);
ny = ceil(nt/12);

% Generate dates using mid-month days
t=1;
for y =y1:y1+ny-1
    for m = 1:12
        e = eomday(y,m);
        mid = round(e/2);
        time(t,:) = [y,m,mid];
        t = t+1;
    end
end

if m1 ==1 
    time2 = datenum(time(1:nt,:));
else
    time2 = datenum(time(m1:nt+m1-1,:));
end

display(['Mask T runs from ' datestr(time2(1)) ' to ' datestr(time2(end))])

tileList=unique(convert2vector(tileno));

if strcmp(mode,'nctiles')||strcmp(mode,'both')
    
    display('===========================================================================================')
    display('***The following options have been selected, please edit OF_netcdf_options.m to change:***')
    display('===========================================================================================')
    display(sprintf('Writing objective function to %s.[ind].nc, %.0f tiles',fout,max(tileList)))
    display(descr(:)')
    for r = 1:length(readme)
        display(readme{r})
    end
    
    % Write C mask
    overwrite = 1;
    coords = 'lat lon dep tim';
    dimsize = [90 90 50 nt];
    fieldname = '';
    dimlist=mywrite2nctiles(fout,[],overwrite,{'descr',descr},{'fldName',fieldname},{'tileNo',tileno},{'rdm',readme},{'coord',coords},{'dimsize',dimsize});
    
    dimIn2D = cell(1,length(dimlist));
    dimInT = dimIn2D;
    dimIn3D = dimIn2D;
    dimInZ = dimIn2D;
    
    for ff=1:length(dimlist)
        dimIn2D{ff} = flipdim({dimlist{ff}{1:2}},2);
        dimInT{ff} = {dimlist{ff}{end}};
        dimIn3D{ff} = flipdim({dimlist{ff}{1:3}},2);
        dimInZ{ff}={dimlist{ff}{3}};
    end
    
    overwrite = 0;
    fieldname = 'maskC';
    longname = 'Spatial Mask on C grid';
    coords = 'lat lon dep';
    mywrite2nctiles(fout,maskC,overwrite,{'descr',descr},{'tileNo',tileno},{'rdm',readme},{'fldName',fieldname},{'longName',longname},{'coord',coords},{'dimIn',dimIn3D});
    
    overwrite = 0;
    fieldname = 'maskT';
    longname = 'Temporal Mask defined monthly';
    coords = 'tim';
    mywrite2nctiles(fout,maskT,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'longName',longname},{'coord',coords},{'dimIn',dimInT});
    
    fieldname = 'lon';
    units = 'degrees east';
    mywrite2nctiles(fout,mygrid.XC,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'dimIn',dimIn2D});
    
    fieldname = 'lat';
    units = 'degrees north';
    mywrite2nctiles(fout,mygrid.YC,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'dimIn',dimIn2D});
    
    fieldname = 'tim';
    units = 'Days since January 0, 0000';
    mywrite2nctiles(fout,time2,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'dimIn',dimInT});
    
    area = mygrid.RAC;
    fieldname = 'area';
    units = 'm^2';
    longname='grid cell area';
    mywrite2nctiles(fout,area,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'longName',longname},{'dimIn',dimIn2D});
    
    mask = mygrid.mskC;
    dep = -mygrid.RC;
    dz = mygrid.DRF;
    
    fieldname = 'land';
    units = '1';
    longname='land mask';
    mywrite2nctiles(fout,mask,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'longName',longname},{'dimIn',dimIn3D});
    
    fieldname = 'thic';
    units = 'm';
    mywrite2nctiles(fout,dz,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'dimIn', dimInZ});
    
    fieldname = 'dep';
    units = 'm';
    mywrite2nctiles(fout,dep,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'dimIn', dimInZ});
    
end

if strcmp(mode,'interp')||strcmp(mode,'both')
    
    readme{length(readme)+1}='Interpolated from LLC grid to 1/2 degree grid using gcmfaces_interp_coeffs.m';
    
    lon=-179.75:0.5:179.75; lat=-89.75:0.5:89.75;
    [lat2,lon2] = meshgrid(lat,lon);
    display('Calculating Interpolation Coeffs')
    interp=gcmfaces_interp_coeffs(lon2(:),lat2(:));
    dep = -mygrid.RC;
    
    display('===========================================================================================')
    display('***The following options have been selected, please edit OF_netcdf_options.m to change:***')
    display('===========================================================================================')
    display(sprintf('Writing interpolated objective function to %s.05deg.nc',fout))
    display(descr(:)')
    for r = 1:length(readme)
        display(readme{r})
    end    
    
    % Make Netcdf file
    fout = [fout '.05deg.nc'];
    ncid=nccreateFile(fout,'CLOBBER');
    nc_global=netcdf.getConstant('NC_GLOBAL');
    
    % Add readme and other attributes
    if ~isempty(readme);
        descr2=[descr ' -- ' readme{1}];
    else
        descr2=descr;
    end
    ncputAtt(ncid,'','description',descr2);
    for pp=2:length(readme);
        tmp1=char(pp+63);
        netcdf.putAtt(ncid,nc_global,tmp1,readme{pp});
    end;
    ncputAtt(ncid,'','date',date);
    netcdf.putAtt(ncid,nc_global,'Conventions','CF-1.6')
    
    %Define and write dimensions    
    ncdefDim(ncid,'lon',720)
    ncdefDim(ncid,'lat',360)
    ncdefDim(ncid,'tim',nt)
    ncdefDim(ncid,'dep',50)
    
    myncdefVar(ncid,'lon','double',{'lon'});
    ncputAtt(ncid,'lon','long_name','longitude');
    ncputAtt(ncid,'lon','units','degrees_east');
    
    myncdefVar(ncid,'lat','double',{'lat'});
    ncputAtt(ncid,'lat','long_name','latitude');
    ncputAtt(ncid,'lat','units','degrees_north');
    
    myncdefVar(ncid,'tim','double',{'tim'});
    ncputAtt(ncid,'tim','long_name','time');
    ncputAtt(ncid,'tim','units','Days since January 0, 0000');
    
    myncdefVar(ncid,'dep','double',{'dep'});
    ncputAtt(ncid,'dep','long_name','depth');
    ncputAtt(ncid,'dep','units','m');
    
    myncdefVar(ncid,'thic','double',{'dep'});
    ncputAtt(ncid,'thic','long_name','cell thickness');
    ncputAtt(ncid,'thic','units','m');
    
    ncclose(ncid);
    
    ncid=ncopen(fout,'write');
    
    dz = mygrid.DRF;
    
    myncputvar(ncid,'lon',lon)
    myncputvar(ncid,'lat',lat)
    myncputvar(ncid,'tim',time2)
    myncputvar(ncid,'dep',dep)
    myncputvar(ncid,'dep',dz)
    
    display('Writing main variables')
    
    units = 1;
    
    netcdf.reDef(ncid);
    fieldname = 'maskT';
    longname = 'Temporal Mask defined monthly';
    coords = {'tim'};
    myncdefVar(ncid,fieldname,'double',coords);
    ncputAtt(ncid,fieldname,'long_name',longname);
    ncputAtt(ncid,fieldname,'units',units);
    netcdf.endDef(ncid);
    ncclose(ncid);
    
    ncid=ncopen(fout,'write');
    myncputvar(ncid,'maskT',maskT)
    
    netcdf.reDef(ncid);
    fieldname = 'maskC';
    longname = 'Spatial Mask on C grid';
    coords = {'lat','lon','dep'};
    myncdefVar(ncid,fieldname,'double',coords);
    ncputAtt(ncid,fieldname,'long_name',longname);
    ncputAtt(ncid,fieldname,'units',units);
    netcdf.endDef(ncid);
    
    % Interpolate maskC
    var_t1 = convert2vector(maskC);
    var_t0 = 1*~isnan(var_t1);
    var_t1=interp.SPM*var_t1;
    var_t0=interp.SPM*var_t0;
    varout = reshape(var_t1./var_t0,[720 360 50]);
    varout(isnan(varout))=0;
    
    vv = netcdf.inqVarID(ncid,'maskC');
    netcdf.putVar(ncid,vv,varout);
    
    ncclose(ncid);
end

end