function write_adj2netcdf(expt,var,mygrid,mode)
% Writes adjoint experiment sensitivity fields to netcdf format
%INPUTS:
%        expt- Location of var files
%        var - Variable name(s) (from file name) to convert to netcdf format
%        mygrid - the mygrid object created by grid_load
%        mode - 'nctiles','interp' or 'both' for output in nctiles format,
%        interpolated to 1/2 degree grid, or both.
%        ### Define variable attributes etc in *adj_netcdf_options.m* ###


nvars = length(var);
for v = 1:nvars
    % Load options
    display('Loading Options')
    adj_netcdf_options;
    
    
    fieldname=var{v};
    longname=['Sensitivity of Objective Function to ' longn(fieldname)];
    units = [unit_of ' per ' unit_v(fieldname)];
    
    tileList=unique(convert2vector(tileno));
    
    display('Loading Variable')
    fin = [expt var{v}];
    fout = [expt var{v}];
    if strcmp(var{v}(1:4),'adxx')
        varin = rdmds2gcmfaces(fin,12);
    elseif strcmp(var{v}(1:3),'ADJ')
        varin = rdmds2gcmfaces(fin,nan);
    else
        display('Expecting variables of the type adxx or ADJ')
    end
    its = load([expt 'its_ad.txt']);
    nt = length(its);
    nd = ndims(varin.f1);
    
    time = date0_num+its*timestep/60/60/24;
    time_lag = time-date_lag0;
    
    if nd ==3
        display('Found 3D variable')
        varin = varin(:,:,1:nt);
        coords = 'lon lat tim';
    elseif nd ==4
        display('Found 4D variable')
        varin = varin(:,:,:,1:nt);
        coords = 'lon lat dep tim';
    else
        display('Expecting 3 or 4d variable input')
    end
    
    
    if strcmp(mode,'nctiles')||strcmp(mode,'both')
        
        readme{length(readme)+1}='Experiment run on LLC grid';
        
        display('===========================================================================================')
        display('***The following options have been selected, please edit adj_netcdf_options.m to change:***')
        display('===========================================================================================')
        display(sprintf('Writing var to %s.[ind].nc, %.0f tiles',fout,max(tileList)))
        display(sprintf('Variable %s:%s',fieldname,longname))
        display(sprintf('Units %s; Coords, %s; Miss val %s; Fill val %s',units,coords,missval,fillval))
        display(descr(:)')
        for r = 1:length(readme)
            display(readme{r})
        end
        
        structIn=[];
        structIn.vars.(var{v})=varin;  
        structIn.vars.lon=mygrid.XC;
        structIn.vars.lat=mygrid.YC;
        structIn.vars.tim=time;
        structIn.vars.tim_lag=time_lag;
        structIn.vars.date0=date0_num;
        structIn.vars.date_lag0=date_lag0;
        structIn.vars.area=mygrid.RAC:
        if nd==3
            structIn.vars.land=squeeze(mygrid.mskC(:,:,1));
        elseif nd==4
            structIn.vars.land=mygrid.mskC;
            structIn.vars.thic=mygrid.DRF;
            structIn.vars.dep=-mygrid.RC;         
        end
        
        structIn.descr=descr;
        
        vars=[]; nv=length(vars)+1;
        vars(nv).fldName=varin; vars(nv).longName=longname; vars(nv).units=units; nv=length(vars)+1;
        vars(nv).fldName='lon'; vars(nv).longName='longitude'; vars(nv).units='degrees_east'; nv=length(vars)+1;
        vars(nv).fldName='lat'; vars(nv).longName='latitude'; vars(nv).units='degrees_north'; nv=length(vars)+1;
        vars(nv).fldName='tim'; vars(nv).longName='time'; vars(nv).units='Days since January 0, 0000'; nv=length(vars)+1;
        vars(nv).fldName='tim_lag'; vars(nv).longName='time lag'; vars(nv).units='Lag in Days'; nv=length(vars)+1;
        vars(nv).fldName='date0'; vars(nv).longName='Start of Simulation'; vars(nv).units='Days since January 0, 0000'; nv=length(vars)+1;
        vars(nv).fldName='date_lag0'; vars(nv).longName='Start of Objective Function : Lag 0'; vars(nv).units='Days since January 0, 0000'; nv=length(vars)+1;
        vars(nv).fldName='area'; vars(nv).longName='grid cell area'; vars(nv).units='m^2'; nv=length(vars)+1;
        if nd==3
            vars(nv).fldName='land'; vars(nv).longName='land mask'; vars(nv).units='1'; nv=length(vars)+1;
        elseif nd==4
            vars(nv).fldName='land'; vars(nv).longName='land mask'; vars(nv).units='1'; nv=length(vars)+1;
            vars(nv).fldName='thic'; vars(nv).longName='grid cell thickness'; vars(nv).units='m'; nv=length(vars)+1;           
            vars(nv).fldName='dep'; vars(nv).longName='depth'; vars(nv).units='m'; nv=length(vars)+1;
        end
        
        structIn.defs=vars;
        
        struct2nctiles(expt,fout,structIn,[90 90]);
            
%         % Write main variable
%         overwrite = 1;
%         dimlist=mywrite2nctiles(fout,varin,overwrite,{'missval',missval},{'fillval',fillval},{'descr',descr},{'tileNo',tileno},{'rdm',readme},{'fldName',fieldname},{'longName',longname},{'units',units},{'coord',coords});
%         
%         % Get the right dimensions for other vars
%         dimIn2D = cell(1,length(dimlist));
%         dimInT = dimIn2D;
%         dimIn3D = dimIn2D;
%         dimInZ = dimIn2D;
%         
%         for ff=1:length(dimlist)
%             dimIn2D{ff} = {dimlist{ff}{end-1:end}};
%             dimInT{ff} = {dimlist{ff}{1}};
%             if nd == 4
%                 dimIn3D{ff} = {dimlist{ff}{end-2:end}};
%                 dimInZ{ff}={dimlist{ff}{end-2}};
%             end
%         end
%         
%         % Add other variables from grid, time
%         overwrite = 0;
%         
%         fieldname = 'lon';
%         units = 'degrees east';
%         mywrite2nctiles(fout,mygrid.XC,overwrite,{'missval',missval},{'fillval',fillval},{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'dimIn',dimIn2D});
%         
%         fieldname = 'lat';
%         units = 'degrees north';
%         mywrite2nctiles(fout,mygrid.YC,overwrite,{'missval',missval},{'fillval',fillval},{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'dimIn',dimIn2D});
%         
%         fieldname = 'tim';
%         units = 'Days since January 0, 0000';
%         mywrite2nctiles(fout,time,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'dimIn',dimInT});
%         
%         fieldname = 'tim_lag';
%         units = 'Lag in Days';
%         mywrite2nctiles(fout,time_lag,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'dimIn',dimInT});
%         
%         fieldname = 'date0';
%         longname = 'Start of Simulation';
%         units = 'Days since January 0, 0000';
%         mywrite2nctiles(fout,date0_num,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'longName',longname},{'dimIn',[]});
%         
%         fieldname = 'date_lag0';
%         longname = 'Start of Objective Function : Lag 0';
%         units = 'Days since January 0, 0000';
%         mywrite2nctiles(fout,date_lag0,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'longName',longname},{'dimIn',[]});
%         
%         area = mygrid.RAC;
%         fieldname = 'area';
%         units = 'm^2';
%         longname='grid cell area';
%         mywrite2nctiles(fout,area,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'longName',longname},{'dimIn',dimIn2D});
%         
%         mask = mygrid.mskC;
%         dep = -mygrid.RC;
%         dz = mygrid.DRF;
%         
%         if nd == 3
%             
%             fieldname = 'land';
%             units = '1';
%             longname='land mask';
%             mywrite2nctiles(fout,squeeze(mask(:,:,1)),overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'longName',longname},{'dimIn',dimIn2D});
%             
%         elseif nd == 4
%             
%             fieldname = 'land';
%             units = '1';
%             longname='land mask';
%             mywrite2nctiles(fout,mask,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'longName',longname},{'dimIn',dimIn3D});
%             
%             fieldname = 'thic';
%             units = 'm';
%             mywrite2nctiles(fout,dz,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'dimIn', dimInZ});
%             
%             fieldname = 'dep';
%             units = 'm';
%             mywrite2nctiles(fout,dep,overwrite,{'tileNo',tileno},{'fldName',fieldname},{'units',units},{'dimIn', dimInZ});
%             
        end
        
    end
    
    if strcmp(mode,'interp')||strcmp(mode,'both')
        
        adj_netcdf_options;
        fieldname=var{v};
        longname=['Sensitivity of Objective Function to ' longn(fieldname)];
        units = [unit_of ' per ' unit_v(fieldname)];
        readme{length(readme)+1}='Interpolated from LLC grid to 1/2 degree grid using gcmfaces_interp_coeffs.m';
        
        if v == 1
            lon=-179.75:0.5:179.75; lat=-89.75:0.5:89.75;
            [lat2,lon2] = meshgrid(lat,lon);
            display('Calculating Interpolation Coeffs')
            interp=gcmfaces_interp_coeffs(lon2(:),lat2(:));
            dep = -mygrid.RC;
        end
        
        display('===========================================================================================')
        display('***The following options have been selected, please edit adj_netcdf_options.m to change:***')
        display('===========================================================================================')
        display(sprintf('Writing interpolated var to %s.05deg.nc',fout))
        display(sprintf('Variable %s: %s',fieldname,longname))
        display(sprintf('Units %s; Coords, %s; Miss val %s; Fill val %s',units,coords,missval,fillval))
        display(descr(:)')
        for r = 1:length(readme)
            display(readme{r})
        end
        
        % Make Netcdf file
        fout = [fout '.05deg.nc'];
        ncid=nccreateFile(fout,'CLOBBER');
        nc_global=netcdf.getConstant('NC_GLOBAL');
        
        % Add readme
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
        
        ncputAtt(ncid,'','_FillValue',fillval);
        ncputAtt(ncid,'','missing_value',missval);
        
        ncdefDim(ncid,'lon',720)
        ncdefDim(ncid,'lat',360)
        ncdefDim(ncid,'tim',nt)
        if nd ==4
            ncdefDim(ncid,'dep',50)
        end
        
        myncdefVar(ncid,'lon','double',{'lon'});
        ncputAtt(ncid,'lon','long_name','longitude');
        ncputAtt(ncid,'lon','units','degrees_east');
        
        myncdefVar(ncid,'lat','double',{'lat'});
        ncputAtt(ncid,'lat','long_name','latitude');
        ncputAtt(ncid,'lat','units','degrees_north');
        
        myncdefVar(ncid,'tim','double',{'tim'});
        ncputAtt(ncid,'tim','long_name','time');
        ncputAtt(ncid,'tim','units','Days since January 0, 0000');
        
        if nd == 4
            myncdefVar(ncid,'dep','double',{'dep'});
            ncputAtt(ncid,'dep','long_name','depth');
            ncputAtt(ncid,'dep','units','m');
            
            myncdefVar(ncid,'thic','double',{'dep'});
            ncputAtt(ncid,'thic','long_name','cell thickness');
            ncputAtt(ncid,'thic','units','m');
        end
        
        ncclose(ncid);
        
        ncid=ncopen(fout,'write');
        
        myncputvar(ncid,'lon',lon)
        myncputvar(ncid,'lat',lat)
        myncputvar(ncid,'tim',time)
        if nd == 4
            myncputvar(ncid,'dep',dep)
            myncputvar(ncid,'dep',dz)
        end
        
        display('Writing time variables')
        
        fieldname = 'tim_lag';
        units = 'Lag in Days';
        longname = 'time_lag';
        coords = {'tim'};
        netcdf.reDef(ncid);
        myncdefVar(ncid,fieldname,'double',coords);
        ncputAtt(ncid,fieldname,'long_name',longname);
        ncputAtt(ncid,fieldname,'units',units);
        netcdf.endDef(ncid);
        myncputvar(ncid,fieldname,time_lag)
        
        fieldname = 'date0';
        longname = 'Start of Simulation';
        units = 'Days since January 0, 0000';
        coords = [];
        netcdf.reDef(ncid);
        myncdefVar(ncid,fieldname,'double',coords);
        ncputAtt(ncid,fieldname,'long_name',longname);
        ncputAtt(ncid,fieldname,'units',units);
        netcdf.endDef(ncid);
        myncputvar(ncid,fieldname,date0_num)
        
        fieldname = 'date_lag0';
        longname = 'Start of Objective Function : Lag 0';
        units = 'Days since January 0, 0000';
        coords = [];
        netcdf.reDef(ncid);
        myncdefVar(ncid,fieldname,'double',coords);
        ncputAtt(ncid,fieldname,'long_name',longname);
        ncputAtt(ncid,fieldname,'units',units);
        netcdf.endDef(ncid);
        myncputvar(ncid,fieldname,date_lag0)
        
        display('Writing main variable')
        
        fieldname=var{v};
        longname=['Sensitivity of Objective Function to ' longn(fieldname)];
        units = [unit_of ' per ' unit_v(fieldname)];
        
        netcdf.reDef(ncid);
        if nd == 3
            coords={'lon','lat','tim'};
        else
            coords={'lon','lat','dep','tim'};
        end
        myncdefVar(ncid,fieldname,'double',coords);
        ncputAtt(ncid,fieldname,'long_name',longname);
        ncputAtt(ncid,fieldname,'units',units);
        netcdf.endDef(ncid);
        
        vv = netcdf.inqVarID(ncid,fieldname);
        
        for t = 1:nt
            if nd == 3
                var_t1 = convert2vector(squeeze(varin(:,:,t)));
            else
                var_t1 = convert2vector(squeeze(varin(:,:,:,t)));
            end
            var_t0 = 1*~isnan(var_t1);
            var_t1=interp.SPM*var_t1;
            var_t0=interp.SPM*var_t0;
            if nd == 3
                varout = reshape(var_t1./var_t0,[720 360]);
                varout(isnan(varout))=0;
                netcdf.putVar(ncid,vv,[0 0 t-1],[720 360 1],varout);
            else
                varout = reshape(var_t1./var_t0,[720 360 50]);
                varout(isnan(varout))=0;
                netcdf.putVar(ncid,vv,[0 0 0 t-1],[720 360 50 1],varout);
            end
            
        end
        
        ncclose(ncid);
    end
end