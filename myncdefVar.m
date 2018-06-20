function varid = myncdefVar(ncid,varname,xtype,dimlist)
% add a variable in a netcdf file.
% Edited 2018 by E Boland from the gcmfaces original version 

global useNativeMatlabNetcdf; 
if isempty(useNativeMatlabNetcdf); useNativeMatlabNetcdf = ~isempty(which('netcdf.open')); end;

if useNativeMatlabNetcdf;
    iDim=[];
    if isempty(dimlist)
	%error('ncdefVar error: no dimension allocated'); 
    else
    	for ii=1:length(dimlist),
        	iDim(ii)=netcdf.inqDimID(ncid,dimlist{ii});
    	end
    end
    netcdf.defVar(ncid,varname,xtype,iDim);
else%try to use old mex stuff
    % inverse the order of list dimensions
    dimlist=fliplr(dimlist);
    switch length(dimlist)
        case 1,
            eval(sprintf('ncid{''%s''}=nc%s(''%s'');',varname,xtype,dimlist{1}));
        case 2,
            eval(sprintf('ncid{''%s''}=nc%s(''%s'',''%s'');',varname,xtype,dimlist{1},dimlist{2}));
        case 3,
            eval(sprintf('ncid{''%s''}=nc%s(''%s'',''%s'',''%s'');',varname,xtype,dimlist{1},dimlist{2},dimlist{3}));
        otherwise
            error('ncdefVar: number of dimension > 3');
    end
end


