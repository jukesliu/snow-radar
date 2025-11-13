% function to save an rd FMCWprofile9 object to a .mat file
function write_rd_to_netcdf(filepath)
load(filepath);
newfilename = [filepath(1:end-4) '.nc']; % target file name

PDATA = rd.PDATA; TWT = rd.TWT; % grab PDATA & TWT
TWT = TWT(size(PDATA,1)+1:end); % grab positive half of TWT

nccreate(newfilename,"PDATA","Dimensions",{"x",size(PDATA,2),"y",size(PDATA,1)}); % create PDATA variable
ncwrite(newfilename,"PDATA",PDATA'); % write PDATA

nccreate(newfilename,"TWT","Dimensions",{"y",size(TWT,2)}); % TWT
ncwrite(newfilename,"TWT",TWT); % positive half of TWT

if length(rd.xyz) > 0
    x = rd.xyz(:,1); y = rd.xyz(:,2); z = rd.xyz(:,3);
    skycal_idx = rd.S.SkycalTraces; trace_idx = rd.S.ProfileTraces;

    % write to file
    nccreate(newfilename,"x","Dimensions",{"x",length(x)}); ncwrite(newfilename,"x", x);
    nccreate(newfilename,"y","Dimensions",{"x",length(y)}); ncwrite(newfilename,"y", y);
    nccreate(newfilename,"z","Dimensions",{"x",length(z)}); ncwrite(newfilename,"z", z);
    
    if length(skycal_idx) > 0
        nccreate(newfilename,"skycal_idx","Dimensions",{"idx",length(skycal_idx)}); ncwrite(newfilename,"skycal_idx", skycal_idx);
    end
    if length(trace_idx) > 0
        nccreate(newfilename,"trace_idx","Dimensions",{"x",length(trace_idx)}); ncwrite(newfilename,"trace_idx", trace_idx);
    end
end
% save(newfilename);
