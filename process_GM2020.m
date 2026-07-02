%% FMCW data process
datadisk='/Users/jukesliu/Documents/POSTDOC/snow-radar/GM2020/FMCW_raw/';
codedisk='/Users/jukesliu/Documents/GitHub/FMCWradarGround/';
cd /Users/jukesliu/Documents/GitHub/FMCWradarGround/; % where FMCW code lives
datefolder='GM20200212/'; % select the date folder to process
year=str2double(datefolder(3:6)); month=str2double(datefolder(7:8));  day=str2double(datefolder(9:10)); 
dir([datadisk datefolder 'FMCW/']).name % print transect names
% path to the geode gps file 
gpsdatafile = [datadisk datefolder 'RADAR_GPS/gm12Feb2020_geode.txt']; % geode
% gpsdatafile = [datadisk datefolder 'RADAR_GPS/RadarGPS_26Jan2020.txt'];

transectfolder=['base_to_1C1/']; % select the transect
rd=FMCWprofile9; % make a FMCW profile object
fs = 12; % standard fontsize for all figures
velocity=2.44e10; % enter avg radar velocity for the site [cm/s]
outputfilepath = ['/Users/jukesliu/Documents/POSTDOC/snow-radar/GM2020/FMCW_processed/' transectfolder]; % path to output folder
mkdir(outputfilepath);

% processing parameters
rd.data_dir=[datadisk datefolder 'FMCW/' transectfolder];
rd.M.frange=[6 18]; % low freq osc used for this survey
rd.M.Fs=100000; % sample frequency [Hz]
rd.P.frange=[7 17]; % process away from endpoints to reduce noise
% rd.P.channel=2; % time domain signal is on this raw data channel (HH?)
% rd.P.channel=4; % 4 indicates HH+HV
rd.P.channel=2; % 5 indicates HH*HV
% rd.P.channel=10; % use 10 to indicate coherence matrix
rd.P.nfft=2^14; % number of points in FFT
% grab number of files in folder
filelist=dir(rd.data_dir);filenum=0;
for f=1:length(filelist)
    if startsWith(filelist(f).name, 'FMCW')
        filenum=filenum+1;
    end
end

% start with skycal
rd=subdivide_daq(rd); % subdivide raw daq files
skycal=rd.TDATA; % store TDATA from this skycal
rd.TDATA=[]; % remove from rd object

% process
batchsize=5; % set batchsize in number of files
if filenum > batchsize % batch process
    rd.P.files=1:batchsize; % for splitting into batches, must be less than the number of file in folder
    nbatches=round(filenum/batchsize); % grab number of batches

    % grab digit format
    if length(num2str(nbatches)) == 1
        format = '%01.f';
    elseif length(num2str(nbatches)) == 2
        format = '%02.f';
    elseif length(num2str(nbatches)) == 3
        format = '%03.f';
    end

    % now load data files
    PDATA=[]; CPUtime=[];
    for n=1:nbatches % loop through number of batches
        s1=1+batchsize*(n-1);  % batch start index
        if n==nbatches % on the last batch
            s2=filenum; % batch end index is the number of files
        else
            s2=s1+(batchsize-1); % batch end index
        end
        rd.P.files=s1:s2; % grab file indices
        rd.data_dir=[datadisk datefolder 'FMCW/' transectfolder]; % raw data location
        rd=subdivide_daq(rd); % subdivide raw daq files into batches
        % rd.TDATA=rd.TDATA-skycal; % perform sky calibration
        rd=cal_psd_radar(rd); % process for freq domain (calculate power spectral density)

        % vertical subset
        hmax_idx = 2^12; % subset index
        rd.PDATA=rd.PDATA(1:hmax_idx,:); % subset PDATA
        % plot(rd.TWT); hold on; plot(rd.TWT(8193:8193+hmax_idx-1));
        rd.TWT=rd.TWT(8193:8193+hmax_idx-1); % subset TWT from 2^13 fft
        
        % save to variable
        PDATA = rd.PDATA; CPUtime = (rd.CPUtime+7/24);
        % PDATA=[PDATA rd.PDATA]; % add next file
        % CPUtime=[CPUtime (rd.CPUtime+7/24)]; % convert to UTC (+6 for MDT / +7 for MST) and store

        % keep only traces without NaNs
        nan_cols_logical = any(~isnan(PDATA), 1); col_indices = find(nan_cols_logical);
        PDATA = PDATA(:, col_indices); 
        CPUtime = CPUtime(:, col_indices);

        % plot result
        x=1:length(rd.CPUtime);
        depth=rd.TWT*velocity/2;
        % figure(1);clf;  
        % imagesc(x,depth,PDATA,[min(PDATA,[],'all'), max(PDATA,[],'all')]); % may need to adjust cmax and cmin
        % xlabel('xidx','FontSize',fs); ylabel('Depth [cm]','FontSize',fs);
        % title(transectfolder(1:end-1),'Interpreter', 'none','FontSize',fs);
        % colorbar; ax = gca; ax.FontSize = fs;

        % now lets geolocate traces
        D=readtable(gpsdatafile, 'ReadVariableNames',false);
        Ix=find(startsWith(string(D.Var1),'$GPGGA'));
        lat=zeros(length(Ix),1);
        lon=zeros(length(Ix),1);
        UTC=zeros(length(Ix),1);
        for m=1:length(Ix)
            S=D(Ix(m),1).Var1{1};
            D2=split(S,',');
            if length(D2{3}) > 1
                lat(m)=str2num(D2{3}(1:2))+(str2num(D2{3}(3:end)))/60;
                lon(m)=-(str2num(D2{5}(1:3))+(str2num(D2{5}(4:end))/60));
            else
                lat(m)=NaN;lon(m)=NaN;
            end
            if isempty(D2{2}) % if no timestamp
                UTC(m,1)=NaN;
            else
                UTC(m,1)=datenum(year,month,day,str2num(D2{2}(1:2)),str2num(D2{2}(3:4)),str2num(D2{2}(5:end))); % ENTER THE DATE
            end
        end
        % grab unique and finite values only
        UTC = UTC(~isnan(UTC)); % remove empty timestamps
        UTC = UTC(isfinite(UTC)); % grab only finite values
        [UTCu, unique_idxs, i0] = unique(UTC);
        latu = lat(unique_idxs); lonu = lon(unique_idxs); 

        % plot geolocated figure
        [x,y,zone]=ll2utm(latu,lonu); % convert to UTM
        Rx=interp1(UTCu,x,CPUtime); % get easting for CPUtimes
        Ry=interp1(UTCu,y,CPUtime); % get northing for CPUtimes
        % figure(2); clf; % figure params
        % ax=gca; ax.XAxis.FontSize = fs; ax.YAxis.FontSize = fs; axis equal; 
        % plot(Rx,Ry,'.'); title('UTM coordinates for profile',FontSize=fs);
        % xlabel('Easting [m]',FontSize=fs); ylabel('Northing [m]','FontSize',fs);

        % write object to netcdf file
        filename = [transectfolder(1:end-1) '_' num2str(n,format) '.nc']; % grab the transect name, add batch number
        newfilename = [outputfilepath filename]; % create filepath and name
        if isfile(newfilename) % File exists, so delete it
            delete(newfilename);
        end
        
        TWT=rd.TWT; % grab PDATA & TWT
        PDATA=cast(PDATA,'double');
        % TWT = TWT(size(PDATA,1)+1:end); % grab positive half of TWT - not
        % necessary for GM

        nccreate(newfilename,"PDATA","Dimensions",{"x",size(PDATA,2),"y",size(PDATA,1)},'Datatype','double'); % create PDATA variable
        ncwrite(newfilename,"PDATA",PDATA'); % write PDATA
        
        nccreate(newfilename,"TWT","Dimensions",{"y",size(TWT,2)}); % TWT
        ncwrite(newfilename,"TWT",TWT); % positive half of TWT
        
        nccreate(newfilename,"UTMx","Dimensions",{"x",length(Rx)}); ncwrite(newfilename,"UTMx", Rx); % interpolated x
        nccreate(newfilename,"UTMy","Dimensions",{"x",length(Ry)}); ncwrite(newfilename,"UTMy", Ry); % interpolated y

        skycal_idx = rd.S.SkycalTraces; trace_idx = rd.S.ProfileTraces;
        if ~isempty(skycal_idx) % if skycal indices are available, create that variable in the netcdf file
            nccreate(newfilename,"skycal_idx","Dimensions",{"idx",length(skycal_idx)}); ncwrite(newfilename,"skycal_idx", skycal_idx);
        end

        % cleanup
        rd.D=[]; rd.TDATA=[];
        rd.PDATA=[];
    end
else % if smaller than batch size
    n=1; % only one batch
    rd.P.files=1:filenum; % grab file indices
    rd.data_dir=[datadisk datefolder 'FMCW/' transectfolder]; % raw data location
    rd=subdivide_daq(rd); % subdivide raw daq files
    % rd.TDATA=rd.TDATA-skycal; % perform sky calibration
    rd=cal_psd_radar(rd); % process for freq domain (calculate power spectral density)

    % vertical subset
    hmax_idx = 2^12; % subset index
    rd.PDATA=rd.PDATA(1:hmax_idx,:); % subset PDATA
    % plot(rd.TWT); hold on; plot(rd.TWT(8193:8193+hmax_idx-1));
    rd.TWT=rd.TWT(8193:8193+hmax_idx-1); % subset TWT from 2^13 fft

    % save to variable
    PDATA = rd.PDATA; CPUtime = (rd.CPUtime+7/24);
    % PDATA=[PDATA rd.PDATA]; % add next file
    % CPUtime=[CPUtime (rd.CPUtime+7/24)]; % convert to UTC (+6 for MDT / +7 for MST) and store

    % keep only traces without NaNs
    nan_cols_logical = any(~isnan(PDATA), 1); col_indices = find(nan_cols_logical);
    PDATA = PDATA(:, col_indices);
    CPUtime = CPUtime(:, col_indices);

    % plot result
    x=1:length(rd.CPUtime);
    depth=rd.TWT*velocity/2;
    % figure(1);clf;
    % imagesc(x,depth,PDATA,[min(PDATA,[],'all'), max(PDATA,[],'all')]); % may need to adjust cmax and cmin
    % xlabel('xidx','FontSize',fs); ylabel('Depth [cm]','FontSize',fs);
    % title(transectfolder(1:end-1),'Interpreter', 'none','FontSize',fs);
    % colorbar; ax = gca; ax.FontSize = fs;

    % now lets geolocate traces
    D=readtable(gpsdatafile);
    Ix=find(startsWith(string(D.Var1),'$GPGGA'));
    lat=zeros(length(Ix),1);
    lon=zeros(length(Ix),1);
    UTC=zeros(length(Ix),1);
    for m=1:length(Ix)
        S=D(Ix(m),1).Var1{1};
        D2=split(S,',');
        if length(D2{3}) > 1
            lat(m)=str2num(D2{3}(1:2))+(str2num(D2{3}(3:end)))/60;
            lon(m)=-(str2num(D2{5}(1:3))+(str2num(D2{5}(4:end))/60));
        else
            lat(m)=NaN;lon(m)=NaN;
        end
        if isempty(D2{2}) % if no timestamp
            UTC(m,1)=NaN;
        else
            UTC(m,1)=datenum(year,month,day,str2num(D2{2}(1:2)),str2num(D2{2}(3:4)),str2num(D2{2}(5:end))); % ENTER THE DATE
        end
    end
    % grab unique values only
    UTC = UTC(~isnan(UTC)); % remove empty timestamps
    UTC = UTC(isfinite(UTC)); % grab only finite values
    [UTCu, unique_idxs, i0] = unique(UTC);
    latu = lat(unique_idxs); lonu = lon(unique_idxs);
    % plot geolocated figure
    [x,y,zone]=ll2utm(latu,lonu); % convert to UTM
    Rx=interp1(UTCu,x,CPUtime); % get easting for CPUtimes
    Ry=interp1(UTCu,y,CPUtime); % get northing for CPUtimes
    % figure(2); clf; % figure params
    % ax=gca; ax.XAxis.FontSize = fs; ax.YAxis.FontSize = fs; axis equal;
    % plot(Rx,Ry,'.'); title('UTM coordinates for profile',FontSize=fs);
    % xlabel('Easting [m]',FontSize=fs); ylabel('Northing [m]','FontSize',fs);

    % write object to netcdf file
    filename = [transectfolder(1:end-1) '_' num2str(n) '.nc']; % grab the transect name, add batch number
    newfilename = [outputfilepath filename]; % create filepath and name
    if isfile(newfilename) % File exists, so delete it
        delete(newfilename);
    end

    TWT=rd.TWT; % grab PDATA & TWT
    PDATA=cast(PDATA,'double');
    % TWT = TWT(size(PDATA,1)+1:end); % grab positive half of TWT - not
    % necessary for GM

    nccreate(newfilename,"PDATA","Dimensions",{"x",size(PDATA,2),"y",size(PDATA,1)},'Datatype','double'); % create PDATA variable
    ncwrite(newfilename,"PDATA",PDATA'); % write PDATA

    nccreate(newfilename,"TWT","Dimensions",{"y",size(TWT,2)}); % TWT
    ncwrite(newfilename,"TWT",TWT); % positive half of TWT

    nccreate(newfilename,"UTMx","Dimensions",{"x",length(Rx)}); ncwrite(newfilename,"UTMx", Rx); % interpolated x
    nccreate(newfilename,"UTMy","Dimensions",{"x",length(Ry)}); ncwrite(newfilename,"UTMy", Ry); % interpolated y

    skycal_idx = rd.S.SkycalTraces; trace_idx = rd.S.ProfileTraces;
    if ~isempty(skycal_idx) % if skycal indices are available, create that variable in the netcdf file
        nccreate(newfilename,"skycal_idx","Dimensions",{"idx",length(skycal_idx)}); ncwrite(newfilename,"skycal_idx", skycal_idx);
    end

    % cleanup
    rd.D=[]; rd.TDATA=[];
    rd.PDATA=[];
end

% keep only traces without NaNs
nan_cols_logical = any(~isnan(PDATA), 1); col_indices = find(nan_cols_logical);
PDATA = PDATA(:, col_indices); 
CPUtime = CPUtime(:, col_indices);
disp([datefolder transectfolder 'complete.'])

%% plot result
x=1:length(rd.CPUtime);
depth=rd.TWT*velocity/2;
figure(4);clf;  
imagesc(x,depth,PDATA,[-110 -80]); % may need to adjust cmax and cmin
xlabel('xidx','FontSize',fs); ylabel('Depth [cm]','FontSize',fs);
title(transectfolder(1:end-1),'Interpreter', 'none','FontSize',fs);
colorbar; ax = gca; ax.FontSize = fs;

%% now lets geolocate traces - how does this connect to bigDriftWofNorthline data???
D=readtable([datadisk 'GM20200211/RADAR_GPS/gm11feb2020_geode.txt']);
Ix=find(startsWith(string(D.Var1),'$GPGGA'));
lat=zeros(length(Ix),1);
lon=zeros(length(Ix),1);
UTC=zeros(length(Ix),1);
for n=1:length(Ix)
    S=D(Ix(n),1).Var1{1};
    D2=split(S,',');
    if length(D2{3}) > 1
        lat(n)=str2num(D2{3}(1:2))+(str2num(D2{3}(3:end)))/60;
        lon(n)=-(str2num(D2{5}(1:3))+(str2num(D2{5}(4:end))/60));
    else
        lat(n)=NaN;lon(n)=NaN;
    end
    UTC(n,1)=datenum(year,month,day,str2num(D2{2}(1:2)),str2num(D2{2}(3:4)),str2num(D2{2}(5:end))); % ENTER THE DATE
end

% plot geolocated figure
[x,y,zone]=ll2utm(lat,lon); % convert to UTM
Rx=interp1(UTC,x,CPUtime); % get easting for CPUtimes
Ry=interp1(UTC,y,CPUtime); % get northing for CPUtimes
figure(1); clf % figure params
ax=gca; ax.XAxis.FontSize = fs; ax.YAxis.FontSize = fs; axis equal; 
plot(Rx,Ry,'.'); title('UTM coordinates for profile',FontSize=fs);
xlabel('Easting [m]',FontSize=fs); ylabel('Northing [m]','FontSize',fs);

%% write object to netcdf file
outputfilepath = '/Users/jukesliu/Documents/POSTDOC/snow-radar/GM2020/FMCW_processed/'; % path to output file
filename = [transectfolder(1:end-1) '.nc']; % grab the transect name
newfilename = [outputfilepath filename]; % create filepath and name
% newfilename = [filepath(1:end-4) '.nc']; % target file name

% File exists, so delete it
if isfile(newfilename)
    delete(newfilename);
end

TWT=rd.TWT; % grab PDATA & TWT
PDATA=cast(PDATA,'double');
% TWT = TWT(size(PDATA,1)+1:end); % grab positive half of TWT - not
% necessary for GM

nccreate(newfilename,"TWT","Dimensions",{"y",size(TWT,2)}); % TWT
ncwrite(newfilename,"TWT",TWT); % positive half of TWT

nccreate(newfilename,"UTMx","Dimensions",{"x",length(Rx)}); ncwrite(newfilename,"UTMx", Rx); % interpolated x
nccreate(newfilename,"UTMy","Dimensions",{"x",length(Ry)}); ncwrite(newfilename,"UTMy", Ry); % interpolated y

nccreate(newfilename,"PDATA","Dimensions",{"x",size(PDATA,2),"y",size(PDATA,1)},'Datatype','double'); % create PDATA variable
ncwrite(newfilename,"PDATA",PDATA'); % write PDATA

skycal_idx = rd.S.SkycalTraces; trace_idx = rd.S.ProfileTraces;

if ~isempty(skycal_idx) % if skycal indices are available, create that variable in the netcdf file
    nccreate(newfilename,"skycal_idx","Dimensions",{"idx",length(skycal_idx)}); ncwrite(newfilename,"skycal_idx", skycal_idx);
end