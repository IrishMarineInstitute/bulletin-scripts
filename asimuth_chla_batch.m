%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script:       asimuth_chla_batch
%
% Description:  This script downloads the merged chlorophyll-a product
%               produced by Ifremer and calls asimuth_chla_process to 
%               calculate latest 60-day medians and anomalies
%
% Inputs:       none
%  
% Outputs:      none
%               
% Author:       KL, July 2014
% Modified By:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dir_data = '';

ftp_host = '';
dir_ftp = '';

%Merged high resolution multi-sensor chlorophyll-a data processed by
%Ifremer, France.


dir_samba='';
fsmbpwd='';

%look at last 20 days just to ensure picking up any stragglers
dend=now;
dstart=dend-20;

%1: download any new data available

%this copes with going between years
for k=dstart:dend
    [iyear,imon,iday,ihr,imin,isec] = datevec(k);
    syear=num2str(iyear);
    diryear=[dir_data syear];
    dir_local=[diryear '/NETCDF'];
    dirimage=[diryear '/IMAGES'];
    if ~(exist(diryear,'dir')==7)
        mkdir(diryear);
        mkdir(dir_local);
        mkdir(dirimage);
    end
    
    %calculate day number
    idy=datenum([iyear imon iday 0 0 0])-datenum([iyear 1 1 0 0 0])+1;
    tdir_ftp=[dir_ftp num2str(iyear) '/' num2str(idy,'%3.3d') '/'];
    sdate = datestr(k,'yyyymmdd');
    nfile_name=[sdate '-EUR-L4-CHL-ATL-v01-fv01-OI.nc'];
    if (exist([dir_local '/' nfile_name],'file')~=2) 
        %download file
        ftp_file = [ftp_host tdir_ftp nfile_name '.bz2'];
        strwget=['wget --no-proxy -t 200 --waitretry=100 --retry-connrefused -c -P' dir_local ' ' ftp_file];
        system(strwget);
        %unzip file
        system (['bunzip2 ' dir_local '/' nfile_name '.bz2']);
        disp(['Downloaded ' nfile_name]);
    end  
    ifile_name=['analysed_chl_' sdate '.png'];
    if (exist([dirimage '/' ifile_name],'file')~=2) 
        %download file
        ftp_file = [ftp_host tdir_ftp ifile_name];
        strwget=['wget --no-proxy -t 200 --waitretry=100 --retry-connrefused -c -P' dirimage ' ' ftp_file];
        system(strwget);
        %put images on OceanSQL for Caroline to access
        dir_win=['' num2str(iyear)];
        copy_to_samba(dir_samba,dir_win,dirimage,ifile_name,fsmbpwd);
    end 
end

disp('Calculating medians and anomalies');

%now call function to create medians and anomalies
asimuth_chla_process(dir_data,dir_samba,fsmbpwd);

clear

quit;