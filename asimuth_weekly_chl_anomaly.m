%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script:       asimuth_weekly_chl_anomaly
%
% Description:  This script creates  ASCII grd file of weekly MEAN chl anomaly
%               Uses the same "week number" logic as asimuth_buoy_sst_anomaly
%
% Inputs:       none
%  
% Outputs:      none
%               
% Author:       KL, July 2014
% Modified By:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dir_data=''; 

dir_out='';

%samba stuff
dir_samba='';
fsmbpwd='';

%The Ifremer chl files are usually delayed by a few days 
%The logic of the script below is that a new weekly anomaly file will be
%created a day or two AFTER the end of the week so what I check for here is
%the existence of the previous week's file and create it when all the data
%is present (usually will happen 1 or 2 days after the "week" ends)

%calculate week number of midnight
dnnow=now();
cwoy=week_number(dnnow);
if cwoy==0
   error('Week number is wrong'); 
end

%calculate the start and end dates of this week
%note I'm checking for existence of PREVIOUS week's file
dv=datevec(dnnow);
iyear=dv(1);
syear=num2str(iyear);
[dnstart,dnend]=week_dates(cwoy-1,iyear);

%does the anomaly netcdf file exist for the end of current week number
dvs=datevec(dnstart);
dve=datevec(dnend);
sdate = [num2str(dvs(1)) num2str(dvs(2),'%2.2d') num2str(dvs(3),'%2.2d')];
edate = [num2str(dve(1)) num2str(dve(2),'%2.2d') num2str(dve(3),'%2.2d')];

afile=[dir_data num2str(iyear) '' edate '.nc'];

wfile=[dir_out 'CA_' sdate '-' edate '.grd'];
wfileshort=['CA_' sdate '-' edate '.grd'];

%stop if this file does not exist or it and the weekly anomaly file exists
if (exist(afile,'file')~=2 || (exist(afile,'file')==2 && exist(wfile,'file')==2))
    exit;
end
    
%if we're here then try to create weekly anomaly
%get lon and lat data from an existing file
fname=[dir_data 'chl_anom_20130101.nc'];
lon=nc_var_read(fname,'longitude');
lat=nc_var_read(fname,'latitude');
lon=double(lon);
lat=double(lat);

ilon=length(lon);
ilat=length(lat);

%create an array for holding a week of data
wchl=zeros(ilon,ilat,7);

icount=0;
%now go through each day and find data if it exists
for k=dnstart:dnend
    cdate=datestr(k,'yyyymmdd');
    syear=datestr(k,'yyyy');
    fname=[dir_data syear '' cdate '.nc'];
    if exist(fname,'file')==2
        chl=nc_var_read(fname,'chl_anomaly');
        chl=double(chl);
        icount=icount+1;
        wchl(:,:,icount)=chl;
    else
        error(['CRITICAL File for ' cdate ' does not exist']);
    end
end

id=wchl<-100;
wchl(id)=nan;

%calculate mean nan using nanmean function
chl=nanmean(wchl,3);

%regrid AFTER averaging
%for grid to ascii need to have data in equal sized cells

nlat=[lat(1):0.015:lat(end)];

[X,Y]=meshgrid(lon,nlat);

glon=repmat(lon,1,ilat);
glat=repmat(lat',ilon,1);
Z=griddata(glon,glat,chl,X,Y);

id=find(isnan(Z));
Z(id)=-999;

%must go from North to South for SaveAsciiRaster function
flon=flipud(X); 
flat=flipud(Y);
fchl=flipud(Z);

sz=size(flon);
sz=sz(1)*sz(2);
fchl=reshape(fchl,sz,1);
flon=reshape(flon,sz,1);
flat=reshape(flat,sz,1);
out=[flon flat fchl];

%write data to raster
saveasciiraster_MI(out, wfile);

dir_win=['' syear ];
copy_to_samba(dir_samba,dir_win,dir_out,wfileshort,fsmbpwd);

quit;