function asimuth_bantry_particle_stats(ffname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:     asimuth_bantry_particle_stats
%
% Description:  Create a netcdf file containing float particle stats from
%               operational forecast releases in Bantry model domain. 
%
% Inputs:       fname: float output file name
%
% Outputs:      none
%               
% Author:       KL, Aug 2014
% Modified By:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%constants
istart = 1;
istride = 350; %
istop = istride;
idaylimit = 3*24*3; %3days * 24hr  * 3 (20min intervals)
nrel=6; %number of release "lines"
outdir='';

%operational grid
gname='bantry_bay.nc';
lon_rho=double(nc_var_read(gname,'lon_rho'));
lat_rho=double(nc_var_read(gname,'lat_rho'));

sz1=size(lon_rho);

%older cartesian grid which will interpolate to for ERDDAP
gname2='bantry_bay_v2.nc';
lon=double(nc_var_read(gname2,'lon_rho'));
lat=double(nc_var_read(gname2,'lat_rho'));

sz2=size(lon);

%1: Read  float  file data
xgrid = nc_var_read(ffname,'Xgrid');
ygrid = nc_var_read(ffname,'Ygrid');
ot = nc_var_read(ffname,'ocean_time');

%2: create Netcdf file
dn=(ot(1)/(24*3600))+datenum(1968,5,23,0,0,0);
dstub=datestr(dn,'yyyymmddHH');
fstub=['BP_' dstub '.nc'];
nc_name = [outdir fstub];
create_bantry_particles_nc(nc_name,sz2(1),sz2(2),nrel);

%3: Process data
%Initialize the array for particle data
particle_hours = zeros([sz2 6]);

%create 1D vectors of lat/lon for TriScatteredInterp
xx=reshape(lon_rho,sz1(1)*sz1(2),1);
yy=reshape(lat_rho,sz1(1)*sz1(2),1);

for m=1:6
    x1 = reshape(xgrid(istart:istop,1:idaylimit),[],1);
    y1 = reshape(ygrid(istart:istop,1:idaylimit),[],1);
    
    ind = find(x1<600 & x1>0); %discard bad data

    x1=max(floor(x1(ind)),1);
    y1=max(floor(y1(ind)),1);

    %aggregate counts in each cell
    pcount = zeros(sz1);
    for n=1:size(ind)
        pcount(x1(n),y1(n)) = pcount(x1(n),y1(n))+1;
    end

    ph=pcount/3;

    %Interpolate onto regular grid
    %resize count matrix to be 1D vector
    ph=reshape(ph,sz1(1)*sz1(2),1);
    TSI = TriScatteredInterp(xx,yy,ph);
    Z = TSI(lon,lat);
    
    particle_hours(:,:,m)=Z;
    
    clear pcount ind x1 y1 ph pho;
           
    istart = istart + istride;
    if(m == 3) istride = 275; end; %change stride value for Bantry releases
    istop = istop + istride;
    
end     

%check for nans and zeros and replace with missing value
id=find(particle_hours==0 | isnan(particle_hours));
particle_hours(id)=-999;

%3: write data to Netcdf file
%write variables to file

%change ocean_time to secs from 1/1/2000 - only use first ot value
ot=24*3600*((ot(1)/(24*3600)+datenum(1968,5,23,0,0,0))-datenum(2000,1,1,0,0,0));
nc_var_write(nc_name,'time',ot);
nc_var_write(nc_name,'lon',lon);
nc_var_write(nc_name,'lat',lat);
nc_var_write(nc_name,'particle_hours',particle_hours);

%4: Publish Netcdf file to Thredds
thredds_dir='';

%delete files other than 30 days
delete_from_thredds(thredds_dir,30);

%copy new file over
thredds_file=[thredds_dir fstub];
copy_to_thredds(nc_name,thredds_file);

return