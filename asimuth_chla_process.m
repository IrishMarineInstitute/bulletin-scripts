function asimuth_chla_process(dir_data,dir_samba,fsmbpwd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:     asimuth_chla_process
%
% Description:  Calculate 60-day running median from merged chlorophyll 
%               product produced by Ifremer and also calculate and plot
%               anomalies
%               Rolling 60-day median written to .mat file for latest date
%
% Inputs:       dir_data: directory where data resides
%               dir_samba: samba directory to copy to
%               fsmbpwd: smabe
% Outputs:      none
%               
% Author:       KL, July 2014
% Modified By:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fstub='-EUR-L4-CHL-ATL-v01-fv01-OI.nc';

dir_anom='';

%	ftp parameters
%
host='';
user='';
pass='';
dir_ftp='';

%file with location info (contents of this would need to change in Ifremer
%grid changes
grdname=[dir_data 'chl_grid.mat'];
load(grdname); %lon,lat,mask

%only use data west of 0 degrees
idlon=find(lon<0);
idlat=find(lat>=43 & lat<58); %smaller area than 2005 (start Nth Spain)

%convert to double for pcolor
lon=double(lon(idlon));
lat=double(lat(idlat));

datemidn=floor(now);

%operational median file - created initially by chl_medians function
%NB: the median file should not be under any year folder
medname=[dir_data 'chl_median.mat'];

%load rolling median file and get last date processed
load(medname);

%now go from last date processed to current date, searching for data

%There are 2 stages here
%       (i) Calculate 60 medians using new data
%       (ii) Calculate anomalies using medians 14 DAYS OLDER than current        

%KL: March 2015 - changing way I'm doing this to cope with file gaps
dnlast=dn_chl(end); %last good file we have
dstart=dn_chl(end)+1; %first day since last file
dnnew=dnlast;

%do we have new files BUT no gaps
blgap=0;
for k=dstart:datemidn
    sdate = datestr(k,'yyyymmdd');
    syear=datestr(k,'yyyy');
    fname=[dir_data syear '/NETCDF/' sdate fstub];
    if (exist(fname,'file')==2)
        %if there is gap then get out and email a warning
        if blgap
            %if data is missing the only thing to do is make a copy of
            %previous file and it will be used in median calc
            %missing files only seem to occur in winter so it's not a big
            %issue I think to use this kludge
            sMsg=['GAP EXISTS IN ASIMUTH CHL NETCDF DATA:: ' gapfname];
            email_message('ASIMUTH CHL DATA GAP',sMsg,'hpc.support@marine.ie');
            return; 
        else
            dnnew=k;
        end
    else
        blgap=1;
        gapfname=fname;
    end
end

%exit function if nonew file found
if dnnew==dnlast
   return; 
end

%Now we create the new median array (and time array) first
dn_chl=zeros(60,1);
chl_medians=zeros(60,800,1500);
    
icount=1;
for k=(dnnew-59):dnnew

    sdate = datestr(k,'yyyymmdd');
    syear=datestr(k,'yyyy');
    fname=[dir_data syear '/NETCDF/' sdate fstub];
    if (exist(fname,'file')==2)
        chl_time=nc_var_read(fname,'time');
        chl=nc_var_read(fname,'analysed_chl_a');
        sfac=ncreadatt(fname,'analysed_chl_a','scale_factor');
        fillval=ncreadatt(fname,'analysed_chl_a','_FillValue');
        mask=nc_var_read(fname,'mask');

        %subset
        chl=single(chl(idlon,idlat));
        mask=mask(idlon,idlat);

        %mask: sea=0,land=1,cloud=2
        %value of chl for land/cloud will be -32768 (prior to using scale_factor)
        %set chl values with mask of land/cloud to NaN
        id=find(mask>0 | chl==fillval); %also check for chl FillValue
        chl(id)=NaN;

        chl=chl*sfac; %apply scale factor

        chl_medians(icount,:,:)=chl;
        dn_chl(icount)=k;
    else
        %if the file does not exist (i.e. there is a gap of a day or 2)
        %use the previous day's chl data
        %WARNING: OBVIOUSLY A LARGE GAP WILL SKEW STATS
        %ths chl array is from previous file but "k" has changed
        chl_medians(icount,:,:)=chl;
        dn_chl(icount)=k;
        disp([fname ' missing: using previous available file for data'])
    end

    icount=icount+1;
end


for k=dstart:datemidn
    sdate = datestr(k,'yyyymmdd');
    syear=datestr(k,'yyyy');
    fname=[dir_data syear '/NETCDF/' sdate fstub];
    
    %NB: does not expect gaps in data so first date that's missing is assumed
    %to be latest date
    
    if (exist(fname,'file')==2)
        chl_time=nc_var_read(fname,'time');
        chl=nc_var_read(fname,'analysed_chl_a');
        sfac=ncreadatt(fname,'analysed_chl_a','scale_factor');
        fillval=ncreadatt(fname,'analysed_chl_a','_FillValue');
        mask=nc_var_read(fname,'mask');
        
        %subset
        chl=single(chl(idlon,idlat));
        mask=mask(idlon,idlat);

        %mask: sea=0,land=1,cloud=2
        %value of chl for land/cloud will be -32768 (prior to using scale_factor)
        %set chl values with mask of land/cloud to NaN
        id=find(mask>0 | chl==fillval); %also check for chl FillValue
        chl(id)=NaN;

        chl=chl*sfac; %apply scale factor

        %CALCULATE MEDIAN from time array
        %If an array has NaN values in it median will return NaN!!!
        %This code is going to be slow but I can't think of another way to do it
        %Minimize work by avoiding land altogether - bizarrely the land
        %mask only sees Spain and southern France as land!!!!
        %so use max function to isolate time series that are not composed
        %only of NaNs
        mx=squeeze(max(chl_medians));
        idsea=find(~isnan(mx));
        %create an array of -1's
        chl_median60=ones(length(idlon),length(idlat))*-1.0;
        for isea=1:length(idsea)
            chl_ts=squeeze(chl_medians(:,idsea(isea)));
            chl_median60(idsea(isea))=median(chl_ts(~isnan(chl_ts)));
        end
        
        disp(['Created median for ' sdate]);
        %put median file in correct year folder
        dir_output=[dir_data syear '/NETCDF/Median/'];
        if ~(exist(dir_output,'dir')==7)
            mkdir(dir_output);
        end        
        
        mfile=[dir_output sdate '.mat'];
        save(mfile,'chl_median60');
        clear chl_median60; %prevent any chance of using wrong median

        %KLC
        %to cope with looping through more than 1 day
        %move time series forward in 60-day array and add in latest data for next median calc
        chl_medians(1:59,:,:)=chl_medians(2:60,:,:);
        chl_medians(60,:,:)=chl;
        
        
        %(ii) Calculate anomaly using median from 14 DAYS AGO
        sdate2=datestr(k-14,'yyyymmdd');
        syear2=datestr(k-14,'yyyy');
        dir_output=[dir_data syear2 '/NETCDF/Median/'];
        fname=[dir_output sdate2 '.mat'];
        disp(['Using median from ' sdate2 ' for anomaly']);
        
        load(fname);
        
        %CREATE pic and save anomaly + median in netcdf file
        
        %values of -1 in median_chl_60 are land/nodata
        id0=find(chl_median60==-1);
    
        chl_anomaly=chl-chl_median60;
    
        %set median land/nodata to 0
        chl_anomaly(id0)=0;
        chl(id0)=0;
        
        %convert to double for pcolor
        chl_anomaly=double(chl_anomaly);
        chl=double(chl);
        
        %just look at ireland for now and use a tight clims to highlight
        %anomalies
        idlonI=find(lon<-4);
        idlatI=find(lat>48);       
        
        %make plot
        %LOG SCALE
        chl_a=chl_anomaly(idlonI,idlatI);
        chll=chl(idlonI,idlatI);
        %need to get rid of 0 and -ve for log scale
        chl_a(chl_a<=0)=0.0000001; 
        chll(chll<=0)=0.0000001;
        chl_a=log10(chl_a);
        chll=log10(chll);
        clims=[log10(0.1) log10(20)];
        dir_pic=[dir_data syear '/IMAGES/Anomaly/'];
        if ~(exist(dir_pic,'dir')==7)
            mkdir(dir_pic);
        end            
        ifile_name=['chl_anom_' sdate '.png']; 
        iname=[dir_pic ifile_name];
        cb_yticks=[log10(0.1) log10(1) log10(2) log10(5) log10(10) log10(20)];
        cb_ylabels=[0.1 1 2 5 10 20];

        strdate=datestr(k);

        chl_anomaly_map(lon(idlonI),lat(idlatI),chll,chl_a,...
            clims,cb_yticks,cb_ylabels,strdate,iname);

        dir_win=['' syear];
        copy_to_samba(dir_samba,dir_win,dir_pic,ifile_name,fsmbpwd);
        
        %put image on FTP site for others to access
        f=ftp(host,user,pass);
        cd(f,dir_ftp);
        %put new file on site
        a=mput(f,iname);
        close(f);
        
        %write data to netcdf file
        chl_a=chl_anomaly(idlonI,idlatI); %includes +ve and -ve
        chl_a(isnan(chl_a))=-9999;
        chl_a(chl_a==0)=-9999;
        dir_output=[dir_data syear '/NETCDF/Anomaly/'];
        if ~(exist(dir_output,'dir')==7)
            mkdir(dir_output);
        end   
        ncname=[dir_output 'chl_anom_' sdate '.nc'];
        create_chla_anomaly_nc(ncname,chl_time,lon(idlonI),lat(idlatI),chl_a);
        
        %write chl data to grd file for surfer
        lonS=double(lon(idlonI));
        latS=double(lat(idlatI));
        chl=chl(idlonI,idlatI);
        ilon=length(lonS);
        ilat=length(latS);
        chl(chl==0)=nan;
        %for grid to ascii need to have data in equal sized cells
        nlat=[lat(1):0.015:lat(end)];

        [X,Y]=meshgrid(lonS,nlat);

        glon=repmat(lonS,1,ilat);
        glat=repmat(latS',ilon,1);
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
        wfile=[dir_anom '/chl_' sdate '.grd']; 
        saveasciiraster_MI(out, wfile);
        dir_win=['' syear];
        gname=['chl_' sdate '.grd'];
        copy_to_samba(dir_samba,dir_win,dir_anom,gname,fsmbpwd);
        
        
    end
end

%after processing save the time of the latest medians file for use next time
%NOTE this is only saved if no errors above 
%IF there are errors then the script will start next time from the last processed day
%recorded by the last error-free running of this script
save(medname,'dn_chl');

return