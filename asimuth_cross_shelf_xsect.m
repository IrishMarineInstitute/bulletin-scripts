function asimuth_cross_shelf_xsect(aname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script:     asimuth_cross_shelf_xsect
%
% Description:  This function plots a cross section for a transect from
%               shelf to inner Bantry. Section created for T,S & RHO
%
% Inputs:       aname: average file name
%
% Outputs:      none
%               
% Author:       KL, Aug 2015
% Modified By:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirplots='';

%fixed scales and contours
ncol=120;
tclim=[8 20];
tcontour=[8:1:20];
sclim=[34.4 35.6];
scontour=[34.4:0.1:35.6];
rclim=[25 27.4];
rcontour=[25:0.2:27.4];

ot=nc_var_read(aname,'ocean_time');
lonr=nc_var_read(aname,'lon_rho');
latr=nc_var_read(aname,'lat_rho');
h=nc_var_read(aname,'h');


%Below is for defining string to uniquely name the figure (need to subtract
%86400 from ot(1) as the averaging time period is 67.5hrs and thus the
%ot(1) being in the middle is equal 07:10hrs on the next day
dn1=floor(roms_time_to_datenum(ot(1)-86400));

datestub=datestr(dn1,'yyyymmdd');

%because this is a high res grid can use start and end L & M positions for 
%section line. This makes it much easier/quicker to process data
iL0=71;
iM0=124;
iL1=189;
iM1=384;
nparts=200;
cL=zeros(nparts,1);
cM=zeros(nparts,1);
for k=1:nparts
    t=k/nparts;
    cL(k)=iL0*(1-t)+iL1*t;
    cM(k)=iM0*(1-t)+iM1*t;
end
%round the values to make them integers
cL=round(cL);
cM=round(cM);


theta_s=4.0;
theta_b=0.85;
hc=3.0;
N=20;
igrid=1; 

zeta=nc_var_read(aname,'zeta');
t=nc_var_read(aname,'temp');
s=nc_var_read(aname,'salt');
rho=rho_pot(t,s);

lons=zeros(nparts,1);
lats=zeros(nparts,1);
zs=zeros(nparts,N);
ts=zeros(nparts,N);
ss=zeros(nparts,N);
rs=zeros(nparts,N);

for k=1:nparts
    iL=cL(k);
    iM=cM(k);
    lons(k)=lonr(iL,iM);
    lats(k)=latr(iL,iM);
    %loop through points and create a profile for each
    z=squeeze(set_depth(1,1,theta_s,theta_b,hc,N,igrid,h(iL,iM),zeta(iL,iM),0));
    zs(k,:)=z;
    ts(k,:)=t(iL,iM,:);
    ss(k,:)=s(iL,iM,:);
    rs(k,:)=rho(iL,iM,:);
end

%this gives longitude to distance conversion for this latitude
%need to provide distances so that the interpolation works properly
dist=spheric_dist(lats(1),lats(1),lons(1),lons(2));
fac=dist/abs(lons(1)-lons(2));

%create the data matrices
lons2D=repmat(lons,1,20)*fac;
id=find(ts>-273);
%temperature
TSI = TriScatteredInterp(lons2D(id),zs(id),ts(id));
x=lons*fac;
y=floor(min(min(zs))):0.5:0;
[qx,qy]=meshgrid(x,y);
Tx=TSI(qx,qy);
%salt
TSI = TriScatteredInterp(lons2D(id),zs(id),ss(id));
Sx=TSI(qx,qy);
%density
TSI = TriScatteredInterp(lons2D(id),zs(id),rs(id));
Rx=TSI(qx,qy)-1000.0;

%this is kind of a crude way of blanking out data that extends below bottom
%depth at each location - but it works!
for k=1:nparts
    z=zs(k,:);
    qz=qy(:,k);
    id=find(qz<min(z));
    Tx(id,k)=nan;
    Sx(id,k)=nan;
    Rx(id,k)=nan;
end


%Create Plots
%convert qx back to longitude for plotting
qx=qx./fac;

%temperature
strTitle = ['Cross Shelf Section Temperature (Forecast Date: ' datestr(dn1,'dd-mmm-yyyy') ')'];
clim=tclim;

hf = figure('Color','w','visible','off','Position',[100 100 1000 500]);
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');

%use colorbrewer for colormap creation
cmap=cbrewer('div', 'RdYlBu', ncol);

%Want to go from blue to red
%to reverse colours use colormap(flipud(cmap))
colormap(flipud(cmap));

pcolor(qx,qy,Tx)
shading flat
caxis(clim);
hold on
[C,h1]=contour(qx,qy,Tx,tcontour,'k');
%the following makes seabed look better
plot(lons,squeeze(zs(:,1)),'w','Linewidth',7)
th=clabel(C,h1,'LabelSpacing',1000,'Rotation',0,'FontSize',12);
for c=1:length(th)
    cl=get(th(c),'UserData');
    cl=char(num2str(cl,'%5.1f'));
    set(th(c),'String', cl);
end

set(gca,'ytick',-140:20:0);
xlabel('Longitude')
ylabel('Depth (m)')
title(strTitle,'FontSize',14)
colorbar('YTick',tcontour,'YTickLabel',sprintf('%5.1f|',tcontour));
box on
fname=[dirplots 'Bantry_Shelf_Section_T_' datestub '.png'];
print ('-dpng', fname);
close(hf);


%salinity
strTitle = ['Cross Shelf Section Salinity (Forecast Date: ' datestr(dn1,'dd-mmm-yyyy') ')'];
clim=sclim;

hf = figure('Color','w','visible','off','Position',[100 100 1000 500]);
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');

%use colorbrewer for colormap creation
cmap=cbrewer('div', 'RdBu', ncol);
colormap(flipud(cmap));

pcolor(qx,qy,Sx)
shading flat
caxis(clim);
hold on
[C,h1]=contour(qx,qy,Sx,scontour,'k');
plot(lons,squeeze(zs(:,1)),'w','Linewidth',7)
th=clabel(C,h1,'LabelSpacing',1000,'Rotation',0,'FontSize',12);
for c=1:length(th)
    cl=get(th(c),'UserData');
    cl=char(num2str(cl,'%5.1f'));
    set(th(c),'String', cl);
end
set(gca,'ytick',-140:20:0);
xlabel('Longitude')
ylabel('Depth (m)')
title(strTitle,'FontSize',14)
colorbar('YTick',scontour,'YTickLabel',sprintf('%5.1f|',scontour));
box on
fname=[dirplots 'Bantry_Shelf_Section_S_' datestub '.png'];
print ('-dpng', fname);
close(hf);


%density
strTitle = ['Cross Shelf Section Density (Forecast Date: ' datestr(dn1,'dd-mmm-yyyy') ')'];
clim=rclim;

hf = figure('Color','w','visible','off','Position',[100 100 1000 500]);
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');

%use colorbrewer for colormap creation
cmap=cbrewer('div', 'RdBu', ncol);

colormap(flipud(cmap));
pcolor(qx,qy,Rx)
shading flat
caxis(clim);
hold on
[C,h1]=contour(qx,qy,Rx,rcontour,'k');
plot(lons,squeeze(zs(:,1)),'w','Linewidth',7)
th=clabel(C,h1,'LabelSpacing',1000,'Rotation',0,'FontSize',12);
for c=1:length(th)
    cl=get(th(c),'UserData');
    cl=char(num2str(cl,'%5.1f'));
    set(th(c),'String', cl);
end
set(gca,'ytick',-140:20:0);
xlabel('Longitude')
ylabel('Depth (m)')
title(strTitle,'FontSize',14)
colorbar('YTick',rcontour,'YTickLabel',sprintf('%5.1f|',rcontour));
box on
fname=[dirplots 'Bantry_Shelf_Section_R_' datestub '.png'];
print ('-dpng', fname);
close(hf);

return;