function asimuth_plot_2D_fluxes_plan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:     asimuth_plot_2D_fluxes_plan
%
% Description:  Create plan view plots of Bantry Bay volumetric fluxes
%               through pre-defined section lines
%
% Inputs:       none
%
% Outputs:      none
%               
% Author:       TD, 2013
% Modified By:  KL, Aug 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%base location for writing plots
dirplots='';

%Input arguments to ultimate function
filenm = 'bantry_bay_avg_0001.nc';
hzslice = -20;   %depth at which w velocities to be plotted in plan view.

%Color ranges for contouring various images.
Wvel_cmax = 10; %mm/s
Flux_cmax = 500; %m3/s
BantryFlux_cmax = 100; %m3/s

% cross-sections requried (hardwired).
off = [1 1 1 199; 72 72 1 199; 118 118 1 199; 170 170 1 199; 207 207 1 199; 283 283 1 199; 1 1 201 319; 72 72 201 230; 118 118 201 292; 170 170 201 292; 207 207 201 267; 283 283 201 313];
along = [1 72 200 200; 72 118 200 200; 118 170 200 200; 170 207 200 200; 207 283 200 200; 118 169 290 290; 160 190 360 360];

%Get dimensions of model domain from netCDF file.
srho_len=nc_dimlen(filenm,'s_rho');
etarho_len=nc_dimlen(filenm,'eta_rho');

%get general data from netCDF file
%note: this script was origonally written using Mex netcdf functions which
%returns matrices cols and ros opposite to Matlab native function. Rather
%than rewriting all the code later on, we transpose the arrays here to have
%in same orientation as previous
masku=nc_var_read(filenm,'mask_u')';
huon=nc_var_read(filenm,'Huon'); %u-flux
huon=permute(huon,[3 2 1]);
h=nc_var_read(filenm,'h')';
srho=nc_var_read(filenm,'s_rho')';
ot=nc_var_read(filenm,'ocean_time');

%Below is for defining string to uniquely name the figure (need to subtract
%86400 from ot as the averaging time period is 67.5hrs and thus ot falls on
%07:10hrs of the 1st day following the start of the forecast
dt1=roms_time_to_datenum(ot(1)-86400);
dt2 = dt1+(12/24);
dt3 = dt1+(72/24);
dtstr(:,1) = datevec(dt1);
dtstr(:,2) = datevec(dt2);
dtstr(:,3) = datevec(dt3);

%Open a matfile for archiving data
datestub=[num2str(dtstr(3,3),'%2.2u') '_' num2str(dtstr(2,3),'%2.2u') '_' num2str(dtstr(1,3),'%4.4u')];

mfile=['' datestub '.mat'];
save(mfile,'ot','huon');

%create figure 1 and figure 2 for plotting into
fig1 = figure('Name','Bantry Flux Sections 1','NumberTitle','off','Units','pixels','Position',[100 100 900 900],'Color',[0.9 0.9 0.9],'visible','off');
fig2 = figure('Name','Bantry Flux Plans','NumberTitle','off','Units','pixels','Position',[100 100 900 900],'Color',[0.9 0.9 0.9],'visible','off');
fig3 = figure('Name','Bantry Flux Sections 2','NumberTitle','off','Units','pixels','Position',[100 100 1200 300],'Color',[0.9 0.9 0.9],'visible','off');
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');

%load default colormap (red -> blue)
colormap('default');
load('FluxColorMap','flux_cmap');

%for each eta plane cross section  
for m = 1:size(off,1);
    
    %extract bathymetry of x-section from [h]
	h_1= squeeze(h(off(m,3):off(m,4),off(m,1)))'; 
    %extract u-velocity mask of x-section from [masku]
    mask = squeeze(masku(off(m,3):off(m,4),off(m,1)))';
    
    %for each sigma level
    for n=1:srho_len
        %extract masked u-flux through x-sect from [huon]
    	flux(n,:) = squeeze(huon(n,off(m,3):off(m,4),off(m,1))).*mask;
        %create array of depths at each grid point in x-sect
        z(n,:) = h_1.*srho(1,n);
        %create array of eta points at each grid point in x-sect
        x(n,:) = off(m,3):off(m,4);
    end   
        
    %define logical array for all u-fluxs >0
    pos = flux>0;
    %calculate sum of positive and negative fluxes and net flux.
    sum_pos = sum(sum(flux(pos==1)));
    sum_neg = sum(sum(flux(pos==0)));
    sum_net = sum_pos-sum_neg;  
     
    %x / y location along x-sect for plotting u-flux arrows
    loc_x = round((off(m,2)-off(m,1))/2)+off(m,1);            %x
    loc_y = round((off(m,4)-off(m,3))/2)+off(m,3);            %y

    %create title text for each x-sect plot
    string = ['Section A' num2str(m)];

   
    figure(fig1);   %Switch to fig1 plot

    if(m < 7)
        subplot(6,2,(m*2)-1);
        left = 0.05;
        bottom = 1.0-(m*0.15);
        width = 0.4; 
        height = 0.075;
    else
        subplot(6,2,(m-6)*2);
        left = 0.6;
        bottom = 1-((m-6)*0.15);
        width = 0.4*((max(max(x))-min(min(x)))/200);
        height = 0.075;
    end;
        
    %plot the x-sect
    pcolor(x,z,double(flux));shading interp;hold on;
    %add the title to the x-sect
    text(off(m,3),z(srho_len,1)+30,string,'rotation',0,'fontsize',10,'fontweight','bold');
    %plot bathymetry profile along x-sect
    line(x(1,:),z(1,:),'Color','k','LineWidth',1);hold on;
    %plot Hzslice line on x-sect
    line([off(m,3) off(m,4)] ,[hzslice hzslice],'LineStyle','--','Color','k','LineWidth',1);hold on;
    %constrain the x-sect axes
    axis([min(min(x)) max(max(x)) -max(max(h)) 0]);     
    %constrain the x-sect color ramping
    caxis([-Flux_cmax Flux_cmax]);
    %apply the colormap to the x-sect u-fluxes
    set(gcf,'Colormap',flux_cmap);
    %apply the x-sects plot coordinates 
    set(gca,'position',[left  bottom  width height]);
    %add the color bar to the x-sect plot
    colorbar('position',[left+width+0.01  bottom  0.01 height],'Ytick',[-Flux_cmax 0 Flux_cmax]);

    figure(fig2)    %Switch to fig 2 plot
    
    %add the line denoting the x-sect
    line([off(m,1),off(m,2)],[off(m,3),off(m,4)],'Color','k','LineWidth',4);hold on;     
    %add the x-sect title text
    text(off(m,1)-5,off(m,3)+5,string,'rotation',90,'fontsize',10,'fontangle','italic');
    %add the positive and negative flux arrows
    quiver(loc_x,loc_y,sum_pos/10000,0,'r','LineWidth',4,'maxheadsize',1);hold on;
    quiver(loc_x,loc_y,sum_neg/10000,0,'b','LineWidth',4,'maxheadsize',1);hold on;

    %clear out variables for use in next x-sect
    clear mask flux z x;
end
%All eta-section x-sects now plotted

%clear arrays for use in xi-section
clear mask_u Huon m;

% get v data from netCDF file
maskv=nc_var_read(filenm,'mask_v')';
hv=nc_var_read(filenm,'Hvom'); %v-flux
hv=permute(hv,[3 2 1]);

% Save to a mat file
save(mfile,'hv','-append');

%rearrange hvom array indices
hvom = permute(hv,[1 3 2]);

%for each xi plane cross section
for m = 1:size(along,1);
    
    %extract bathymetry of x-section from [h]
    h_1 = h(along(m,3),along(m,1):along(m,2))';
    %extract v-velocity mask of x-section from [maskv]
    mask = maskv(along(m,3),along(m,1):along(m,2));

    %for each sigma level
    for n=1:srho_len      
        %extract masked v-flux through x-sect from [hvom]
        flux(n,:) = hvom(n,along(m,1):along(m,2),along(m,3)).*mask;      
        %create array of depths at each grid point in x-sect
        z(n,:) = h_1.*srho(1,n);
        %create array of xi points at each grid point in x-sect
        x(n,:) = along(m,1):along(m,2);
    end
    
    %define logical array for all v-fluxs >0
    pos = flux>0;
    %calculate sum of positive and negative fluxes and net flux.
    sum_pos = sum(sum(flux(pos==1)));
    sum_neg = sum(sum(flux(pos==0)));
    sum_net = sum_pos-sum_neg;
    
    %x / y location along x-sect for plotting v-flux arrows
    loc_x = round(((along(m,2)-along(m,1))/2)+along(m,1));      %x
    loc_y = round((along(m,4)-along(m,3))/2)+along(m,3);        %y

    %create title text for each x-sect plot
    string = ['Section B' num2str(m)];       

    figure(fig3);   %Switch to fig1 plot
    %for all x-sects (except Bantry mouth) set up plot coordinates
    if(m < 6)  
        subplot(2,5,m+5);   
        width = 0.3*((max(max(x))-min(min(x)))/200);
        left = 0.1*(m*1.75)-0.1;
        bottom = 0.15;     
        height = 0.2;
    else
        %set up plot coordinates for Bantry mouth x-sect
        if(m == 6)
            subplot(2,5,1:3);
            width = 1.75*((max(max(x))-min(min(x)))/200);
            left = 0.075;
        elseif(m==7)
            subplot(2,5,4:5);
            width = 1.75*((max(max(x))-min(min(x)))/200);
            left = 0.6;
        end
        bottom = 0.5;     
        height = 0.4;
    end
        
    %plot the x-sect
    pcolor(x,z,double(flux));shading interp;hold on;       
    %plot bathymetry profile along x-sect
    line(x(1,:),z(1,:),'Color','k','LineWidth',1);hold on;           
    %plot Hzslice line on x-sect
    line([along(m,1) along(m,2)] ,[hzslice hzslice],'LineStyle','--','Color','k','LineWidth',1);hold on;

    if(m<6)
        %constrain the x-sect axes
        axis([min(min(x)) max(max(x)) -max(max(h)) 0]);
        %add the title to the x-sect
        text(along(m,1),z(srho_len,1)+30,string,'rotation',0,'fontsize',10,'fontweight','bold');
        %constrain the x-sect color ramping
        caxis([-Flux_cmax Flux_cmax]);
        %add the color bar to the x-sect plot
        colorbar('position',[left+width+0.01  bottom  0.01 height],'Ytick',[-Flux_cmax 0 Flux_cmax]);              
   else
        %constrain the Bantry mouth x-sect axes
        axis([min(min(x)) max(max(x)) -65 0]);

        %add the title to the Bantry mouth x-sect
        if(m==6)text(along(m,1),z(srho_len,1)+5,[string ': Mouth of Bantry Bay'],'rotation',0,'fontsize',10,'fontweight','bold');end
        if(m==7)text(along(m,1),z(srho_len,1)+5,[string ': Inner Bantry Bay'],'rotation',0,'fontsize',10,'fontweight','bold');end

        %constrain the x-sect color ramping
        caxis([-BantryFlux_cmax BantryFlux_cmax]);
        colorbar('position',[left+width+0.01  bottom  0.01 height],'Ytick',[-BantryFlux_cmax 0 BantryFlux_cmax]);              
    end
        
    %apply the colormap to the x-sect v-fluxes
    set(gcf,'Colormap',flux_cmap);
    %apply the x-sects plot coordinates 
    set(gca,'position',[left  bottom  width height]);          
            
    figure(fig2);   %Switch to fig 2 plot
                      
    %add the line denoting the x-sect
    line([along(m,1),along(m,2)],[along(m,3),along(m,4)],'Color','k','LineWidth',4);hold on; 
    %add the x-sect title text
    text(along(m,1)+5,along(m,3)+5,20,string,'rotation',0,'fontsize',10,'fontangle','italic');
    %add the positive and negative flux arrows
    quiver(loc_x,loc_y,0,sum_pos/10000,'r','LineWidth',4,'maxheadsize',1);hold on;
    quiver(loc_x,loc_y,0,sum_neg/10000,'b','LineWidth',4,'maxheadsize',1);hold on;

    %clear out variables        
    clear mask flux z x;
end
%All xi-section x-sects now plotted

%clear arrays 
clear mask_v hv;

%Move to plotting the horizontal colored contour of w velocity
%get rho mask
mask_rho=nc_var_read(filenm,'mask_rho')';
%mask bathymetry
h = h.*mask_rho;
%create 3d maskrho array and rearrange indices to 'normal' ROMS indices
maskrho = permute(repmat(mask_rho,[1 1 21]),[3 1 2]);

% get general data from netCDF file.
zeta=nc_var_read(filenm,'zeta')';
theta_s=nc_var_read(filenm,'theta_s');
theta_b=nc_var_read(filenm,'theta_b');
hc=nc_var_read(filenm,'hc');
w=nc_var_read(filenm,'w');
w=permute(w,[3,2,1]);

%save to a mat file
save(mfile,'zeta','w','-append');

%calculate 3d arrays of z-levels at each sigma level
z = zlevs(h,zeta,theta_s,theta_b,hc,srho_len,'w',1);

%convert w to mm/s
w = w.*maskrho*1000;        %convert from m/s to mm/s

%interpolate to find wvel at chosen Hzslice depth
wvel = vinterp(w,z,hzslice);

figure(fig2);   %Switch to fig2 plot
   
%Flux Legend

%create coordinate arrays for plotting legend
a(1:5) = 40;
b(1:5) = [340:((380-340)/4):380];
c(1:5) = [(Flux_cmax/5):(Flux_cmax-(Flux_cmax/5))/4:Flux_cmax];
d(1:5) = 0;
    
%plot legend arrows
quiver(a,b,c,d,'k','LineWidth',2,'maxheadsize',1);hold on;

%plot legend text
for n = 1:5;
    text(a(n)-40,b(n),[num2str(c(n),3) ' m^3/s'],'rotation',0,'fontsize',10,'fontweight','bold');
end;

%plot the 2D wvelocity color contour at the chosen Hzslice depth 
pcolor(double(wvel));shading interp;

%plot the land/sea contour
contour(h,[-1:0],'k');hold on;

%%add the title to the contour
text(0,etarho_len-30,'Volumetric Fluxes','rotation',0,'fontsize',10,'fontweight','bold');

%Add colorbar & title
colbar = colorbar('Ytick',[-Flux_cmax:1:Flux_cmax]);
title(colbar, {'Vertical velocity',' at 20m (mm/sec)'});

%Axis manipulation
axis([-25 315 1 etarho_len]);

%constrain the contour color ramping
caxis([-Wvel_cmax Wvel_cmax]);set(gca,'DataAspectRatio',[1 1 1]);  

%apply the colormap to the w velocity contour plot
set(gcf,'Colormap',flux_cmap);


figure (fig1)
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');
iname=[dirplots 'Bantry_A_Sections_' datestub '.png'];
print ('-dpng',iname);

figure (fig3)
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');
iname=[dirplots 'Bantry_B_sections/Bantry_B_sections_' datestub '.png'];
print ('-dpng',iname);

close(fig1);
close(fig2);
close(fig3);
    
end