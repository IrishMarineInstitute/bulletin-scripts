function asimuth_plot_vslice_fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:     asimuth_plot_vslice_fluxes
%
% Description:  This function is similar to Plot_2D_fluxes_plan.m, but it 
%               will only plot fluxes through the Bantry mouth and inner 
%               Bantry, both on separate figures. It is kept separate from 
%               the above function so that the other one can be easily 
%               switched on/off if required
%
%
% Inputs:       none
%
% Outputs:      none
%               
% Author:       TD, 2013
% Modified By:  KL, Aug 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE IMPORTANT!!!!!:  If Plot_2D_fluxes_plan.m is ever switched off then
% archiving data to mfiles has to be copied over to this function

%base location for writing plots
dirplots='';
dirdata='';

%Input arguments to ultimate function
filenm = 'bantry_bay_avg_0001.nc';

load 'inflow_Bantry_mouth.mat';
load 'inflow_Bantry_inner.mat';

%Define matfiles
mfile1='inflow_Bantry_mouth.mat';
mfile2='inflow_Bantry_inner.mat';

hzslice = -20;   %depth at which w velocities to be plotted in plan view.

%Color ranges for contouring various images.
BantryFlux_cmax = 100; %m3/s

% cross-sections requried (hardwired).
along = [118 169 290 290; 160 190 360 360];

%Get dimensions of model domainfrom netCDF file.
srho_len=nc_dimlen(filenm,'s_rho');

% get general data from netCDF file.
%note: this script was originally written using Mex netcdf functions which
%returns matrices cols and ros opposite to Matlab native function. Rather
%than rewriting all the code later on, we transpose the arrays here to have
%in same orientation as previous
h=nc_var_read(filenm,'h')';
srho=nc_var_read(filenm,'s_rho')';
ot=nc_var_read(filenm,'ocean_time')';

%Below is for defining string to uniquely name the figure (need to subtract
%86400 from ot(1) as the averaging time period is 67.5hrs and thus the
%ot(1) being in the middle is equal 07:10hrs on the next day
dn1=floor(roms_time_to_datenum(ot(1)-86400));
datestub=datestr(dn1,'yyyymmdd');

%create figure 1 and figure 2 for plotting into
fig1 = figure('Name','Bantry Mouth Flux','NumberTitle','off','Units','pixels','Color',[0.9 0.9 0.9],'visible','off');
fig2 = figure('Name','Bantry Inner Flux','NumberTitle','off','Units','pixels','Color',[0.9 0.9 0.9],'visible','off');

%figures showing the timeseries of inflows:
fig3 = figure('Name','Inflow to Bantry Bay - time series','NumberTitle','off','Units','pixels','Color',[0.9 0.9 0.9],'visible','off');
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');

%load default colormap (red -> blue)
colormap('default');
load('FluxColorMap','flux_cmap');

% get flux data from netCDF file.
maskv=nc_var_read(filenm,'mask_v')';
hv=nc_var_read(filenm,'Hvom');

%rearrange hvom array indices taking into account using matlab netcdf funcs
hvom = permute(hv,[3 1 2]);

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
    sum_net = sum_pos+sum_neg;
    
    if(m==1)
        %Calculate how the current inflow relates to long term mean [%]:
        mean_sum_pos=mean(sum_pos1);
        rel_sum_pos=round(((sum_pos - mean_sum_pos)/mean_sum_pos)*100); 
        %Concatenate
        otime1=[otime1;ot];        
        sum_pos1=[sum_pos1;sum_pos];
        sum_neg1=[sum_neg1;sum_neg];
        sum_net1=[sum_net1;sum_net];
        save(mfile1,'otime1','sum_pos1','sum_neg1','sum_net1');
        %KLC 12 aug 2015 - save flux data too
        mfile=[dirdata  'Bantry_mouth_' datestub '_flux.mat'];
        save(mfile,'dn1','flux');
    end
    
    if(m==2)
        %Calculate how the current inflow relates to long term mean [%]:
        mean_sum_pos=mean(sum_pos2);
        rel_sum_pos=round(((sum_pos - mean_sum_pos)/mean_sum_pos)*100);
        %Concatenate
        otime2=[otime2;ot];        
        sum_pos2=[sum_pos2;sum_pos];
        sum_neg2=[sum_neg2;sum_neg];
        sum_net2=[sum_net2;sum_net];
        
        save(mfile2,'otime2','sum_pos2','sum_neg2','sum_net2');
        %KLC 12 aug 2015 - save flux data too
        mfile=[dirdata  'Bantry_inner_' datestub '_flux.mat'];
        save(mfile,'dn1','flux');
    end
    
    %create title text for each x-sect plot
    string = ['Section B' num2str(m)];
    strFCDate=['(Forecast Date: ' datestr(dn1,'dd-mmm-yyyy') ')'];
    
    if(m==1)
        figure(fig1);   %Switch to fig1 plot
        %plot the x-sect
        pcolor(x,z,double(flux));shading interp;hold on;       
        %plot bathymetry profile along x-sect
        line(x(1,:),z(1,:),'Color','k','LineWidth',1);hold on;           
        %plot Hzslice line on x-sect
        line([along(m,1) along(m,2)] ,[hzslice hzslice],'LineStyle','--','Color','k','LineWidth',1);hold on;
        set(gca,'XTick',[118 169])
        set(gca,'XTickLabel',['N';'S'])
        set(gca,'YTickLabel',[60 50 40 30 20 10 0])
        %constrain the Bantry mouth x-sect axes
        axis([min(min(x)) max(max(x)) -65 0]);
        %add the title to the Bantry mouth x-sect
        text(along(m,1),z(srho_len,1)+3,[string ': Mouth of Bantry Bay ' strFCDate],'rotation',0,'fontsize',10,'fontweight','bold');
        ylabel('Depth [m]')
        text(along(m,2),z(srho_len,1)+3,['Flow [m3/s]'],'rotation',0,'fontsize',12,'fontweight','bold');
        text(along(m,2)+7,z(srho_len,1)-10,['IN'],'rotation',0,'fontsize',15,'fontweight','bold');
        text(along(m,2)+7,z(srho_len,1)-50,['OUT'],'rotation',0,'fontsize',15,'fontweight','bold');
        %Write how the flow relates to long term mean:
        if(rel_sum_pos<0)
            rel_sum_pos=abs(rel_sum_pos);
            text(along(m,1)+3,z(srho_len,1)-62,['Current inflow is lower by ',num2str(rel_sum_pos,'%5.1u'),'% than average'],'rotation',0,'fontsize',13,'fontweight','bold','Color',[0 .5 0]);
        else
            text(along(m,1)+3,z(srho_len,1)-62,['Current inflow is greater by ',num2str(rel_sum_pos,'%5.1u'),'% than average'],'rotation',0,'fontsize',13,'fontweight','bold','Color','r');
        end
        %constrain the x-sect color ramping
        caxis([-BantryFlux_cmax BantryFlux_cmax]);
        colorbar('Ytick',[-BantryFlux_cmax 0 BantryFlux_cmax]);
        %apply the colormap to the x-sect v-fluxes
        set(gcf,'Colormap',flux_cmap);
    end
    if(m==2)
        figure(fig2);   %Switch to fig2 plot
        %plot the x-sect
        pcolor(x,z,double(flux));shading interp;hold on;       
        %plot bathymetry profile along x-sect
        line(x(1,:),z(1,:),'Color','k','LineWidth',1);hold on;           
        %plot Hzslice line on x-sect
        line([along(m,1) along(m,2)] ,[hzslice hzslice],'LineStyle','--','Color','k','LineWidth',1);hold on;
        set(gca,'XTick',[160 190])
        set(gca,'XTickLabel',['N';'S'])
        set(gca,'YTickLabel',[60 50 40 30 20 10 0])
        %constrain the Bantry mouth x-sect axes
        axis([min(min(x)) max(max(x)) -65 0]);
        %add the title to the Bantry mouth x-sect
        text(along(m,1),z(srho_len,1)+3,[string ': Inner Bantry Bay ' strFCDate],'rotation',0,'fontsize',10,'fontweight','bold');
        ylabel('Depth [m]')
        text(along(m,2),z(srho_len,1)+3,['Flow [m3/s]'],'rotation',0,'fontsize',12,'fontweight','bold');
        text(along(m,2)+5,z(srho_len,1)-10,['IN'],'rotation',0,'fontsize',15,'fontweight','bold');
        text(along(m,2)+5,z(srho_len,1)-50,['OUT'],'rotation',0,'fontsize',15,'fontweight','bold');
        %Write how the flow relates to long term mean:
        if(rel_sum_pos<0)
        rel_sum_pos=abs(rel_sum_pos);
        text(along(m,1)+1,z(srho_len,1)-62,['Current inflow is lower by ',num2str(rel_sum_pos,'%5.1u'),'% than average'],'rotation',0,'fontsize',13,'fontweight','bold','Color',[0 .5 0]);
        else
        text(along(m,1)+1,z(srho_len,1)-62,['Current inflow is greater by ',num2str(rel_sum_pos,'%5.1u'),'% than average'],'rotation',0,'fontsize',13,'fontweight','bold','Color','r');
        end
        %constrain the x-sect color ramping
        caxis([-BantryFlux_cmax BantryFlux_cmax]);
        colorbar('Ytick',[-BantryFlux_cmax 0 BantryFlux_cmax]);
        %apply the colormap to the x-sect v-fluxes
        set(gcf,'Colormap',flux_cmap);
    end
    
        
    %clear out variables        
    clear mask flux z x;
end
%All xi-section x-sects now plotted

%Now plot the timeseries
%get the date:
dn=roms_time_to_datenum(otime1);
figure(fig3)
plot(dn,sum_pos1,dn,sum_pos2);hold on;
datetick('x',19);
xlabel('Date')
ylabel('Flow [m3/s]')
legend('Bantry Bay mouth','Bantry Bay inner')

%clear arrays 
clear mask_v hv;

figure (fig1)
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');
print ('-dpng', [dirplots 'Bantry_Mouth_Section_' datestub '.png']);
filename1=['Bantry_Mouth_Section_' datestub '.png'];

figure (fig2)
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');
print ('-dpng', [dirplots 'Bantry_Inner_Section_' datestub '.png']);
filename2=['Bantry_Inner_Section_' datestub '.png'];

figure (fig3)
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');
print ('-dpng', [dirplots 'Bantry_Inflow_Timeseries.png']);
filename3=['Bantry_Inflow_Timeseries.png'];

close(fig1);
close(fig2);
close(fig3);


end