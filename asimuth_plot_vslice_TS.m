function asimuth_plot_vslice_TS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:     asimuth_plot_vslice_TS
%
% Description:  Plot T & S sections for Bantry mouth and inner Bantry
%
% Inputs:       none
%
% Outputs:      none
%               
% Author:       TD, 2013
% Modified By:  KL, Aug 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirplots='';
dirdata='';

%Input arguments to ultimate function
filenm = 'bantry_bay_avg_0001.nc';

hzslice = -20;   %depth at which w velocities to be plotted in plan view.

%Color ranges for contouring various images.
BantryFlux_cmax = 100; %m3/s

% cross-sections requried (hardwired).
along = [118 168 290 290; 160 189 359 359];

%limits and contours
ncol=120; %T,S,R all have 12 colours
tclim=[8 20];
tcontour=[8:1:20];
sclim=[34.4 35.6];
scontour=[34.4:0.1:35.6];
rclim=[25 27.4];
rcontour=[25:0.2:27.4];

%Get dimensions of model domain from netCDF file.
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
fig1 = figure('Name','Bantry Mouth Salinity','NumberTitle','off','Units','pixels','Color',[0.9 0.9 0.9],'visible','off');
fig2 = figure('Name','Bantry Inner Salinity','NumberTitle','off','Units','pixels','Color',[0.9 0.9 0.9],'visible','off');
fig3 = figure('Name','Bantry Mouth Temperature','NumberTitle','off','Units','pixels','Color',[0.9 0.9 0.9],'visible','off');
fig4 = figure('Name','Bantry Inner Temperature','NumberTitle','off','Units','pixels','Color',[0.9 0.9 0.9],'visible','off');
fig5 = figure('Name','Bantry Mouth Density','NumberTitle','off','Units','pixels','Color',[0.9 0.9 0.9],'visible','off');
fig6 = figure('Name','Bantry Inner Density','NumberTitle','off','Units','pixels','Color',[0.9 0.9 0.9],'visible','off');

% get flux data from netCDF file.
maskr = nc_var_read(filenm,'mask_rho')';
salt=nc_var_read(filenm,'salt');
temp=nc_var_read(filenm,'temp');

salt=permute(salt,[3 2 1]);
temp=permute(temp,[3 2 1]);

%Calculate rho (sigma-t):
rho=rho_pot(temp,salt)-1000.;

%rearrange salt array indices
salt = permute(salt,[1 3 2]);
temp = permute(temp,[1 3 2]);
rho = permute(rho,[1 3 2]);

%for each xi plane cross section
for m = 1:size(along,1);
    
    %extract bathymetry of x-section from [h]
    h_1 = h(along(m,3),along(m,1):along(m,2))';
    %extract v-velocity mask of x-section from [maskr]
    mask = maskr(along(m,3),along(m,1):along(m,2));

    %for each sigma level
    for n=1:srho_len      
        %extract masked v-flux through x-sect from [salt]
        fluxS(n,:) = salt(n,along(m,1):along(m,2),along(m,3)).*mask;
        fluxT(n,:) = temp(n,along(m,1):along(m,2),along(m,3)).*mask;
        fluxR(n,:) = rho(n,along(m,1):along(m,2),along(m,3)).*mask;
        %create array of depths at each grid point in x-sect
        z(n,:) = h_1.*srho(1,n);
        %create array of xi points at each grid point in x-sect
        x(n,:) = along(m,1):along(m,2);
    end
    
    
    %create title text for each x-sect plot
    string = ['Section B' num2str(m)];
    strFCDate=['(Forecast Date: ' datestr(dn1,'dd-mmm-yyyy') ')'];
    
    if(m==1)
        figure(fig1);   %Switch to fig1 plot
        %plot the x-sect
        cmap=cbrewer('div', 'RdBu', ncol);
        colormap(flipud(cmap));
        pcolor(x,z,double(fluxS));
        shading interp;
        hold on;
        [C,h1]=contour(x,z,fluxS,scontour,'k');
        th=clabel(C,h1,'LabelSpacing',1000,'Rotation',0,'FontSize',12);
        for c=1:length(th)
            cl=get(th(c),'UserData');
            cl=char(num2str(cl,'%5.1f'));
            set(th(c),'String', cl);
        end
        %plot bathymetry profile along x-sect
        line(x(1,:),z(1,:),'Color','k','LineWidth',1);hold on;           
        %plot Hzslice line on x-sect
        line([along(m,1) along(m,2)] ,[hzslice hzslice],'LineStyle','--','Color','k','LineWidth',1);hold on;
        set(gca,'XTick',[118 168])
        set(gca,'XTickLabel',['N';'S'])
        set(gca,'YTickLabel',[60 50 40 30 20 10 0])
        %constrain the Bantry mouth x-sect axes
        axis([min(min(x)) max(max(x)) -65 0]);
        %add the title to the Bantry mouth x-sect
        text(along(m,1),z(srho_len,1)+3,[string ': Mouth of Bantry Bay ' strFCDate],'rotation',0,'fontsize',10,'fontweight','bold');
        ylabel('Depth [m]')
        text(along(m,2),z(srho_len,1)+3,['Salinity [psu]'],'rotation',0,'fontsize',12,'fontweight','bold');
        caxis(sclim);
        colorbar('YTick',scontour,'YTickLabel',sprintf('%5.1f|',scontour));
        
        figure(fig3);
        cmap=cbrewer('div', 'RdYlBu', ncol);
        colormap(flipud(cmap));
        pcolor(x,z,double(fluxT));
        shading interp;
        hold on;
        [C,h1]=contour(x,z,fluxT,tcontour,'k');
        th=clabel(C,h1,'LabelSpacing',1000,'Rotation',0,'FontSize',12);
        for c=1:length(th)
            cl=get(th(c),'UserData');
            cl=char(num2str(cl,'%5.1f'));
            set(th(c),'String', cl);
        end
        %plot bathymetry profile along x-sect
        line(x(1,:),z(1,:),'Color','k','LineWidth',1);hold on;           
        %plot Hzslice line on x-sect
        line([along(m,1) along(m,2)] ,[hzslice hzslice],'LineStyle','--','Color','k','LineWidth',1);hold on;
        set(gca,'XTick',[118 168])
        set(gca,'XTickLabel',['N';'S'])
        set(gca,'YTickLabel',[60 50 40 30 20 10 0])
        %constrain the Bantry mouth x-sect axes
        axis([min(min(x)) max(max(x)) -65 0]);
        %add the title to the Bantry mouth x-sect
        text(along(m,1),z(srho_len,1)+3,[string ': Mouth of Bantry Bay ' strFCDate],'rotation',0,'fontsize',10,'fontweight','bold');
        ylabel('Depth [m]')
        text(along(m,2)-2,z(srho_len,1)+3,['Temperature [degC]'],'rotation',0,'fontsize',12,'fontweight','bold');
        caxis(tclim);
        colorbar('YTick',tcontour,'YTickLabel',sprintf('%5.1f|',tcontour));
 
        figure(fig5); 
        cmap=cbrewer('div', 'RdBu', ncol);
        colormap(flipud(cmap));
        pcolor(x,z,double(fluxR));
        shading interp;
        hold on; 
        [C,h1]=contour(x,z,fluxR,rcontour,'k');
        th=clabel(C,h1,'LabelSpacing',1000,'Rotation',0,'FontSize',12);
        for c=1:length(th)
            cl=get(th(c),'UserData');
            cl=char(num2str(cl,'%5.1f'));
            set(th(c),'String', cl);
        end
        %plot bathymetry profile along x-sect
        line(x(1,:),z(1,:),'Color','k','LineWidth',1);hold on;           
        %plot Hzslice line on x-sect
        line([along(m,1) along(m,2)] ,[hzslice hzslice],'LineStyle','--','Color','k','LineWidth',1);hold on;
        set(gca,'XTick',[118 168])
        set(gca,'XTickLabel',['N';'S'])
        set(gca,'YTickLabel',[60 50 40 30 20 10 0])
        %constrain the Bantry mouth x-sect axes
        axis([min(min(x)) max(max(x)) -65 0]);
        %add the title to the Bantry mouth x-sect
        text(along(m,1),z(srho_len,1)+3,[string ': Mouth of Bantry Bay ' strFCDate],'rotation',0,'fontsize',10,'fontweight','bold');
        ylabel('Depth [m]')
        text(along(m,2),z(srho_len,1)+3,['Sigma-t'],'rotation',0,'fontsize',12,'fontweight','bold');
        caxis(rclim);
        colorbar('YTick',rcontour,'YTickLabel',sprintf('%5.1f|',rcontour));
        %KLC: aug 2015 - save section data for future intercomparison
        mfile=[dirdata 'Bantry_mouth_' datestub '_TSR.mat'];
        save(mfile,'dn1','fluxT','fluxS','fluxR');
    end
    
     if(m==2)
        figure(fig2);   
        cmap=cbrewer('div', 'RdBu', ncol);
        colormap(flipud(cmap));
        pcolor(x,z,double(fluxS));
        shading interp;
        hold on;
        [C,h1]=contour(x,z,fluxS,scontour,'k');
        th=clabel(C,h1,'LabelSpacing',1000,'Rotation',0,'FontSize',12);
        for c=1:length(th)
            cl=get(th(c),'UserData');
            cl=char(num2str(cl,'%5.1f'));
            set(th(c),'String', cl);
        end
        %plot bathymetry profile along x-sect
        line(x(1,:),z(1,:),'Color','k','LineWidth',1);hold on;           
        %plot Hzslice line on x-sect
        line([along(m,1) along(m,2)] ,[hzslice hzslice],'LineStyle','--','Color','k','LineWidth',1);hold on;
        set(gca,'XTick',[160 189])
        set(gca,'XTickLabel',['N';'S'])
        set(gca,'YTickLabel',[60 50 40 30 20 10 0])
        %constrain the Bantry mouth x-sect axes
        axis([min(min(x)) max(max(x)) -65 0]);
        %add the title to the Bantry mouth x-sect
        text(along(m,1),z(srho_len,1)+3,[string ': Inner Bantry Bay ' strFCDate],'rotation',0,'fontsize',10,'fontweight','bold');
        ylabel('Depth [m]')
        text(along(m,2),z(srho_len,1)+3,['Salinity [psu]'],'rotation',0,'fontsize',12,'fontweight','bold');
        caxis(sclim);
        colorbar('YTick',scontour,'YTickLabel',sprintf('%5.1f|',scontour));
        
        figure(fig4);
        cmap=cbrewer('div', 'RdYlBu', ncol);
        colormap(flipud(cmap));
        pcolor(x,z,double(fluxT));
        shading interp;
        hold on;
        [C,h1]=contour(x,z,fluxT,tcontour,'k');
        th=clabel(C,h1,'LabelSpacing',1000,'Rotation',0,'FontSize',12);
        for c=1:length(th)
            cl=get(th(c),'UserData');
            cl=char(num2str(cl,'%5.1f'));
            set(th(c),'String', cl);
        end
        %plot bathymetry profile along x-sect
        line(x(1,:),z(1,:),'Color','k','LineWidth',1);hold on;           
        %plot Hzslice line on x-sect
        line([along(m,1) along(m,2)] ,[hzslice hzslice],'LineStyle','--','Color','k','LineWidth',1);hold on;
        set(gca,'XTick',[160 189])
        set(gca,'XTickLabel',['N';'S'])
        set(gca,'YTickLabel',[60 50 40 30 20 10 0])
        %constrain the Bantry mouth x-sect axes
        axis([min(min(x)) max(max(x)) -65 0]);
        %add the title to the Bantry mouth x-sect
        text(along(m,1),z(srho_len,1)+3,[string ': Inner Bantry Bay ' strFCDate],'rotation',0,'fontsize',10,'fontweight','bold');
        ylabel('Depth [m]')
        text(along(m,2)-2,z(srho_len,1)+3,['Temperature [degC]'],'rotation',0,'fontsize',12,'fontweight','bold');
        caxis(tclim);
        colorbar('YTick',tcontour,'YTickLabel',sprintf('%5.1f|',tcontour));

        figure(fig6);
        cmap=cbrewer('div', 'RdBu', ncol);
        colormap(flipud(cmap));
        pcolor(x,z,double(fluxR));
        shading interp;
        hold on;
        [C,h1]=contour(x,z,fluxR,rcontour,'k');
        th=clabel(C,h1,'LabelSpacing',1000,'Rotation',0,'FontSize',12);
        for c=1:length(th)
            cl=get(th(c),'UserData');
            cl=char(num2str(cl,'%5.1f'));
            set(th(c),'String', cl);
        end
        %plot bathymetry profile along x-sect
        line(x(1,:),z(1,:),'Color','k','LineWidth',1);hold on;           
        %plot Hzslice line on x-sect
        line([along(m,1) along(m,2)] ,[hzslice hzslice],'LineStyle','--','Color','k','LineWidth',1);hold on;
        set(gca,'XTick',[160 189])
        set(gca,'XTickLabel',['N';'S'])
        set(gca,'YTickLabel',[60 50 40 30 20 10 0])
        %constrain the Bantry mouth x-sect axes
        axis([min(min(x)) max(max(x)) -65 0]);
        %add the title to the Bantry mouth x-sect
        text(along(m,1),z(srho_len,1)+3,[string ': Inner Bantry Bay ' strFCDate],'rotation',0,'fontsize',10,'fontweight','bold');
        ylabel('Depth [m]')
        text(along(m,2),z(srho_len,1)+3,['Sigma-t'],'rotation',0,'fontsize',12,'fontweight','bold');
        caxis(rclim);
        colorbar('YTick',rcontour,'YTickLabel',sprintf('%5.1f|',rcontour));
        %KLC: aug 2015 - save section data for future intercomparison
        mfile=[dirdata 'Bantry_inner_' datestub '_TSR.mat'];
        save(mfile,'dn1','fluxT','fluxS','fluxR');
    end
     
     %clear out variables        
    clear mask fluxS fluxT fluxR z x;
end
%All xi-section x-sects now plotted


figure (fig1)
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');
print ('-dpng', [dirplots '/Bantry_Mouth_Section_S_' datestub '.png']);
filename1=['Bantry_Mouth_Section_S_' datestub '.png'];

figure (fig2)
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');
print ('-dpng', [dirplots 'Bantry_Inner_Section_S_' datestub '.png']);
filename2=['Bantry_Inner_Section_S_' datestub '.png'];

figure (fig3)
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');
print ('-dpng', [dirplots 'Bantry_Mouth_Section_T_' datestub '.png']);
filename3=['Bantry_Mouth_Section_T_' datestub '.png'];

figure (fig4)
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');
print ('-dpng', [dirplots 'Bantry_Inner_Section_T_' datestub '.png']);
filename4=['Bantry_Inner_Section_T_' datestub '.png'];

figure (fig5)
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');
print ('-dpng', [dirplots 'Bantry_Mouth_Section_RHO_' datestub '.png']);
filename5=['Bantry_Mouth_Section_RHO_' datestub '.png'];

figure (fig6)
set(gcf,'PaperPositionMode','auto','Renderer','painters')
set(gcf,'InvertHardcopy','off');
print ('-dpng', [dirplots 'Bantry_Inner_Section_RHO_' datestub '.png']);
filename6=['Bantry_Inner_Section_RHO_' datestub '.png'];

close(fig1);
close(fig2);
close(fig3);
close(fig4);
close(fig5);
close(fig6);

end
    
