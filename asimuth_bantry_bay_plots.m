%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script:       asimuth_bantry_bay_plots
%
% Description:  This is the control script for creating plots from Bantry
%               Bay model FORECAST runs for the ASIMUTH project
%
% Inputs:       none
%  
% Outputs:      none
%               
% Author:       TD, 2013
% Modified By:  KL, Aug 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
aname='bantry_bay_avg_0001.nc';
%
%
% Call function to create section plot (volumetric flows)
% NOTE: No need to pass the filename, as ..._avg_0001.nc is always used
%
asimuth_plot_2D_fluxes_plan;

% This one creates Bantry Mouth and Inner sections on two separate plots
asimuth_plot_vslice_fluxes;

% This one creates Bantry Mouth and Inner sections through T & S on two separate plots
asimuth_plot_vslice_TS;

%KLC: 14/8/2015 Create cross-shelf section plots
asimuth_cross_shelf_xsect(aname);

% Now call the script to plot particle-hours (requires passing on the
% proper filename)
% Generate the filename here:

%get the current float file name from the postproc.in file
infile='bantry_post_process_fc.in';
[varid,param_name]=textread(infile,'%s %s','delimiter',',','commentstyle','matlab');
for j=1:length(varid)
    switch varid{j}
        case 'float_file'
            fname=param_name{j};
    end;
end;

if exist(fname,'file')==2
    asimuth_plot_floats(fname);
    disp('Now create Netcdf file');
    asimuth_bantry_particle_stats(fname);
else
    error('Bantry float file does not exist');
end

%KLC- Sep 2016
%copy and rename files into WEEK_ARCHIVE
[fair,fwind,sfile,ifirstfile,nfiles] = get_FC_model_run_details(infile);
indir='';
outdir='';
ofstub='bantry_bay_his_';
nfstub='BANT';
ilastfile=ifirstfile+nfiles-1;
roms_files_copy_rename(indir,outdir,ifirstfile,ilastfile,ofstub,nfstub);

quit;
