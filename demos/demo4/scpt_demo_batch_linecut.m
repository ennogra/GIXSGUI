% ***********************************************
% Copyright (c) 2017 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% This example illustrates the usage of gixsdata.linecut function to do
% linecut on a batch of the images.


%% collect names of the image files
filepath = pwd;         % give the fold of the image files
filename = dir('giwaxs1_s087_*.tif'); % specify the template of your images
filename = {filename(:).name}; 
filename = filename(:);
filename = sort(filename);  % the total number of images should match that of th angles.

%% define the alpha_i (incident angle) points; make sure the incident angles match the
% sequence of filenames
alpha_i = 0.8052:0.1:2.4053;     % angles in deg;


%% get incident angles
scanvarlist = alpha_i(:);       % assign values of the scanned variable
scanvarname = 'IncidentAngle';  % scanned variable name (to math the field in gixsdata

%% linecut parameters
% --- gixsdata parameter file
loaded_params = load('params_demo_roiscan');
params = loaded_params.params;
xflag = 1;              % 1-11. it is 1 if the linecut is against q
dataflag = 2;           %1/2: use MaskedData or Solidanglecorreted data (keep it for corrected data)
nofp    = 100;          % number of points in the linecut
% --- constraint for linecut
            %       1 - q (total q)
            %       2 - phi (azimuthal angle with respect to direct beam)
            %       3 - qz (for reflection only)
            %       4 - qx (for reflection only)
            %       5 - qy (for reflection only)
            %       6 - qr (parallel q in surface plane; for reflection only)
            %       7 - 2Theta (scattering angle in surface plane; for
            %           reflection only)
            %       8 - Alpha_f (exit angle; for reflection only)
            %       9 - Chi (polar angle; for reflection only)
            %       10 - x pixel
            %       11 - y pixel
const = cell(1,length(filename));       % initialize constraint
figure
for ii=1:length(filename)
   % constraint for each image;
    const{ii} = [...
        1 9  -90           -80
        1 1  0.5           2.2];
    % perform linecut
    params.IncidentAngle = alpha_i(ii);
    params.ImFile = fullfile(filepath,filename{ii});
    [xdata,ydata] = linecut(params,xflag,const{ii},nofp,dataflag);
    % plot 
    plot(xdata,ydata,'o-');
    title(titlestr(filename{ii}));
    % save data
    savefilename = ['linecut_',filename{ii}(1:end-4),'.dat'];
    filesave = fullfile(filepath,savefilename);
    datasave = [xdata,ydata];
    save(filesave,'datasave','-ASCII','-TABS');
    disp(['Linecut data for ',filename{ii},' is saved to ',savefilename]);
    %pause;      % press anykey (e.g. space bar) to continue
    pause(0.2);
end

%% clean up
clear params;