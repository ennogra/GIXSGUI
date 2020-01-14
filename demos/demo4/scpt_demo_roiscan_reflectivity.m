% ***********************************************
% Copyright (c) 2017 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% This example illustrates the usage of roiscan function to do ROI scan

%% define the alpha_i (incident angle) points in the image scan
alpha_i = 0.8052:0.1:2.4053;     % angles in deg

%% collect names of the image files
filepath = pwd;         % give the fold of the image files
filename = dir('giwaxs1_s087_*.tif'); 
filename = {filename(:).name}; 
filename = filename(:);
filename = sort(filename);  % the total number of images should match that of th angles.

%% get incident angles
scanvarlist = alpha_i(:);            % assign values of the scanned variable
scanvarname = 'IncidentAngle';  % scanned variable name (to math the field in gixsdata

%% roiscan parameters
% --- gixsdata parameter file
loaded_params = load('params_demo_roiscan');
params = loaded_params.params;
% --- constraint for roiscan
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
const = cell(1,length(filename));
% constraints (a cell) are identical for all angles in rocking scan. In
% reflectivity scan, the constraint floats and Alphaf in each constraint 
% is twice of the IncidentAngle
for ii=1:length(filename)
   % constraint for reflectivity: define a virtual detector so that 
   % alpha_f is equal to incident angle;
    const{ii} = [...
        1 8  alpha_i(ii)-0.1   alpha_i(ii)+0.1
        1 7  -0.5           0.5];
   % constraint for rocking scan: define a virtual detector so that 
   % alpha_i+alpha_f = constant 
    %     const{ii} = [...
%         1 8  ???   ???;
%         1 7  ???   ???];  
end
% --- data flag
dataflag = 2;           %1/2: use MaskedData or Solidanglecorreted data
plotflag = 1;           % 0/1 plot during ROI scan?
cmaskdisplayflag = 1;   % 0/1 display cmark box during ROI scan?

%% start roiscan
scanresult = roiscan(filepath,filename,scanvarname,scanvarlist,params,const,dataflag,plotflag,cmaskdisplayflag);


%% plot scan
figure
plot(scanresult.X,scanresult.I,'bo-');
xlabel('IncidentAngle');
ylabel('Intensity');
set(gca,'yscale','log');

%% delete params
delete params;