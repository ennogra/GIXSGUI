% ***********************************************
% Copyright (c) 2017 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% This example illustrates the usage of linecut method in the script mode

%% Load parameter file. To creat a parameter file, one can either set the 
% parameters in GUI mode and save the settings to file, or build from
% scratch using class gixsdata.
loaded_params = load('params_demo_linecut.mat');   

%% Give image file
filepath = pwd;     % path for image file
filename = 'image_demo_superlattice.tif';     % image file name

%% Create a new gixsdata object by duplicating the one in loaded_params
obj = copyhobj(loaded_params.params);

%% Load image to obj (two ways to do it)
% --- use file name
obj.ImFile = fullfile(filepath,filename);
% % --- directly pass image data (this method does not keep file names and 
% % image information)
% obj.RawData = imread(fullfile(filepath,filename));

%% Display the image with qz qy axis
obj.PlotAxisLabel = 2; 
figure
imagesc(obj);

%% Decide linecut constraint. One can define a mask to do that or more 
% precised using constraints.
% --- define constraint
constr = [1 3 0.025 0.03];      % use 0.025<qz<=0.03 as the constraint
% --- perform linecut. 
xflag = 5;          % qy linecut
nofpts = 1000;      % number of points in the linecut
dataflag = 2;       % 2 for corrected data; 1 for masked rawdata
[x,y] = linecut(obj,xflag,constr,nofpts,dataflag); %

%% plot linecut
figure
plot(x,y,'bo-')
xlabel('q_y (A^{-1})');
ylabel('Intensity (a.u.)');

%% save linecut to ASCII file
% linecut_data = [x,y];
% save linecut_data.dat linecut_data -ASCII -TABS;

%% Further linecut analysis using your own program or Matlab Curve Fitting toolbox
% cftool(x,y); % to Curve Fitting Toolbox

%% delete gixsdata object
delete(obj)