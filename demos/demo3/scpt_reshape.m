% ***********************************************
% Copyright (c) 2017 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% This example illustrates the usage of reshape_image method in the script mode

%% Load parameter file. To creat a parameter file, one can either set the 
% parameters in GUI mode and save the settings to file, or build from
% scratch using class gixsdata.
loaded_params = load('params_demo_reshape.mat');   

%% Give image file
filepath = pwd;     % path for image file
filename = 'image_demo_raw_gapfilled.tif';     % image file name

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

%% Reshape
% --- Define reshaping parameter for (qz vs qr) 
param_reshape.X = 6;  % qr
param_reshape.Y = 3;  % qz
param_reshape.XNOfPts = 600; % number of points for X axis
param_reshape.YNOfPts = 600; % number of points for y axis
param_reshape.XRange = [-2,0.5];  % range for x
param_reshape.YRange = [0,2];  % range for y

% % --- Define reshaping parameter for (qz vs chi) 
% param_reshape.X = 9;  % chi
% param_reshape.Y = 3;  % qz
% param_reshape.XNOfPts = 300; % number of points for X axis
% param_reshape.YNOfPts = 300; % number of points for y axis
% param_reshape.XRange = [-90,90];  % range for x
% param_reshape.YRange = [0,2];  % range for y

% --- reshape 
dataflag = 2;       % 2 for corrected data; 1 for masked rawdata
[x,y,img_reshaped,countdata] = reshape_image(obj,param_reshape,dataflag);

% --- plot reshaped image
figure
imagesc(x,y,log10(img_reshaped),[1,4]);
axis ij;
set(gca,'ydir','norm');
xlabel('q_r (A^{-1})');
ylabel('q_z (A^{-1})');
title('Reshaped image (log scale)');

% --- plot counter data
figure
imagesc(x,y,countdata)
axis ij
set(gca,'ydir','norm');
xlabel('q_r (A^{-1})');
ylabel('q_z (A^{-1})');
title('Counter for reshaped image')
%% save reshaped image to tif file
% imwrite2tif(img_reshaped,[],'your_image_file_name_here.tif','single');

%% delete gixsdata object
delete(obj)