%% This example demonstrate how one can perform multiple linecut and fit
% with LineFit packages.

%% Load parameter file. To creat a parameter file, one can either set the 
% parameters in GUI mode and save the settings to file, or build from
% scratch using class gixsdata.
loaded_params = load('giwaxs_params_20180204.mat');   

%% Give image file
filepath = pwd;     % path for image file
filename = 'H3__th0.130_sec005_up_001_GapFilled.tif';     % image file name

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

%% Decide linecut constraint. One can define a mask to do that, or to be 
% precised use constraint matrix.
% --- linecut paramters
xflag = 1;          % q linecut
nofpts = 300;      % number of points in the linecut
dataflag = 2;       % 2 for corrected data; 1 for masked rawdata

% --- define constraint matrix
constr_q = [1 1 0.71 0.96];     % apply q 
chi_edges = [-16:4:-4, 4:4:16]; % boundaries for chi contraints
nofc = length(chi_edges)-1;     % number of chi values (i.e linecut curves)
chi = nan(nofc,1);              % initialize the list of chi for each linecut
data_linecut = cell(nofc,1);    % initialize a cell to store all linecut
for ii=1:nofc
    constr_chi = [1 9 chi_edges(ii:ii+1)];  % apply chi constraint
    constr = [constr_q; constr_chi];        % combine constraints
    [x,y,~,mapdata] = linecut(obj,xflag,constr,nofpts,dataflag); %
    chi(ii) = mean(mapdata.Data(:,9),'omitnan');   % mean chi value
    data_linecut{ii} = [x,y];
end

%% Plot linecut
figure
hold on;
for ii=1:nofc
    plot(data_linecut{ii}(:,1),data_linecut{ii}(:,2));
end
hold off; box;
legend(num2str(chi(:)));

%% Define fitting model (see LineFit manual for instructions)
objfit = linefit;               % create of linefit object
objfit.BkgdModelIndex = 2;      % use linear background
objfit.CurveModelIndex = 13;    % use voigt model
fresult = repmat(objfit,nofc,1);     % initialize a list to store fit result
for ii=1:nofc
    objfit.Data = data_linecut{ii};
    % --- initialize model parameters
    objfit.Model.BkgdModel.StartParams = [min(objfit.Data(:,2)); 0];            % use min value as background
    line_features = objfit.lineprops(objfit.Data);      % get data features
    objfit.Model.CurveModel.StartParams(1) = line_features.area;        % area
    objfit.Model.CurveModel.StartParams(2) = line_features.peak(1);     % peak position
    objfit.Model.CurveModel.StartParams(3) = line_features.fwhm;        % Gaussian width parameter of voigt
    objfit.Model.CurveModel.StartParams(4) = line_features.fwhm;        % Cauchy width parameter of voigt
    objfit = startfit(objfit);
    fresult(ii) = acceptfit(objfit);
    plot(objfit,'fit');
    pause(0.5);
    %pause
end
% --- Save fit result to MAT file
[~,fnamestem,~]=fileparts(filename);
save(fullfile(filepath,['lfresult_',fnamestem,'.mat']),'fresult','data_linecut','chi');

%% Collect fit result parameters
fresult_params_table = [];
for ii=1:nofc
    fresult_params_table = [...
        fresult_params_table
        struct2table(getmodelprops(fresult(ii)))];
end
fresult_params_table.chi = chi(:);
display(fresult_params_table);
% --- write fit parameters to txt file
writetable(fresult_params_table,fullfile(filepath,['lfresult_',fnamestem,'__params.txt']),'Delimiter','tab');

%% Save raw and fit data to txt file
% One file per curve
for ii=1:nofc
    [y,iy,ibk] = evalmodel(fresult(ii));
    save_data_table = array2table([data_linecut{ii},y,iy,ibk],'VariableNames',{'X','Raw','Fit','FitModel','FitBackground'});
    writetable(save_data_table,fullfile(filepath,['lfresult_',fnamestem,'__',num2str(ii),'.txt']),'Delimiter','tab');
end

