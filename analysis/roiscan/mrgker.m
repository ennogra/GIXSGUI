function mergedata = mrgker(varargin)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% MRGKER Merge kernel function.
%   MERGEDDATA = MRGKER(SCANDATA) 
%   MERGEDDATA = MRGKER(SCANDATA,OPTS)
%
%   Input format:
%       SCANDATA must be in a cell structure with format of each element:
%           col_1      col_2      col_3        
%           xdata     ydata    ydataError(absolute error)
%       OPTS is a structure to define the merging options. It has two
%           fields:
%           OPTS.interpMethod: 1/2/3/4 for spline/linear/nearest/cubic
%           OPTS.mode: 1/2 for intensity based or number of points based
%               vertical shiftting of the data sections.
%       
%   Output format:
%       N x 3 matrix
%
% Copyright 2004, Zhang Jiang

scanData = varargin{1};
scanData = scanData(:)';
merge = varargin{2};

% --- determine merging method
switch merge.interpMethod
    case 1
        interpMethod = 'spline';
    case 2
        interpMethod = 'linear';
    case 3
        interpMethod = 'nearest';
    case 4
        interpMethod = 'cubic';
end

% if the 3d column (error) does not exit, add zeros
for iScan = 1:length(scanData)
    if size(scanData{iScan},2) == 2
        scanData{iScan}(:,3) = 0;
    end
end 

% --- convert the third column to relative error because in merging scans,
% the relatie error of each point is unchanged, although the absolute
% error changs with the merging factors
for iScan = 1:length(scanData)
    scanData{iScan}(:,3) = scanData{iScan}(:,3)./scanData{iScan}(:,2);
end

% --- sort scans in the order of increasing x-axis values
for iScan=1:length(scanData)
    scanData{iScan} = sortrows(scanData{iScan},1);
end
for iScan=1:length(scanData)
    firstAngle(iScan) = scanData{iScan}(1,1);
end
tempData = scanData;
for iScan=1:length(scanData)
    minAngleIndex = find(firstAngle==min(firstAngle));
    tempData{1,iScan} = scanData{1,minAngleIndex};
    firstAngle(minAngleIndex)=inf;
end

% --- remove extra data points that are at the same x positions
for iScan=1:length(tempData)
    [b,m,n] = unique(tempData{iScan}(:,1));
    tempData{iScan} = tempData{iScan}(m,:);    
end

% --- merge scanData
data=tempData{length(tempData)};    % start merging from rightmost scan
for iDataSet=length(tempData):-1:2
    left=tempData{iDataSet-1};      % get left-side scan
    right=data;                     % get right-side scan
    if min(right(:,1)) < max(left(:,1))     % if there is an overlapping region
        right_index = find(right(:,1)<=max(left(:,1)));
        left_index = find(left(:,1)>=right(1,1));        
        interp_x = (linspace( right(1,1),left(end,1),max(length(right_index),length(left_index))*5 ))';
        right_y = interp1(right(:,1),right(:,2),interp_x,interpMethod);
        left_y = interp1(left(:,1),left(:,2),interp_x,interpMethod); 
        num_points = length(interp_x);
        right_area = (sum(right_y)-(right_y(1)+right_y(num_points))/2.0)*(interp_x(num_points)-interp_x(1))/(num_points-1);
        left_area = (sum(left_y)-(left_y(1)+left_y(num_points))/2.0)*(interp_x(num_points)-interp_x(1))/(num_points-1);
        ratio = right_area/left_area;
        if ratio >= 1           % use intensity based
            left(:,2)   = left(:,2)*ratio;        
        else
            right(:,2)  = right(:,2)/ratio;
        end
        removeFlag = [];        % 1: remove left 2:remove right
        switch merge.mode
            case 1      % intensity based: for the overlapping part, keep the one with larger intensity
                if ratio >= 1, removeFlag = 1; else removeFlag = 2; end 
            case 2      % # of point based: keep the one with more data points. But if both sides have the same number use intensity based
                if length(right_index) > length(left_index),        removeFlag = 1;               
                elseif length(right_index) < length(left_index),    removeFlag = 2;
                else
                    if ratio >= 1,  removeFlag = 1; else removeFlag = 2; end
                end
        end
        switch removeFlag
            case 1, left(left_index(1):end,:)=[];
            case 2, right(1:right_index(end),:)=[];
        end
        data=[left;right];
    elseif min(right(:,1)) == max(left(:,1))        % if there is only one overlapping point
        ratio = right(1,2)/left(end,2);             
        if ratio >= 1       % if right side has larger intensity
            left(:,2)   = left(:,2)*ratio;
            left(end,:) = [];
        else                % left side has larger intensity        
            right(:,2) = right(:,2)/ratio;
            right(1,:) = [];
        end
        data=[left;right];        
    else                                            % if there is no overlapping point, just merge assuming original ratio is 1
        ratio = right(1,2)/tempData{iDataSet}(1,2); % But since right has been modified, the ratio has to be recalculated)
        left(:,2)=left(:,2)*ratio;
        data=[left; right];
%        uiwait(msgbox('Curves do not overlap ! Will be merged anyway. ','Merge Warning','warn','modal'));
    end
end
data(:,3) = data(:,3).*data(:,2);       % convert the relative error back to absolute error
mergedata = data;