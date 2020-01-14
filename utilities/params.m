function [scanPeak,scanCOM,scanFWHM] = params(scan)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% PARAMS Calculate peak,COM and FWHM of a scan
%
% Copyright 2004, Zhang Jiang

% --- peak
[peak,peakIndex] = max(scan(:,2));
scanPeak.X = scan(peakIndex,1);
scanPeak.Y = peak;

% --- COM
scanCOM = trapz(scan(:,1),scan(:,1).*scan(:,2))/trapz(scan(:,1),scan(:,2));

% --- FWHM
try
    leftError = 0;
    halfLeftIndex = find(scan(1:peakIndex,2)<peak/2);
    left_scan = scan(halfLeftIndex(end):peakIndex,:);
    [~,m,~] = unique(left_scan(:,2));
    left_scan = left_scan(m,:);
    left_scan = sortrows(left_scan,1);
    left = interp1(...
        left_scan(:,2),left_scan(:,1),peak/2);
catch
    leftError = 1;
end
try
    rightError = 0;
    halfRightIndex = find(scan(peakIndex:end,2)<peak/2)+peakIndex-1;
    right_scan = scan(peakIndex:halfRightIndex(1),:);
    [~,m,~] = unique(right_scan(:,2));
    right_scan = right_scan(m,:);
    right_scan = sortrows(right_scan,1);  
    right = interp1(...
        right_scan(:,2),right_scan(:,1),peak/2);
catch
    rightError = 1;
end
if leftError == 0 && rightError == 0
    scanFWHM.center = (left+right)/2;
    scanFWHM.FWHM = right-left;
elseif leftError == 0 && rightError == 1
    scanFWHM.center = left;
    scanFWHM.FWHM = scan(end,1)-left;
elseif leftError == 1 && rightError == 0
    scanFWHM.center = right;
    scanFWHM.FWHM = right-scan(1,1);
else
    scanFWHM.center = (scan(end,1)+scan(1,1))/2;
    scanFWHM.FWHM = scan(end,1)-scan(1,1);
end