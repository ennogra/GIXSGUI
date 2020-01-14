function gapmask = gapmask_pilatus1mf(dm,dn)
% ***********************************************
% Copyright (c) 2017 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% GAPMASK_PILATUS Create gap mask for Pilatus 1MF.
%   GAPMASK = GAPMASK_PILATUS1MF(DM,DN) returns a new gap mask with 2*DM 
%   more gap rows and 2*DN more gap columns added to default gap mask.
%
%   Use IMWRITE2TIF to export GAPMASK to TIF file.Specify 'logical' as the
%   export data type.

%   Zhang Jiang and Joe Strzalka @8ID/APS/ANL
%   $Revision: 1.0 $  $Date: 2011/03/02 $

% default (hardware) gap
gaprows = [196 212;
           408 424;
           620 636;
           832 848]; % 4 gaps of 17 pixels
gapcols = [488 494];    % 1 gap of 7 pixels

% create gap mask with all true elements
m = 1043;
n = 981;
gapmask = true(1043,981);

% add extra gap
gaprows(:,1) = gaprows(:,1)-dm;
gaprows(:,2) = gaprows(:,2)+dm;
gapcols(:,1) = gapcols(:,1)-dn;
gapcols(:,2) = gapcols(:,2)+dn;
gaprows(gaprows<1) = 1;
gaprows(gaprows>m) = m;
gapcols(gapcols<1) = 1;
gapcols(gapcols>n) = n;

% creat new gap mask
for ii=1:size(gaprows,1)
    gapmask(gaprows(ii,1):gaprows(ii,2),:) = 0;
end
for ii=1:size(gapcols,1);
    gapmask(:,gapcols(ii,1):gapcols(ii,2)) = 0;
end

