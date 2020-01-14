function gfitgui(varargin)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%


% --- get figure handle
if nargin == 1
    hLine = varargin{1};
    if ~isgraphics(hLine,'line')
        error('Invalid line handle ...');
    end
elseif nargin == 0
    hFig = gcf;
else
    error('Invalid input argument ...');    
end
figure(hFig);

udata.hFig = hFig;
udata.hLines = getlinehandles(hFig);    % class: matlab.graphics.chart.primitive.Line  
%keyboard

return;


function hLines = getlinehandles(hFig)
hLines = findall(hFig,'type','line');