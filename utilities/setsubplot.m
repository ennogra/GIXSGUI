function setsubplot(varargin)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% SETSUBPLOT Set the properties of all the axes in the figure.
%   SETSUBPLOT(H,'PropertyName',PropertyValue) sets the valus of the
%   specified properties for all the axes in the figure with handle H. H
%   can be vector.

% Zhang Jiang $2006/07/06$

if nargin < 3 || mod(nargin+1,2) == 1
    error('Wrong input arguments.');
    return;
end

h = varargin{1};
h = h(:);
for ii=1:length(h)
    if ~ishandle(h(ii)) || ~strcmpi(get(h(ii),'type'),'figure')
        error('Invalid figure handle.');
        return;
    end
end

h_axes = findall(h,'type','axes','visible','on');
if isempty(h_axes)
    error('No plot exist in the figure.');
    return;
end

for ii=1:length(h_axes)
    for i=2:2:nargin
        set(h_axes(ii),varargin{i},varargin{i+1});
    end
end