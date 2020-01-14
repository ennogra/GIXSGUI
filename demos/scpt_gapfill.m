% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% this code demonstrate how to gapfill in the script batch mode

fpath_summed =  'N:\2014-2\wchen_giwaxs_201406p1\summed';
fpath_gapfilled =  'N:\2014-2\wchen_giwaxs_201406p1\gapfilled';

flist = dir(fullfile(fpath_summed,'*.tif'));
flist = {flist().name}';

for ii = 1:2; %numel(flist)/2
    fup = fullfile(fpath_summed,flist{ii*2});
    fdn = fullfile(fpath_summed,flist{ii*2-1});
    im = single(mgapfill(fup,fdn));
    info = imfinfo(fup);
    fexport = ['gapfilled_',flist{ii*2}]
    % --- get img type to save
    hFigGapFill = findall(0,'Tag','gapfill_fig');
    udata = get(hFigGapFill,'UserData');
    hPopupmenuImageType = findall(hFigGapFill,'tag','gapfill_PopupmenuImageType');
    type_str = get(hPopupmenuImageType,'string');
    type_value = get(hPopupmenuImageType,'value');
    % --- save
    imwrite2tif(im,info,fullfile(fpath_gapfilled,fexport),type_str{type_value});
end

