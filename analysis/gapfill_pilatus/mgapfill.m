function imfdata = mgapfill(fup,fdn)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% MGAPFILL Manually gapfill up and down images.
%   IMFDATA = MGAPFILL(FUP,FDN) FUP and FDN are full file names of the up 
%       and down images. IMFDATA is the ouput gapfilled image. MGAPFILL
%       uses current settings in GAPFILL GUI window; otherwise, it uses
%       GAPFILL default settings.
%
%   Zhang Jiang @8ID/APS/ANL
% $Revision: 1.0 $  $Date: 2014/10/28 $


hFigGapFill = findall(0,'Tag','gapfill_fig');
if ~isempty(hFigGapFill)
    figure(hFigGapFill);
else
    gapfillgui;
end


hFigGapFill = findall(0,'Tag','gapfill_fig');
udata = get(hFigGapFill,'UserData');

% --- load base image
hEditBase = findall(hFigGapFill,'tag','gapfill_EditBase');
set(hEditBase,'string',fup);
feval(get(hEditBase,'callback'));

% --- load image1
hEditImage1 = findall(hFigGapFill,'tag','gapfill_EditImage1');
set(hEditImage1,'string',fdn);
feval(get(hEditImage1,'callback'),[],[],1);

% --- fill
feval(get(findall(hFigGapFill,'tag','gapfill_PushbuttonFill'),'callback'));

% --- output and close figure
udata = get(hFigGapFill,'UserData');
imfdata = udata.imfdata;
feval(get(findall(hFigGapFill,'tag','gapfill_PushbuttonCloseAll'),'callback'));

