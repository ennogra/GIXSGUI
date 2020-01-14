function [hImage,hAxes,hFig] = gu_imhandles(h)
%IMHANDLES Get all image handles.  
%   HIMAGE = IMHANDLES(H) takes a graphics handle H as an input and
%   returns all of the image handles whose ancestor is H. H can
%   be an array of valid figure, axes, image, or uipanel handles.
%
%   HIMAGE is an array of image handles. 
%
%   IMHANDLES ignores a colorbar in H and does not include its handle in
%   HIMAGE.
%
%   Note
%   ----
%   IMHANDLES errors if the image objects in HIMAGE do not have the same
%   figure as their parent.
%
%   Examples
%   --------
%       figure,imshow('moon.tif');
%       hImage = imhandles(gca)
%
%       subplot(1,2,1),imshow('autumn.tif');
%       subplot(1,2,2),imshow('glass.png');
%       hImage = imhandles(gcf)
%
%   See also IMGCA, IMGCF.

%   These undocumented syntaxes may change in the future.  
%   [HIMAGE,HAXES] = IMHANDLES(H) returns an array of axes handles
%   corresponding to HIMAGE. HAXES is a vector of axes handles.
%
%   [HIMAGE,HAXES,HFIG] = IMHANDLES(H) returns a figure handle that
%   contains HIMAGE and HAXES.

%   Copyright 1993-2010 The MathWorks, Inc.
%   $Revision: 1.1.8.7 $  $Date: 2010/04/15 15:18:30 $

if ishghandle(h)
  type = get(h,'Type');
  
  % correctHGobj will be an array contain true if type is
  % part of the valid list and false if type is not part of valid it.
  correctHGobj = ismember(type, ...
  {'image','figure','uicontainer','axes','uipanel'});
  
  % correctHGobj would have the same length as h, which can be > 1. 
  if all(correctHGobj)
    hImage = findobj(h,'Type','image');
    
    % remove the colorbars contained in hImage
    colorbarHandle = findobj(hImage, 'Type', 'Image', ...
                             'Tag', 'TMW_COLORBAR');
    if ~isempty(colorbarHandle)
      hImage(ismember(hImage,colorbarHandle)) = [];
    end
    
    hAxes = ancestor(hImage,'Axes');
    if iscell(hAxes)
      hAxes = [hAxes{:}]';
    end
    hFig = ancestor(h,'Figure'); %use h b/c hImage and hAxes could be empty

    % an array of images must belong to the same figure.
    if iscell(hFig) 
      figs = [hFig{:}];
      if any(figs ~= figs(1))
        eid = sprintf('Images:%s:invalidImageHandleArray',mfilename);
        msg = 'HIMAGE must belong to the same figure.';
        error(eid,'%s',msg);
      else
        hFig = hFig{1};
      end
    end

  else
    eid = sprintf('Images:%s:wrongKindOfGraphicsHandle',mfilename);
    msg = 'The handle must be a valid figure, axes, image, uicontainer, ';
    msg2 = 'or uipanel handle.';
    error(eid,'%s%s',msg,msg2);
  end
else
  eid = sprintf('Images:%s:invalidGraphicsHandle',mfilename);
  msg = 'The handle must be a valid graphics handle.';
  error(eid,'%s%s',msg);
end
