function scanresult = roiscan(filepath,filename,scanvarname,scanvarlist,params,const,dataflag,plotflag,cmaskdisplayflag)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% ROISCAN Extract integrated intensity over a user-defined ROI (or virtual 
%   detector slit) over a set of images.
%   SCANRESULT=ROISCAN(FILEPATH,FILENAME,SCANVARNAME,SCANVARLIST,PARAMS,
%   CONST,DATAFLAG,PLOTFLAG,CMASKDISPLAYFLAG) extracts ROI intensity for 
%   a set of images. It calls roicount method in gixsdata. 
%   
%   FILEPATH (string) and FILENAME are the path and list of names (in a nx1
%   cell array). 
%   
%   SCANVARNAME is the variable name during the scan (e.g., 'IncidentAngle'),
%   which has to math the property name of gixsdata if it is actively
%   assigned with values during the ROI scan.
%
%   SCANVARLIST is the corresponding list (nx1) of the scan variable (e.g, 
%   incident angle values), the total number of  which has to match that of
%   FILENAME on a one-to-one basis.
%
%   PARAMS is the gixsdata parameter. 
%
%   CONST is the constraint to define the ROI, and it is an nx1 cell array
%   following the usage of linecut method.
%
%   DATAFLAG=1/2 specifies either MaskedData or SolidAngleCorrectedData to
%   be used.
% 
%   PLOTFLAG=0/1 specifies whether the scanned image and reduced scan data 
%   are displayed and updated during the ROI scan data reduction.
%
%   CMASKDISPLAY=0/1 specifies whether a box inclosing ROI is drawn.
% 
%   Output SCANRESULT is a structure consisting of fields:
%       X and I: lists (nx1) of the scan dependent variable SCANVARLIST, 
%           and ROI intensity, respectively. 
%       MAPNAME: list of map names for the experiment geometry. It is 4x1 
%           for transmission or 11x1 for reflection. 
%       ROIMAP: list (nx4 or nx11) of averaged map values in the ROI for 
%           each image.
%       PCOUNTS: list (nx1) of the number of pixels in the ROI for each image.


%   Zhang Jiang @8ID/APS/ANL
%   $Revision: 1.0 $  $Date: 2013/01/24 $

n = length(filename);
I = NaN*ones(n,1);
pcounts = NaN*ones(n,1);
if params.Geometry == 1
    roimap = NaN*ones(n,4);
elseif params.Geometry == 2
    roimap = NaN*ones(n,11);
end
obj = params;
% if plotflag == 1
%     hfig = figure;
% end
for ii=1:n
    obj.ImFile = fullfile(filepath,filename{ii});
    if isprop(obj,scanvarname)
        obj.(scanvarname) = scanvarlist(ii);
        qmaps(obj);
    elseif ~isprop(obj,scanvarname) && ii==1
        qmaps(obj);
    end
    [I(ii),roimap(ii,:),mapname,pcounts(ii),cmask] = roicount(obj,const{ii},dataflag);
    if plotflag == 1
        % plot image
        if ii == 1
            hfig = figure('tag','roiscan_plottedfigures');
            hfig_pos = get(hfig,'position');
            set(hfig,'position',[hfig_pos(1)-hfig_pos(3)/2,hfig_pos(2:4)]);
            imagesc(obj);
            himg = findall(obj.FigHandle,'tag',['gixsdata:img:',obj.ImFileName]);
            himg_axes = get(himg,'Parent');
            hfig_scan = figure('position',[hfig_pos(1)+hfig_pos(3)/2+15,hfig_pos(2:4)]);
            haxes_scan = axes('parent',hfig_scan);
            set(haxes_scan,'xlim',[0,n+1]);
            xlabel(haxes_scan,'image #');
            ylabel(haxes_scan,'Intensity');   
            box(haxes_scan,'on');
        else
            set(himg,'tag',['gixsdata:img:',obj.ImFileName]);
            obj.post_imagesc(obj);
%             % update axis limits so the roi is in the middle of the plot
%             xlims = get(himg_axes,'xlim');
%             ylims = get(himg_axes,'ylim');         
%             xlims_new = xlims+ (roimap(ii,10)-mean(xlims));
%             ylims_new = ylims+ (roimap(ii,11)-mean(ylims));
%             if xlims_new(2)>obj.ImDim(1)
%                 xlims_new = [xlims(1),obj.ImDim(1)+0.5];
%             end
%             if xlims_new(1)<0
%                 xlims_new = [0.5,xlims(2)];                
%             end
%             if ylims_new(2)>obj.ImDim(2)
%                 ylims_new = [ylims(1),obj.ImDim(2)+0.5];
%             end
%             if ylims_new(1)<0
%                 ylims_new = [0.5,ylims(2)];                
%             end            
%             set(himg_axes,'xlim',xlims_new);
%             set(himg_axes,'ylim',ylims_new);            
        end
        if cmaskdisplayflag == 1
              displaycmask(cmask,himg);
        end
        cla(haxes_scan);
        %plot(1:n,I,'bo-','parent',haxes_scan);
        line(1:n,I,'marker','o','parent',haxes_scan);
        pause(0.05);
    end
end
scanresult.X = scanvarlist;
scanresult.I = I;
scanresult.roimap = roimap;
scanresult.mapname = mapname;
scanresult.pcounts = pcounts;

%%  --- add a slider to figure 
if plotflag == 1
    set(hfig,'unit','pixel',...
    'MenuBar','figure',...
    'Toolbar','auto');    
    figPos = get(hfig,'position');
    figSize = figPos(3:4);
    uicontrol('Parent',hfig,...
    'style','slider',...
    'Min',1,'Max',n,'Value',n,...
    'SliderStep',[1/n,10/n],...
    'unit','pixel',...
    'Position',[figSize(1)-35,figSize(2)/8,20,figSize(2)*6/8],...
    'Tag','roiscan_Slider',...
    'callback',{@roiscan_SliderFcn,filepath,filename,scanvarname,scanvarlist,obj,const,cmaskdisplayflag,himg});
     set(hfig,...
         'Toolbar','figure',...
         'Unit','pixel',...
         'ResizeFcn',@roiscan_ImageFigureResizeFcn);
end

function roiscan_ImageFigureResizeFcn(~,~)
hfig = gcbf;
hSlider = findall(hfig,'tag','roiscan_Slider');
figPos = get(hfig,'Position');
figSize = figPos(3:4);
set(hSlider,'position',[figSize(1)-35,figSize(2)/8,20,figSize(2)*6/8]);

function roiscan_SliderFcn(~,~,filepath,filename,scanvarname,scanvarlist,obj,const,cmaskdisplayflag,himg)
hSlider = gcbo;
ii = round(get(hSlider,'Value'));
set(hSlider,'Value',ii);
obj.ImFile = fullfile(filepath,filename{ii});
if isprop(obj,scanvarname)
    obj.(scanvarname) = scanvarlist(ii);
    qmaps(obj);
end
set(himg,'tag',['gixsdata:img:',obj.ImFileName]);
obj.post_imagesc(obj);
if cmaskdisplayflag == 1
    cmask = obj.get_cmask(obj,const{ii});
    displaycmask(cmask,himg);  
end

function displaycmask(cmask,himg)
v=ver;
checkimagetoolbox = any(strcmp('Image Processing Toolbox', {v.Name})); %1/0: yes/no
if checkimagetoolbox
    B = bwboundaries(cmask);
    delete(findall(get(himg,'Parent'),'tag','roiscan_LineROI'));
    for iLine = 1:length(B)
        line(B{iLine}(:,2),B{iLine}(:,1),'Parent',get(himg,'Parent'),'linestyle','-','color','r','tag','roiscan_LineROI');
    end
else
    cdata = get(himg,'CData');
    cdata(~cmask) = NaN;
    set(himg,'CData',cdata);
end