function varargout =supertitle(varargin)
%SUPTITLE Puts a title above all subplots.
%	SUPTITLE('text') adds text to the top of the figure above all subplots 
%   (a "super title"). Use this function after all subplot commands.
%
%	SUPTITLE(H,'text') adds title to figure with handle H.
%
%	SUPTITLE(H,'text','Property1',PropertyValue1,'Property2',PropertyValue2,...)
%   and sets the values of the specified properties of the super title.
%
%   HT = SUPERTITLE (...) returns handle HT for the title text.
%
%   [HT,HAX] = SUPERTITLE (...) also returns the handle HAX for the title axes.
%
% Drea Thomas 6/15/95 drea@mathworks.com
% John Cristion 12/13/00 modified
% Mark Histed 03/13/04 histed@mit.edu: fix disappearing legend on last plot
% Zhang Jiang 08/19/10 modified to accept various input argument
% 
% $Id: suptitle.m,v 1.2 2004/03/13 22:17:47 histed Exp $
%
% Warning: If the figure or axis units are non-default, this
% will break.

% --- determine the input argument
try
    if nargin == 1
        h_fig = gcf;
        str = varargin{1};
    elseif nargin == 2
        h_fig = varargin{1};
        str = varargin{2};
    elseif mod(nargin,2)==0
        h_fig = varargin{1};
        str = varargin{2};
        for i=1:(nargin-2)/2
            st_property{i}         = varargin{i*2+1};
            st_property_value{i}   = varargin{i*2+2};            
        end
    end
catch
    error('Incorrect input argument');
    return;
end

% --- Parameters used to position the supertitle.
plotregion = .99;                           % Amount of the figure window devoted to subplots
titleypos  = .965;                          % Y position of title in normalized coordinates
%fs = get(h_fig,'defaultaxesfontsize')+4;     % Fontsize for supertitle
fs = get(h_fig,'defaultaxesfontsize');
fudge=0.01;                                 % Fudge factor to adjust y spacing between subplots

% --- % Get the (approximate) difference between full height (plot + title
% + xlabel) and bounding rectangle.
set(0,'CurrentFigure',h_fig);
haold = get(h_fig,'CurrentAxes');
figunits = get(h_fig,'units');
if (~strcmp(figunits,'pixels')),
    set(h_fig,'units','pixels');
    pos = get(h_fig,'position');
    set(h_fig,'units',figunits);     
else
    pos = get(h_fig,'position');
end
ff = (fs-4)*1.27*5/pos(4)*fudge;
% The 5 here reflects about 3 characters of height below
% an axis and 2 above. 1.27 is pixels per point.

% --- Determine the bounding rectange for all the plots and if a title exists
h = findobj(h_fig,'Type','axes');  % Change suggested by Stacy J. Hills
max_y=0;
min_y=1;
oldtitle =0;
for i=1:length(h),
	if (~strcmp(get(h(i),'Tag'),'suptitle')),
		pos=get(h(i),'pos');
		if (pos(2) < min_y), min_y=pos(2)-ff/5*3;end;
		if (pos(4)+pos(2) > max_y), max_y=pos(4)+pos(2)+ff/5*2;end;
    else
		oldtitle = h(i);
	end
end

% --- if needed rescale the plotregion
if max_y > plotregion,
	scale = (plotregion-min_y)/(max_y-min_y);
	for i=1:length(h),
		pos = get(h(i),'position');
		pos(2) = (pos(2)-min_y)*scale+min_y;
		pos(4) = pos(4)*scale-(1-scale)*ff/5*3;
		set(h(i),'position',pos);
	end
end

np = get(h_fig,'nextplot');
set(h_fig,'nextplot','add');
if (oldtitle),
	delete(oldtitle);
end
ha=axes('pos',[0 1 1 1],'visible','off','Tag','suptitle');
ht=text(.5,titleypos-1,str);
set(ht,'horizontalalignment','center','fontsize',fs,'Interpreter','none');
if nargin>2
    for ii=1:(nargin-2)/2
        set(ht,st_property{ii},st_property_value{ii});
    end
end
set(h_fig,'nextplot',np);                   % reset the nextplot property and make the axes handles current
if ishandle(haold), set(h_fig,'CurrentAxes',haold); end        % switch back to the old gca

% --- fix legend if one exists
legH = legend;
if ~isempty(legH)
    axes(legH);
end

% --- output argument
if nargout == 1
	varargout{1} = ht;
elseif nargout == 2
	varargout{1} = ht;
	varargout{2} = ha;
end