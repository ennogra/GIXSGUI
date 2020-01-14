function sddcalibgui
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% SDDCALIBGUI Calibrate sample-to-detector distance using reflection.

%   Zhang Jiang @8ID/APS/ANL
%   $Revision: 1.0 $  $Date: 2014/10/28 $
%   $Revision: 1.1 $  $Date: 2015/04/23 $
%       (1) Make the 3rd numeric also when creating the table

hFigGIXS = findall(0,'Tag','gixs_fig');
if isempty(hFigGIXS)
    warndlg('Active GIXSGUI window is required.','SDD Calibration Warning','modal');
    return;
end

hFigSDDCalib = findall(0,'Tag','sddcalib_fig');
if ~isempty(hFigSDDCalib)
    figure(hFigSDDCalib);
    return;
end

%% initialize udata field
udata.tdata = [];      % table data: column 1/2/3 for th, x and y pixels
udata.SDD = NaN;    % sdd in mm

%% figure layout
screenSize = get(0,'screensize');
figureSize = [440,620];
figPos = [(screenSize(3)-figureSize(1))/2, (screenSize(4)-figureSize(2))/2, figureSize];
hFigSDDCalib = figure(...
    'DockControls','off',...
    'Resize','off',...
    'position',figPos,...
    'PaperOrient','portrait',...
    'PaperPositionMode','auto',...
    'IntegerHandle','off',...
    'NumberTitle','off',...
    'MenuBar','none',...
    'HandleVisibility','on',...
    'Toolbar','none',...
    'Name','SDD Calibration - 1.0 (8ID/APS/ANL)',...
    'Tag','sddcalib_fig','UserData',udata);
backgroundcolor = get(hFigSDDCalib,'color');
panelcolor = backgroundcolor;
hPanelSDDCalib = uipanel('Parent',hFigSDDCalib,...
    'BackgroundColor',panelcolor,...
    'units','norm',...
    'Position',[0.005 0.045 0.99 0.955],...
    'Title','SDD Calibration Using Speular Reflections',...
    'TitlePosition','centertop',...
    'tag','sddcalib_PanelROI');
SizePanelSDDCalib = get(hPanelSDDCalib,'position'); 
SizePanelSDDCalib = SizePanelSDDCalib(3:4).*figureSize;

% instruction
PosInstr = [5,SizePanelSDDCalib(2)-120,200,15];
str_instr = [...
    'Sample-to-detector distance calibration can be done using the reflections ',...
    'from a thin film or flat substrate (usually silicon wafer). Matlab ',...
    'Optimization Toolbox is required. The calibration needs to load direct ',...
    'beam position [Beam0] and pixel sizes from main GIXSGUI window. You need to ',...
    'enter below the pixel positions of the specular spots for at least two incident angles. ',...
    '"Find COM" can help find the center of mass of the specular spots.'...
    ];
uicontrol('Parent',hPanelSDDCalib,...
    'style','text',...
    'Units','pixel',...
    'backgroundcolor',panelcolor,...
    'String',str_instr,...
    'HorizontalAlignment','left',...
    'Position',[5,PosInstr(2),SizePanelSDDCalib(1)-10,100]);

% table
PosList = [5,SizePanelSDDCalib(2)-155,100,15];
uicontrol('Parent',hPanelSDDCalib,...
    'style','Text',...
    'Units','pixel',...
    'backgroundcolor',panelcolor,...
    'String','Enter at least two incident angles and their specular pixel positions: ',...
    'HorizontalAlignment','left',...
    'Position',[5,PosList(2),300,15]);
colName = {'Theta (deg)','Pixel X', 'Pixel Y'};
nRows = 10;
rowName = {1:nRows};
tdata = cell(nRows,3);
tdata(:,:) = {[]};
hTableList = uitable('Parent',hPanelSDDCalib ,...
    'unit','pixel',...
    'Position',[1,PosList(2)-200,SizePanelSDDCalib(1)-3,195],...
    'Tag','sddcalib_TableList',...
    'UserData',[]);         % userdata for selected cell indices
set(hTableList, 'ColumnName', colName,...
    'RowName',rowName,...
    'columnformat',repmat({'numeric'},[1,2,3]),...
    'ColumnWidth', repmat({133},[1,3]),...
    'ColumnEditable', true(1,length(colName)),...
    'Data',tdata);
% 'CellEditCallBack',@sddcalib_TableList_CellEditCallback,...
   % 'CellSelectionCallback',@sddcalib_TableList_CellSelectionCallback');
set(hFigSDDCalib,'HandleVisibility','callback');         % after adding uimenu, set main figure handlvisibility backto callback

PosCalib = [5,SizePanelSDDCalib(2)-380,100,20];
uicontrol('Parent',hPanelSDDCalib,...
    'style','pushbutton',...
    'String','Clear',...
    'unit','pixel',...
    'Position',[SizePanelSDDCalib(1)-105*4,PosCalib(2)-0,100,20],...
    'Tag','sddcalib_PushbuttonCalibrate',...
    'ToolTipString','Clear table',...   
    'callback',@sddcalib_PushbuttonClearFcn);
uicontrol('Parent',hPanelSDDCalib,...
    'style','pushbutton',...
    'String','Calibrate',...
    'unit','pixel',...
    'Position',[SizePanelSDDCalib(1)-105*1,PosCalib(2)-0,100,20],...
    'Tag','sddcalib_PushbuttonCalibrate',...
    'callback',@sddcalib_PushbuttonCalibrateFcn);

% --- Quit/Hide Figure
uicontrol('Parent',hFigSDDCalib,...
    'style','pushbutton',...
    'String','Quit',...
    'unit','pixel',...
    'Position',[5,5,70,20],...
    'Tag','sddcalib_PushbuttonQuit',...
    'callback',@sddcalib_CloseRequestFcn);
uicontrol('Parent',hFigSDDCalib,...
    'style','pushbutton',...
    'String','Close All',...
    'unit','pixel',...
    'ToolTipString','Close all figures',...
    'Position',[80,5,70,20],...
    'Tag','sddcalib_PushbuttonCloseAll',...
    'callback',@sddcalib_PushbuttonCloseAllFcn);

function sddcalib_PushbuttonCalibrateFcn(~,~)
hFigSDDCalib = gcbf;
udata = get(hFigSDDCalib,'UserData');
hTableList = findall(hFigSDDCalib,'tag','sddcalib_TableList');
% get table data
tdata = get(hTableList,'Data');
udata.tdata = []; 
for ii=1:size(tdata,1)
    tmp_data = cell2mat(tdata(ii,:));
    if length(tmp_data)==3 && nnz(isnan(tmp_data))==0
        udata.tdata = [udata.tdata;tmp_data];
    end
end
if isempty(udata.tdata) || size(udata.tdata,1)<2
    return;
end
% get GIXSGUI data
hFigGIXS = findall(0,'Tag','gixs_fig');
if isempty(hFigGIXS)
    warndlg('Active GIXSGUI window is required.','SDD Calibration Warning','modal');
    return;
end
udata_gixsgui = get(hFigGIXS,'UserData');
udata.Beam0 = udata_gixsgui.paramsDefault.Beam0;
udata.PixelSize = udata_gixsgui.paramsDefault.PixelSize;
if nnz(isnan(udata.Beam0)) ~=0 || nnz(isnan(udata.PixelSize)) ~=0
   warndlg('Invalid pixel size or beam0 position.','SDD Calibration Warning','modal');     
   return;
end
% calibrate
udata.tdata = sortrows(udata.tdata,1);
th = udata.tdata(:,1);
specx = udata.tdata(:,2);
specy = udata.tdata(:,3);
beam0x = udata.Beam0(1);
beam0y = udata.Beam0(2);
xpix = udata.PixelSize(1);
ypix = udata.PixelSize(2);
r = sqrt((specx-beam0x).^2*xpix^2+(specy-beam0y).^2*ypix^2);
fresult = sdd_calibration(th,r);
udata.SDD = fresult.sdd;
set(hFigSDDCalib,'UserData',udata);

function fresult = sdd_calibration(th,r)
% SDD_CALIBRATION(TH,R) Calibrate sample-to-detector distance using reflection.
%   FRESULT = SDD_CALIBRATION(TH,R) calculate SDD and incident angle
%   offset. TH is the nominal incident angle (unit: deg). R is the distance
%   (unit: mm) from the specular relfection to beam0.

%   Zhang Jiang @8ID/APS/ANL
%   $Revision: 1.0 $  $Date: 2013/03/18 $


%% prepare data
th = th(:);
r = r(:);

%% initialize parameters
sdd = 200;
dth = 0;

%% construct fitting structure
x       = [sdd  dth];
lb      = [eps  -Inf];
ub      = [Inf  Inf];
fitFlag = [1     1];
x1 = x(fitFlag==1);
x2 = x(fitFlag==0);
lb1 = lb(fitFlag==1);
ub1 = ub(fitFlag==1);

%% start fitting
options = optimset (...
    'Display','off',...
     'TolX',1e-8,...
     'TolFun',1e-8,...
     'MaxFunEvals',1600,...
     'MaxIter',1200);
[fittedX1,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(@fcn_fit,x1,th,r,lb1,ub1,options,x2,fitFlag);
% --- error estimation
s2=resnorm/(length(residual) - 3);
[~,R]=qr(jacobian,0);
Rinv=inv(R);
sigmaest=(Rinv*Rinv')*s2;

stderrors=full(sqrt(diag(sigmaest)));

%% assign fitted result and plot
x=zeros(1,length(fitFlag));
x(fitFlag==1) = fittedX1;
x(fitFlag==0) = x2;
sdd   = x(1);
dth   = x(2);

r_cal = fcn_fit(fittedX1,th,x2,fitFlag);
th_ext = [0;min(abs(th))];
r_cal_ext = fcn_fit(fittedX1,th_ext,x2,fitFlag);
figure('tag','sddcalib_plottedfigures')
hold on;
plot(th,r,'bo');
plot(th,r_cal,'r-');
plot(th_ext,r_cal_ext,'r:');

hold off; box on;
legend('data','fit','location','northwest');
xlabel('\theta (deg)');
ylabel('Distance to Beam0 (mm)');
text_str = {...
    'Model: r=SDD*tan(\theta+\Delta\theta)';
    sprintf('%s%0.3f%s%0.3f%s','SDD = ',sdd,' \pm ',stderrors(1),' mm');
    sprintf('%s%0.4f%s%0.4f%s','\Delta\theta = ',dth,' \pm ',stderrors(2),' deg')};
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
posx = xlim(1) + (xlim(2)-xlim(1))/20;
posy = ylim(1) + (ylim(2)-ylim(1))/1.5;

text(posx,posy,text_str);    

%% output
fresult.sdd = sdd;
fresult.sdd_std = stderrors(1);
fresult.dth = dth;
fresult.dth_std = stderrors(2);
fresult.resnorm = resnorm;
fresult.residual = residual;
fresult.exitflag = exitflag;

function I = fcn_fit(x1,th,x2,fitFlag)
% --- Generate all the fitting parameters
x=zeros(1,length(fitFlag));
x(fitFlag==1) = x1;
x(fitFlag==0) = x2;

sdd   = x(1);
dth   = x(2);

I = sdd*tand(2*(th+dth));

function sddcalib_PushbuttonClearFcn(~,~)
hFigSDDCalib = gcbf;
udata = get(hFigSDDCalib,'UserData');
udata.SDD = NaN;
set(hFigSDDCalib,'UserData',udata);
hTableList = findall(hFigSDDCalib,'tag','sddcalib_TableList');
tdata = cell(10,3);
tdata(:,:) = {[]};
set(hTableList,'Data',tdata );

function sddcalib_CloseRequestFcn(~,~)
delete(gcbf);

function sddcalib_PushbuttonCloseAllFcn(~,~)
delete(findall(0,'tag','sddcalib_plottedfigures'));

