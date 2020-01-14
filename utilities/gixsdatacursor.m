function output_txt =  gixsdatacursor(~,event_obj,hFig)
% ***********************************************
% Copyright (c) 2017 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
gdata = get(hFig,'UserData');
if isempty(gdata) || ~isvalid(gdata)
    set(hFig,'UserData',[]);
    output_txt = {'Image handle has been removed.';'Please delete all datatips to reset.'};
    hc = datacursormode(gcbf);
    set(hc,'UpdateFcn',[]);
    return;
end
himg = get(event_obj,'Target');
xdata = get(himg,'xdata');
ydata = get(himg,'ydata');
pos = get(event_obj,'Position');
%[~,n]=min(abs(xdata-pos(1))); n=n(1);
%[~,m]=min(abs(ydata-pos(2))); m=m(1);
n = round(gu_axes2pix(gdata.ImDim(1),xdata,pos(1)));
m = round(gu_axes2pix(gdata.ImDim(2),ydata,pos(2)));
% construct output text
output_txt = {['[x,y] = [',num2str(n),',',num2str(m),']']};
output_txt = [output_txt;['I_raw          = ',num2str(gdata.RawData(m,n))]];
output_txt = [output_txt;['I_masked    = ',num2str(gdata.MaskedData(m,n))]];
if ~isempty(gdata.SolidAngleCorrectedData)
    output_txt = [output_txt;['I_corrected = ',num2str(gdata.SolidAngleCorrectedData(m,n))]];
end
if isempty(gdata.QMap),return; end
output_txt = [output_txt;['q          = ',num2str(gdata.QMap(m,n))]];
output_txt = [output_txt;['phi        = ',num2str(gdata.PhiMap(m,n))]];
if isempty(gdata.QzMap), return; end
output_txt = [output_txt;['qz         = ',num2str(gdata.QzMap(m,n))]];
output_txt = [output_txt;['qx         = ',num2str(gdata.QxMap(m,n))]];
output_txt = [output_txt;['qy         = ',num2str(gdata.QyMap(m,n))]];
output_txt = [output_txt;['qr         = ',num2str(gdata.QrMap(m,n))]];
output_txt = [output_txt;['2Theta = ',num2str(gdata.TwoThetaMap(m,n))]];
output_txt = [output_txt;['Alphaf = ',num2str(gdata.AlphafMap(m,n))]];
output_txt = [output_txt;['Chi    = ',num2str(gdata.ChiMap(m,n))]];