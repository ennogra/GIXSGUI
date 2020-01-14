function sresult = rdspec(varargin)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% RDSPEC Non-interactive function to load a single scan from spec file.
%   S = RDSPEC(FILENAME,SCAN)
%
%   Input Argument:
%       FILENAME : Spec file name
%           SCAN : Single scan number
%
%   Output Argument:
%       S.SCAN  : Scan title
%       S.TIME : Scan time
%       S.COLS : Total number of the columns
%       S.HEAD : Column names
%       S.DATA : Scan data
%
% Copyright 2006 Zhang Jiang
% $Revision: 1.0 $  $Date: 2006/01/24 $

if nargin ~=2 
    error('Invalid input argument.');
    return;
end
file    = varargin{1};
scan    = varargin{2};
scanStr = ['#S ',num2str(scan),' '];
[fid,message] = fopen(file);        % open file
if fid == -1                % return if open fails
    uiwait(msgbox(message,'File Open Error','error','modal'));
    fclose(fid);
    restorePath(prePath);
    return;
end

% --- find scan
scanFlag = 0;
while feof(fid) == 0
    scanline    = fgetl(fid);
    matches = findstr(scanline,scanStr);
    if length(matches) > 0
        scanFlag = 1;
        break;
    end
end

% --- return if scan does not exsit
if scanFlag == 0
    error('Scan does not exist.');
    fclose(fid);
    return;
end

% --- get scan information
s.scan = scanline;
while ~strcmp(scanline(1:2),'#D')
    scanline = fgetl(fid);
end
s.time = scanline;
while ~strcmp(scanline(1:2),'#N')
    scanline = fgetl(fid);
end
s.cols = str2double(scanline(4:end));
while ~strcmp(scanline(1:2),'#L')
    scanline = fgetl(fid);
end
scanline = scanline(4:end);
space = findstr(scanline,'  ');
lengthSpace = length(space);
for iSpace = lengthSpace:-1:2
    if space(iSpace) == space(iSpace-1)+1
        space(iSpace) = [];
    end
end
space = [-1 space length(scanline)+1];
s.head = cell(1,s.cols);
for iCol = 1:s.cols
    colHead = scanline(space(iCol)+2:space(iCol+1)-1);
    while colHead(1) == ' '
        colHead(1) = '';
    end
    s.head{iCol} = colHead;
end

% --- read scan data
colData = [];
str_scanline = num2str(scanline);
while ~strcmp(str_scanline,'') & ~strcmp(str_scanline,'-1')
    fidPos = ftell(fid);
    scanline = fgetl(fid);
    str_scanline = num2str(scanline);
    while ~strcmp(str_scanline,'') & ~strcmp(str_scanline,'-1')...
            & ~strcmp(str_scanline(1:1),'#')
        fseek(fid,fidPos,'bof');
        colData = [colData fscanf(fid,'%g',[s.cols,1])];
        scanline = fgetl(fid);
        fidPos = ftell(fid);
        scanline = fgetl(fid);
        str_scanline = num2str(scanline);
    end
end
s.data = colData';

fclose(fid);
sresult = s;