function imwrite2tif(varargin)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% IMWRITE2TIF Write image to TIF file with specified datatype.
%   IMWRITE2TIF(IMGDATA,HEADER,IMFILE,DATATYPE) exports IMGDATA with HEADER
%   to TIF file named IMFILE. HEADER is usally obtained by IMFINFO from 
%   original image file, and it can be left empty []. String DATATYPE 
%   specifies data type to export. Supported data types include logical
%   uint8, int8, uint16, int16, uint32, int32, uint64, int64, single and
%   double.
%
%   IMWRITE2TIF(IMGDATA,HEADER,IMFILE,DATATYPE,TAG NAME 1,TAG VALUE 1,TAG 
%   NAME 2,TAG VALUE 2,...) writes with specified Matlab supported TIF
%   tags. New tag values directly defined in the input override those 
%   already defined in HEADER. Type Tiff.getTagNames() to list all Matlab 
%   supported Tif tags.
%
%   Note 1: 
%       To avoid errors such as '??? Error using ==> tifflib The value for
%       MaxSampleValue must be ...', overide the tag of problem, e.g., 
%       MaxSampleValue by Matlab supported values. Or simply remove the 
%       tag from HEADER and the default value will be assgined automatically.
%
%   Note 2:
%       IMWRITE2TIF does not check existing image files before overwriting.
%       Be cautious with export image file names.
%
%   Example 1:
%       imgdata = imread('ngc6543a.jpg');
%       header  = imfinfo('ngc6543a.jpg');
%       imwrite2tif(imgdata,header,'new_ngc6543a.tif','uint8');
%
%   Example 2:
%       imgdata = imread('mri.tif');
%       imwrite2tif(imgdata,[],'new_mri.tif','int32','Copyright','MRI',
%       'Compression',1);
%
%   More information can be found by searching for 'Exporting Image Data 
%   and Metadata to TIFF Files' in Matlab Help.

%   Zhang Jiang 
%   $Revision: 1.0 $  $Date: 2011/02/23 $

if nargin<4 || mod(nargin,2)==1
    error('Invalid number of input arguments.');
end

% assign input argument
imgdata  = varargin{1};
header   = varargin{2};
imfile   = varargin{3};
datatype = varargin{4};
if nargin>4
    tag_name  = cell((nargin-4)/2,1);
    tag_value = cell((nargin-4)/2,1);
    for ii=1:(nargin-4)/2
        tag_name{ii}  = varargin{2*ii+3};
        tag_value{ii} = varargin{2*ii+4};
    end
end

% check errors
if ~isnumeric(imgdata) && ~islogical(imgdata)
     error('The first input argument (image data) must be a numeric array.');
end
if ~isempty(header) && ~isstruct(header)
    error('The second input argument (header info) must be empty or a structure.');
end
if ~ischar(imfile)
    error('The third input argument (output tif file name) must be string.');
end
switch lower(datatype)
    case 'logical'
        BitsPerSample = 1;
        SampleFormat  = 1;
        imgdata = logical(imgdata);
    case 'uint8'
        BitsPerSample = 8;   
        SampleFormat  = 1;         
        imgdata = uint8(imgdata);    
    case 'int8'
        BitsPerSample = 8;
        SampleFormat  = 2;       
        imgdata = int8(imgdata);          
    case 'uint16'
        BitsPerSample = 16;       
        SampleFormat  = 1;         
        imgdata = uint16(imgdata);           
    case 'int16'
        BitsPerSample = 16;       
        SampleFormat  = 2;         
        imgdata = int16(imgdata);        
    case 'uint32'
        BitsPerSample = 32;       
        SampleFormat  = 1;                
        imgdata = uint32(imgdata);            
    case 'int32'
        BitsPerSample = 32;       
        SampleFormat  = 2;     
        imgdata = int32(imgdata);         
     case 'single'
        BitsPerSample = 32;        
        SampleFormat  = 3;                
        imgdata = single(imgdata);       
    case {'uint64','int64','double'}
        BitsPerSample = 64;        
        SampleFormat  = 3;                
        imgdata = double(imgdata);           
    otherwise
        error('Invalid output data type.');
end

% creat a Tiff object
t = Tiff(imfile,'w');

% duplicate tags from header
tagnamelist = Tiff.getTagNames();
tagnamelist_delete = {...
    'StripByteCounts';...
    'StripOffsets';
    'TileByteCounts';...
    'TileOffsets';...
    'MaxSampleValue';...
    'MinSampleValue';...
    'ResolutionUnit'};
for ii=1:length(tagnamelist_delete)    % remove read only tag names
    tagnamelist(strcmpi(tagnamelist_delete{ii},tagnamelist)) = [];
end
if ~isempty(header)
    for ii=1:length(tagnamelist)
        if isfield(header,tagnamelist{ii}) && ~isempty(header.(tagnamelist{ii}))
           tagstruct.(tagnamelist{ii}) = header.(tagnamelist{ii});
        end
    end
end

% update tags determined from imgdata and datatype
tagstruct.ImageLength = size(imgdata,1);
tagstruct.ImageWidth = size(imgdata,2);
tagstruct.SamplesPerPixel = size(imgdata,3);
tagstruct.SampleFormat = SampleFormat;
tagstruct.BitsPerSample = BitsPerSample;

% update some default tags (these tags can be overriden by user input below)
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.Compression = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';

% update user input tag values
if nargin>4
    for ii=1:length(tag_name)
        tagstruct.(tag_name{ii}) = tag_value{ii};
    end
end

% set tags and write to tif file
t.setTag(tagstruct)
t.write(imgdata);
t.close();   