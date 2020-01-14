%
% Filename: $RCSfile: cbfread.m,v $
%
% $Revision: 1.1 $  $Date: 2008/06/10 17:05:13 $
% $Author: bunk $
% $Tag: $
%
% Description:
% Macro for reading Crystallographic Binary File (CBF) files written by the
% Pilatus detector control program camserver. 
%
% Note:
% Compile the C-program cbf_uncompress using mex (see header of
% cbf_uncompress.c) to use it for uncompression instead of the slower
% Matlab code. 
% Currently this routine supports only the subset of CBF features needed to
% read the Pilatus detector data. 
% Call without arguments for a brief help text.
%
% Dependencies:
% - image_read_set_default
% - fopen_until_exists
% - get_hdr_val
% - compiling cbf_uncompress.c increases speed but is not mandatory
%
%
% history:
%
% May 9th 2008, Oliver Bunk: 1st version
%

function [frame,vararg_remain] = cbfread(filename,varargin)

% 0: no debug information
% 1: some feedback
% 2: a lot of information
debug_level = 0;

% initialize return argument
frame = struct('header',[], 'data',[]);


% check minimum number of input arguments
if (nargin < 1)
    image_read_sub_help(mfilename,'cbf');
    error('At least the filename has to be specified as input parameter.');
end

% accept cell array with name/value pairs as well
no_of_in_arg = nargin;
if (nargin == 2)
    if (isempty(varargin))
        % ignore empty cell array
        no_of_in_arg = no_of_in_arg -1;
    else
        if (iscell(varargin{1}))
            % use a filled one given as first and only variable parameter
            varargin = varargin{1};
            no_of_in_arg = 1 + length(varargin);
        end
    end
end

% check number of input arguments
if (rem(no_of_in_arg,2) ~= 1)
    error('The optional parameters have to be specified as ''name'',value pairs');
end
    
% set default values for the variable input arguments and parse the named
% parameters: 
vararg = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        otherwise
            % pass further arguments on to fopen_until_exists
            vararg{end+1} = name;
            vararg{end+1} = value;
    end
end


% expected maximum length for the text header
max_header_length = 4096;

% end of header signature
eoh_signature = char([ 12 26 4 213 ]);

% CBF file signature
cbf_signature = '###CBF: VERSION';

% Calling an external C routine for uncompressing the data did save about
% 30% time on a specific machine. 
% The C-routine is used if a compiled version of it exists. 
% See the header of cbf_uncompress.c for information on how to compile the
% C file using mex in Matlab. 
c_routine = (length(which('cbf_uncompress')) > 0);

% try to open the data file
if (debug_level >= 1)
    fprintf('Opening %s.\n',filename);
end
[fid,vararg_remain] = fopen_until_exists(filename,vararg);
if (fid < 0)
    return;
end

% read all data at once
[fdat,fcount] = fread(fid,'uint8=>uint8');

% close input data file
fclose(fid);
if (debug_level >= 2)
    fprintf('%d data bytes read\n',fcount);
end

% search for end of header signature within the expected maximum length of
% a header
end_of_header_pos = ...
    strfind( fdat(1:min(max_header_length,length(fdat)))',...
             eoh_signature );
if (length(end_of_header_pos) < 1)
    cbf_error(filename,'no header end signature found');
    return;
end
if (debug_level >= 2)
    fprintf('Header length is %d bytes.\n',end_of_header_pos -1);
end

% return the complete header as lines of a cell array
frame.header = char_to_cellstr( char(fdat(1:(end_of_header_pos-1))') );

% check for CBF signature
if (~strncmp(cbf_signature,frame.header{1},length(cbf_signature)))
    cbf_error(filename,[ 'CBF signature ''' cbf_signature ...
        ''' not found in first line ''' frame.header{1} '''' ]);
end

% extract the mandatory information for decompression from the header
no_of_bin_bytes = get_hdr_val(frame.header,'X-Binary-Size:','%f',1);
dim1 = get_hdr_val(frame.header,'X-Binary-Size-Fastest-Dimension:','%f',1);
dim2 = get_hdr_val(frame.header,'X-Binary-Size-Second-Dimension:','%f',1);
el_type = get_hdr_val(frame.header,'X-Binary-Element-Type: "','%[^"]',1);
switch (el_type)
    case 'signed 32-bit integer'
        bytes_per_pixel = 4;
    otherwise
        cbf_error(filename,[ 'unknown element type ' el_type ]);
end
compr_type = get_hdr_val(frame.header,'conversions="','%[^"]',1);
switch (compr_type)
    case 'x-CBF_BYTE_OFFSET' 
        compression_type = 1;
    otherwise
        cbf_error(filename,[ 'unknown compression type ' compr_type ]);
end
if (debug_level >= 2)
    fprintf('Frame dimensions are %d x %d.\n',dim2,dim1);
end

% uncompress the binary data
[frame.data] = ...
    extract_frame(fdat((end_of_header_pos+length(eoh_signature)):end),...
        dim1,dim2,no_of_bin_bytes,compression_type,...
        filename,...
        c_routine,debug_level);
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = cbf_error(filename,text)

fprintf('cbfread of %s:\n %s\n',filename,text);
return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [frame] = ...
    extract_frame(dat_in,...
        dim1,dim2,no_of_in_bytes,compression_type,...
        filename,...
        c_routine,debug_level)

if (c_routine)    
    if (debug_level >= 2)
        fprintf('C routine called.\n');
    end
    [frame] = ...
        cbf_uncompress(dat_in,dim1,dim2,no_of_in_bytes,compression_type);
    return;
end

if (debug_level >= 2)
    fprintf('Matlab routine called.\n');
end


% initialize return array
frame = zeros(dim1,dim2);

% only byte-offset compression is supported
if (compression_type ~= 1)
    cbf_error(filename,...
        ['extract_frame does not support compression type no. ' ...
        num2str(compression_type)]);
end


% In byte-offset compression the difference to the previous pixel value is
% stored as a byte, 16-bit integer or 32-bit integer, depending on its
% size. 
% The sizes above one byte are indicated by the escape sequence -1 in the
% previous data format, i.e, a 32-bit integer is preceded by the sequence %
% 0x80 (too large for a byte)
% 0x8000 (too large for a 16-bit integer). 
ind_out = 1;
ind_in = 1;
val_curr = 0;
val_diff = 0;
while (ind_in <= no_of_in_bytes)
    val_diff = double(dat_in(ind_in));
    ind_in = ind_in +1;
    if (val_diff ~= 128)
        % if not escaped as -128 (0x80=128) use the current byte as
        % difference, with manual complement to emulate the sign
        if (val_diff >= 128)
            val_diff = val_diff - 256;
        end
    else
        % otherwise check for 16-bit integer value
        if ((dat_in(ind_in) ~= 0) || (dat_in(ind_in+1) ~= 128))
            % if not escaped as -32768 (0x8000) use the current 16-bit integer
            % as difference 
            val_diff = double(dat_in(ind_in)) + ...
                256 * double(dat_in(ind_in+1));
            % manual complement to emulate the sign
            if (val_diff >= 32768)
                val_diff = val_diff - 65536;
            end
            ind_in = ind_in +2;
        else
            ind_in = ind_in +2;
            % if everything else failed use the current 32-bit value as
            % difference
            val_diff = double(dat_in(ind_in)) + ...
                256 * double(dat_in(ind_in+1)) + ...
                65536 * double(dat_in(ind_in+2)) + ...
                16777216 * double(dat_in(ind_in+3));
            % manual complement to emulate the sign
            if (val_diff >= 2147483648)
                val_diff = val_diff - 4294967296;
            end
            ind_in = ind_in +4;
        end
    end
	val_curr = val_curr + val_diff;
    frame(ind_out) = val_curr;
    ind_out = ind_out +1;
end
    
if (ind_out-1 ~= dim1*dim2)
    cbf_error(filename,[ 'mismatch between ' num2str(ind_out-1) ...
        ' bytes after decompression with ' num2str(dim1*dim2) ...
        ' expected ones' ]);
end


% if (~original_orientation)
%     frame.data = frame.data(end:-1:1,end:-1:1)';
% end
