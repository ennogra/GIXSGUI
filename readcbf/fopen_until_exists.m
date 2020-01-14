%
% Filename: $RCSfile: fopen_until_exists.m,v $
%
% $Revision: 1.8 $  $Date: 2010/11/12 15:14:21 $
% $Author: bunk $
% $Tag: $
%
% Description:
% Open a file, in case of failure retry repeatedly if this has been
% specified. 
%
% Note:
% Call without arguments for a brief help text.
%
% Dependencies:
% none
%
%
% history:
%
% September 5th 2009, Oliver Bunk:
% bug fix in the zero file length check
%
% August 28th 2008, Oliver Bunk: 
% use dir rather than fopen to check for the file and check additionally
% that it is not of length zero
%
% May 9th 2008, Oliver Bunk: 1st version
%
function [fid,vararg_remain] = fopen_until_exists(filename,varargin)

if (nargin < 1)
    fprintf('Usage:\n');
    fprintf('[fid] = %s(filename [[,<name>,<value>],...]);\n',...
        mfilename);
    fprintf('filename                             name of the file to open\n');
    fprintf('The optional name value pairs are:\n');
    fprintf('''RetryReadSleep'',<seconds>           if greater than zero retry opening after this time (default: 0.0)\n');
    fprintf('''RetryReadMax'',<0-...>               maximum no. of retries, 0 for infinity (default: 0)\n');
    fprintf('''MessageIfNotFound'',<0-no,1-yes>     display a mesage if not found, 1-yes is default\n');
    fprintf('''ErrorIfNotFound'',<0-no,1-yes>       exit with an error if not found, default is 1-yes\n');
    fprintf('The file ID of the opened file is returned or -1 in case of failure.\n');
    error('Invalid number of input parameters.');
end

% check minimum number of input arguments
if (nargin < 1)
    display_help();
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
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end

% set default values for the variable input arguments and parse the named
% parameters:

% If the file has not been found and if this value is greater than 0.0 than
% sleep for the specified time in seconds and retry reading the file. 
% This is repeated until the file has been successfully read
% (retry_read_max=0) or until the maximum number of iterations is exceeded
% (retry_read_max>0). 
retry_read_sleep_sec = 0.0;
retry_read_max = 0;

% exit with error message if the file has not been found
error_if_not_found = 1;

% display a message once in case opening failed
message_if_not_found = 1;

% parse the variable input arguments
vararg_remain = cell(0,0);
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'RetryReadSleep'
            retry_read_sleep_sec = value;
        case 'RetryReadMax'
            retry_read_max = value;
        case 'MessageIfNotFound'
            message_if_not_found = value;
        case 'ErrorIfNotFound'
            error_if_not_found = value;
        otherwise
            vararg_remain{end+1} = name;
            vararg_remain{end+1} = value;            
    end
end


% try to access the file entry
file_non_empty = 0;
dir_entry = dir(filename);

% if it has not been found or if it is empty
if ((isempty(dir_entry)) || (size(dir_entry,1) == 0) || ...
    (dir_entry.bytes <= 0))
    if (message_if_not_found)
        if (isempty(dir_entry))
            fprintf('%s not found',filename);
        else
            fprintf('%s found but of zero length',filename);            
        end
    end
    % retry, if this has been specified
    if (retry_read_sleep_sec > 0.0)
        if (message_if_not_found)
            fprintf(', retrying\n');
        end
        % repeat until found or the specified number of repeats has been
        % exceeded (zero repeats means repeat endlessly)
        retry_read_ct = retry_read_max;
        while ((~file_non_empty) && ...
               ((retry_read_max <= 0) || (retry_read_ct > 0)))
            pause(retry_read_sleep_sec);
            dir_entry = dir(filename);
            if ((~isempty(dir_entry)) && (dir_entry.bytes > 0))
                file_non_empty = 1;
                % !!! arbitrary pause to increase the likelihood that an
                % HDF5 file is written !!!
                pause(1);
            end
            retry_read_ct = retry_read_ct -1;
        end
    else
        fprintf('\n');
    end
else
    file_non_empty = 1;
end

% open the file for read access
if (file_non_empty)
    fid = fopen(filename,'r');
else
  fid = -1;
end

% exit with an error message, if this has been specified and if the file
% could not be opened 
if (fid < 0)
    if (error_if_not_found)
        error('file ''%s'' not found',filename);
    end
end
