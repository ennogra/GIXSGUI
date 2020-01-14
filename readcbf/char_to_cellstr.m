%
% Filename: $RCSfile: char_to_cellstr.m,v $
%
% $Revision: 1.3 $  $Date: 2008/06/23 13:49:50 $
% $Author: bunk $
% $Tag: $
%
% Description:
% Convert an array of text to a cell array of lines.
%
% Note:
% Used for making file headers accessible. 
%
% Dependencies:
% none
%
%
% history:
%
% June 22nd 2008, Oliver Bunk: bug fix for fliread adding the nl_only 
% parameter, to be replaced by named parameter later on
%
% May 9th 2008, Oliver Bunk: 1st version
%
function [outstr] = char_to_cellstr(inchars,nl_only)

if (nargin < 2)
    nl_only = 0;
end

% get positions of end-of-line signatures
eol_ind = regexp(inchars,'\r\n');
eol_offs = 1;
if ((length(eol_ind) < 1) || (nl_only))
    eol_ind = regexp(inchars,'\n');
    eol_offs = 0;
end
if (length(eol_ind) < 1)
    eol_ind = length(inchars) +1;
end
if (length(eol_ind) < 1)
    outstr = [];
    return;
end

% dimension return array with number of lines
outstr = cell(length(eol_ind),1);

% copy the lines to the return array, suppressing empty lines
start_pos = 1;
ind_out = 1;
for (ind = 1:length(eol_ind))
    end_pos = eol_ind(ind) -1;
    % cut off trailing spaces
    while ((end_pos >= start_pos) && (inchars(end_pos) == ' '))
        end_pos = end_pos -1;
    end
    % store non-empty strings
    if (end_pos >= start_pos)
        outstr{ind_out} = inchars(start_pos:end_pos);
        ind_out = ind_out +1;
    end
    
    start_pos = eol_ind(ind) +1 + eol_offs;
    ind = ind +1;
end

% resize cell array in case of empty lines
if (ind_out <= length(eol_ind))
    outstr = outstr(1:(ind_out-1));
end

