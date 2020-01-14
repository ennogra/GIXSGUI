% COPY puts a string, cell array, or numerical array into the clipboard,
% from where it can easily be imported by other programs.
%
% Usage:
%   copy(x)
%
% x can be a numerical array, a cell array, or a string. The program does
% not handle other variable classes, such as struct etc.
%
% Note: This program was inspired by NUM2CLIP on the Mathworks file
% exchange, see http://www.mathworks.com/matlabcentral/fileexchange/8472
%
% Author : Yvan Lengwiler
% Release: $1.3$
% Date   : $2011-06-22$
%
% See also PASTE, CLIPBOARD

% History:
% 2010-06-25    First version.
% 2010-07-10    Now also covers multiline character arrays.
% 2011-06-05    Using strrep instead of regexprep.
% 2011-06-22    Special treatment of empty cells (Thank you, Joseph Burgel).

function copy(x)
    % *** string argument ************************************************
    if isa(x,'char')
        [r,c] = size(x);
        if r == 1                   % not a multi-line character array ...
            clipboard('copy',x);    % ... so just push the string into the
                                    %     clipboard
        else
            x = [x, repmat(char(10),r,1)];  % append linefeed to each line
            x = reshape(x',1,r*(c+1));      % make it a single line
            clipboard('copy',x);            % push this to the clipboard
        end
    % *** numeric argument ***********************************************
    elseif isa(x,'numeric')
        s = mat2str(x);                 % write as [.. .. ..;.. .. ..]
        if s(1)=='['                    % it's a proper array
            s = s(2:end-1);             % remove '[' and ']'
        end
        s = strrep(s,' ',char(9));      % replace spaces with tabs
        s = strrep(s,';',char(10));     % replace semicolons with linefeeds
        s(end+1) = char(10);            % append a linefeed
        clipboard('copy',s);            % place resulting string in clipboard
    % *** cell argument **************************************************
    elseif isa(x,'cell')
        [nrow, ncol] = size(x);
        str = '';
        for r = 1:nrow
            for c = 1:ncol-1
                str = onecell(str, x{r,c}, char(9));    % treat cell, append a tab
            end
            str = onecell(str, x{r,end}, char(10));     % treat cell, append a linefeed
        end
        clipboard('copy',str);          % copy to clipboard
    % *** anything else **************************************************
    else
        warning('COPY:unsupported_content', ...
            'cannot copy this kind of object.');
    end

% *** convert one cell into a string; append it to str and append a special
%     character (ch = tab or linefeed)
function str = onecell(str,e,ch)
    if isempty(e)
        str = [str, ch];            % copy nothing if cell is empty
    elseif isa(e,'char')
        if size(e,1) == 1           % not a multi-line char array?
            str = [str, e, ch];
        else
            str = [str, mat2str(e), ch];
        end
    elseif isa(e,'numeric')
        str = [str, mat2str(e), ch];
    else
        str = [str, '(cannot copy component)', ch];
    end
