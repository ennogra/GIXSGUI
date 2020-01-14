% PASTE copies the content of the clipboard to a variable; creates a
% string, a cell array, or a numerical array, depending on content.
%
% Usage:
%   x = paste;
%
% The program will try to create an array. For this to succeed, the
% material that is passed from the clipboard must be a tab delimited array,
% such as produced by Excel's copy command. If the content does not have
% this structure, the program simply returns a string. If the content is an
% array, x will be a numerical array if all its components qualify as
% numerical. If not, it will be a cell array.
%
% Examples:
%
% 1) If the clipboard contains 'George's job is to chop wood.', then
%    x = paste produces x = 'George's job is to chop wood.'
%
% 2) If the content of the clipboard is a simple text with multiple lines
%    (copied from Notepad or Word or similar), then x = paste produces a
%    cell array with one column and one row per line of the input so each
%    line of text will be separated in different cells. For example, if you
%    copy the follwing text from some other program,
%
%       Manche meinen lechts und rinks kann man nicht velwechsern.
%       Werch ein Illtum!
%
%    then, in Matlab, x = paste produces a 2x1 cell array with
%
%       x{1} = 'Manche meinen lechts und rinks kann man nicht velwechsern.'
%       x{2} = 'Werch ein Illtum!'
%
%    [Note: x = clipboard('copy') would produce just a string in this case,
%    not an array of stringcells, so choose the code that is most useful for
%    your purpose.]
%
% 3) However, if your text contains an equal number of tabs on each line,
%    for instance because you've copied something like this from Word,
%      1  ->  item 1
%      2  ->  item 2
%      3  ->  item 3
%    where -> denotes TABs, then x = paste produces a 3x2 cell array,
%      x = 
%         [1]    'item 1'
%         [2]    'item 2'
%         [3]    'item 3'
%
% 4) If the content of clipboard is an array of cells, e.g.
%         1  2  3
%         4  5  6
%    for instance by copying these six cells from an Excel spreadsheet,
%    then x = paste makes a 2x3 array of doubles with the same content.
%    The same is true if there are NaN cells. So if the Excel excerpt was
%         1     2     3
%         4   #N/A    6
%    then x =
%         1     2     3
%         4   NaN     6
%
% 5) If the cell collection in the clipboard is
%         A  1.3  NaN
%    then x will not be a numerical array, but a 1x3 cell array, with
%     x = 
%        'A'    [1.3000]    [NaN]
%    so x{1} is a string, but x{2} and x{3} are doubles.
%
%
% Here is a practical example:
% ----------------------------
%   In Excel, select your data, say, a sample of observations organized in
%   a few columns. Place them into the clipboard with Ctrl-C.
%   Now switch to Matlab and say
%       x = paste;
%   This puts the data that you copied in Excel into a variable x in
%   Matlab's workspace.
%   Next, you can analyze the data. For instance, compute the principal
%   components (an analysis that is not readily available in Excel), and
%   push the result back into the clipboard,
%       [c,s] = princomp(x);
%       copy(s)
%   Now, back in Excel, you can paste the result into your spreadsheet with
%   Ctrl-V.
%   
%
% Author : Yvan Lengwiler
% Release: $1.03$
% Date   : $2011-06-22$
%
% See also COPY, CLIPBOARD

% History:
%    Version 1.01 (2010-06-25): correction of a bug that occurred with
%       multiple string cells on a single line.
%    Version 1.02 (2011-06-05): Simplified detection of line feeds.
%    Version 1.03 (2011-06-22): Removal of an unused variable.

function x = paste_clipboard()
    
    % get the material from the clipboard
    p = clipboard('paste');
    
    % get out of here if nothing usable is in the clipboard
    % (Note: MLs 'clipboard' interface supports only text, not images or
    % the like.)
    if isempty(p)
        x = [];
        return;
    end
    
    % find linebreaks
    if p(end) ~= char(10)
        p = [p,char(10)];               % append linefeed if missing
    end
    posLF = find(ismember(p,char(10))); % find linefeeds
    nLF   = numel(posLF);               % count linefeeds
    
    % break into separate lines; parse each line by tab
    lines  = cell(nLF,1);
    posTab = cell(nLF,1);
    numTab = zeros(nLF,1);
    last = 0;
    for i = 1:nLF
        lines{i}  = [p(last+1:posLF(i)-1),char(9)]; % append a tabulator
        last      = posLF(i);
        tabs      = ismember(lines{i},char(9));     % find tabulators
        aux       = linspace(1,numel(lines{i}),numel(lines{i}));
        posTab{i} = aux(tabs);                      % positions of tabs
        numTab(i) = sum(tabs(:));                   % count tabs in line
    end

    % is it an array (i.e. a rectangle of cells)?
    isArray = true;
    i = 1;
    while isArray && i <= nLF
        isArray = (numTab(i) == numTab(1));
        i = i+1;
    end
    
    if ~isArray
        % it's not an array, so just return the raw content of the clipboard
        x = p;
        % Note: A simple single or multi-line text with no tabs *does*
        % qualify as an array, so the program splits such content line-wise
        % into a one-column cell array.
    else
        % it is an array, so put it into a Matlab cell array
        isNum = true;   % will remain true if it is never switched off below
        x = cell(nLF,numTab(1));
        for i = 1:nLF
            last = 0;
            pos = posTab{i};
            for j = 1:numTab(1);
                x{i,j} = lines{i,1}(last+1:pos(j)-1);
                % try to make numerical cells if possible
                if ismember(x{i,j},{'NaN','#N/A'})
                    x{i,j} = NaN;
                else
                    aux = str2double(x{i,j});
                    if isnan(aux)
                        % this cell is not numerical (turn off switch for
                        % later)
                        isNum = false;
                    else
                        % str2double has produced a ligit number
                        x{i,j} = aux;
                    end
                end
                last = pos(j);
            end
        end
        if isNum	% make a numerical array if possible
            x = cell2mat(x);
        end
    end
    
    % remove cell encapsulation if there is only one cell
    if numel(x) == 1
        if iscell(x)
            x = x{1};
        end
    end
