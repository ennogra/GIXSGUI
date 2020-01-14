function b = titlestr(a)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% TITLESTR Convert a string with underscores to a string that can be
%   displayed in figures without recognizing underscores as subscript
%   flags. 
% 
% $Rev 2015/07/27 ZJ
%   input can be a cell

if ischar(a)
    b = convert_title(a);
elseif iscellstr(a)
    b = cell(size(a));
    for ii=1:length(b(:))
        b{ii} = convert_title(a{ii});
    end
else
    error('Invalid input argument.');
end
    
function b = convert_title(a)
b = a;
c = strfind(a,'_');
for i = length(c):-1:1;
    b = [b(1:c(i)-1),'\',b(c(i):end)];
end
    
    
