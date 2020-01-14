function y = latticetype(sg)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% LATTICETYPE Returns lattice type of space group.

switch sg
    case num2cell(1:2)
        y = 'Triclinic';
    case num2cell(3:15)
        y = 'Monoclinic';
    case num2cell(16:74)
        y = 'Orthorhombic';
    case num2cell(75:142)
        y = 'Tetragonal';
    case num2cell(143:167)
        y = 'Trigonal';
    case num2cell(168:194)
        y = 'Hexagonal';
    case num2cell(195:230)
        y = 'Cubic';
    otherwise 
        error('Invalid space group');
end
            