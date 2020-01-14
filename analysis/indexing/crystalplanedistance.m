function d=crystalplanedistance(p,g)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% CRYSTALPLANEDISTANCE Calculate the crystal plane distance.
%   CRYSTALPLANEDISTANCE(P,G) calculates the plane distance of parallel
%   planes with miller index p=[h,k,l]. g is the reciprocal metric tensor.
%
%  Reference:
%   1) Marc De Graef and Michael E. McHenry, Structure of Materials: An
%   Introduction to Crystallography, Diffraction and Symmetry, 1st Edition,
%   Cambridge University Press (2007)
%
%   See also CRYSTALMETRICTENSOR

%   Zhang Jiang 
%   $Revision: 1.0 $  $Date: 2012/07/31 $

p = p(:);
G = sqrt(p'*g*p);
d=1/G;
