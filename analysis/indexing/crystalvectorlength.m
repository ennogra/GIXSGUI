function t=crystalvectorlength(p,g)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% CRYSTALVECTORLENGTH Calculate the length of a vector.
%   CRYSTALVECTORLENGTH(P,G) calculates the length of vector. In real
%   space, p=[u,v,w] is the coordinate of the vector and g is the metric
%   tensor. In reciprocal space, p=[h,k,l] is the miller index and g is the
%   reciprocal metric tensor.
%
%  Reference:
%   1) Marc De Graef and Michael E. McHenry, Structure of Materials: An
%   Introduction to Crystallography, Diffraction and Symmetry, 1st Edition,
%   Cambridge University Press (2007)
%
%   See also CRYSTALDIRECTMETRICTENSOR

%   Zhang Jiang 
%   $Revision: 1.0 $  $Date: 2012/07/31 $

p = p(:);
t = sqrt(p'*g*p);