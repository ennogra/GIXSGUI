function g = crystaldirectmetrictensor(a,b,c,alpha,beta,gamma)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% CRYSTALDIRECTMETRICTENSOR Calculate the crystal metric tensor.
%   CRYSTALDIRECTMETRICTENSOR(A,B,C,ALPHA,BETA,GAMMA) calculates the metric
%   tensor for a three dimensional crystal system. (a,b,c) are the
%   length of three basis vector of a Bravais lattice. alpha is the angle 
%   between b and c; beta is the angle between a and c; gamma is the angle
%   between a and b. All angles are in degree.
%
%  Reference:
%   1) Marc De Graef and Michael E. McHenry, Structure of Materials: An
%   Introduction to Crystallography, Diffraction and Symmetry, Chapter 3
%   and 4, 1st Edition, Cambridge University Press (2007)

%   Zhang Jiang 
%   $Revision: 1.0 $  $Date: 2012/07/31 $

alpha   = alpha*pi/180;
beta    = beta*pi/180;
gamma   = gamma*pi/180;

g = [...
    a^2             a*b*cos(gamma)  a*c*cos(beta);...
    b*a*cos(gamma)  b^2             b*c*cos(alpha);...
    c*a*cos(beta)   c*b*cos(alpha)  c^2];
g(abs(g)<=eps) = 0;  