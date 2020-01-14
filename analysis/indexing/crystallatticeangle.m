function theta = crystallatticeangle(p1,p2,g)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% CRYSTALLATTICEANGLE Calculate the angle between lattice vectors or planes.
%   CRYSTALLATTICEANGLE(P1,P2,G) calculates the angle between two lattice 
%   vectors with coordinates P1 and P2, or lattice planes with plane index
%   P1 and P2. For angle between vectors, p=[u,v,w] and g is the real space
%   metric tensor of the Bravais lattice. For angle between planes, p is the
%   miller index [h,k,l] for the plane, and g is the reciprocal space
%   metric tensor. 
%
%  Reference:
%   1) Marc De Graef and Michael E. McHenry, Structure of Materials: An
%   Introduction to Crystallography, Diffraction and Symmetry, 1st Edition,
%   Cambridge University Press (2007)
%
%   See also CRYSTALDIRECTMETRICTENSOR, CRYSTALRECIPROCALMETRICTENSOR

%   Zhang Jiang 
%   $Revision: 1.0 $  $Date: 2012/07/31 $

p1=p1(:);
p2=p2(:);
theta = 180/pi*acos(p1'*g*p2/sqrt(p1'*g*p1)/sqrt(p2'*g*p2));