function G = gvector(a,miller)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% GVECTOR Calculate the reciprocal lattice vectors
%   G = GVECTOR(A,MILLER) calculates the reciprocal lattice vectors of
%   miller indices (1x3 vector) for a real space lattice with lattice 
%   vectors a=(a1,a2,a3) whose coordinates are defined in a Cartesian 
%   reference system.

%   Zhang Jiang @8ID/APS/ANL
%   $Revision: 1.0 $  $Date: 2012/08/08 $


% --- real space lattce vector
a1 = a(:,1);    % 3x1
a2 = a(:,2);    % 3x1
a3 = a(:,3);    % 3x1

% --- basis of reciprocal lattice
b1 = 2*pi*cross(a2,a3)/dot(a1,cross(a2,a3));        % 3x1
b2 = 2*pi*cross(a3,a1)/dot(a1,cross(a2,a3));        % 3x1
b3 = 2*pi*cross(a1,a2)/dot(a1,cross(a2,a3));        % 3x1

% --- miller index
h = miller(1);
k = miller(2);
l = miller(3);

% --- reciprocal lattice vector (3x1)
G=h*b1+k*b2+l*b3;