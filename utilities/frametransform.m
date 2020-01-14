function [x1,y1,z1] = frametransform(T,x,y,z)
% ***********************************************
% Copyright (c) 2017 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% FRAMETRANSFORM Rotate frame. 
%   [X1,Y1,Z1] = FRAMETRANSFORM(T,X,Y,Z) calculates new vector (X1,Y1,Z1) 
%   in the transformed Cartesian coordinate system. (X,Y,Z) is the vector 
%   in in the old Cartesian coordinate system. R (3x3) is coordinate 
%   transformation matrix. X,Y,Z can be arrays.

%   Zhang Jiang 
%   $Revision: 1.0 $  $Date: 2012/07/27 $

T=inv(T);

x1 = T(1,1)*x + T(1,2)*y + T(1,3)*z;
y1 = T(2,1)*x + T(2,2)*y + T(2,3)*z;
z1 = T(3,1)*x + T(3,2)*y + T(3,3)*z;

