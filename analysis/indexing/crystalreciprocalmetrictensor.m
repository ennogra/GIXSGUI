function g = crystalreciprocalmetrictensor(a,b,c,alpha,beta,gamma)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% CRYSTALRECIPROCALMETRICTENSOR Calculate the reciprocal metric tensor.
%   CRYSTALRECIPROCALMETRICTENSOR(A,B,C,ALPHA,BETA,GAMMA) calculates the
%   reciprocal metric tensor for all three dimensional crystal systems in 
%   the reciprocal space. (a,b,c) are the length of three basis vector of a
%   Bravais lattice in real space. alpha is the angle between b and c; beta
%   is the angle between a and c; gamma is the angle between a and b. All 
%   angles are in degree.
% 
%   Note: coefficient 2*pi is dropped.
%
%  Reference:
%   1) Marc De Graef and Michael E. McHenry, Structure of Materials: An
%   Introduction to Crystallography, Diffraction and Symmetry, 1st Edition,
%   Cambridge University Press (2007)

%   Zhang Jiang 
%   $Revision: 1.0 $  $Date: 2012/07/31 $

alpha   = alpha*pi/180;
beta    = beta*pi/180;
gamma   = gamma*pi/180;

V2 = a^2*b^2*c^2*(1-cos(alpha)^2-cos(beta)^2-cos(gamma)^2+2*cos(alpha)*cos(beta)*cos(gamma));

g = [...
    b^2*c^2*sin(alpha)^2        a*b*c^2*F(alpha,beta,gamma) a*b^2*c*F(gamma,alpha,beta);...
    a*b*c^2*F(alpha,beta,gamma) a^2*c^2*sin(beta)^2         a^2*b*c*F(beta,gamma,alpha);...
    a*b^2*c*F(gamma,alpha,beta) a^2*b*c*F(beta,gamma,alpha) a^2*b^2*sin(gamma)^2];   
g = g/V2;
g(abs(g)<=eps) = 0;     

function F=F(alpha,beta,gamma)
F=cos(alpha)*cos(beta)-cos(gamma);