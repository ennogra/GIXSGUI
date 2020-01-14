function A=crystaldirectstructruematrix(a,b,c,alpha,beta,gamma)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% CRYSTALDIRECTSTRUCTUREMATRIX Calculate the direct structure matrix.
%   A=CRYSTALDIRECTSTRUCTUREMATRIX(A,B,C,ALPHA,BETA,GAMMA) calculates the
%   direct structure matrix for the coordinate transformation from Bravais
%   lattice to Cartesian coordinate frame. (a,b,c) are the length of three 
%   basis vectors of a Bravais lattice. alpha is the angle between b and c;
%   beta is the angle between a and c; gamma is the angle between a and b. 
%   All angles are in degree. Then a lattice vector p=[uvw] in the Bravis
%   lattice can be expressed as A*p' in the Cartesian coordinate frame.
%
%  Reference:
%   1) Marc De Graef and Michael E. McHenry, Structure of Materials: An
%   Introduction to Crystallography, Diffraction and Symmetry, 1st Edition,
%   Cambridge University Press (2007)

%   Zhang Jiang 
%   $Revision: 1.0 $  $Date: 2012/07/31 $

% % --- method 1
% g  = crystaldirectmetrictensor(a,b,c,alpha,beta,gamma);
% gr = crystalreciprocalmetrictensor(a,b,c,alpha,beta,gamma);
% 
% alpha   = alpha*pi/180;
% beta    = beta*pi/180;
% gamma   = gamma*pi/180;
% 
% V = sqrt(a^2*b^2*c^2*(1-cos(alpha)^2-cos(beta)^2-cos(gamma)^2+2*cos(alpha)*cos(beta)*cos(gamma)));
% 
% A = [...
%     sqrt(g(1,1))    g(2,1)/sqrt(g(1,1))     g(3,1)/sqrt(g(1,1));...
%     0               V*sqrt(gr(3,3)/g(1,1))  -V*gr(3,2)/sqrt(gr(3,3)*g(1,1));...
%     0               0                       1/sqrt(gr(3,3))];
% 
% --- method 2 (about 3 times faster than method 1)
alpha   = alpha*pi/180;
beta    = beta*pi/180;
gamma   = gamma*pi/180;

V = sqrt(a^2*b^2*c^2*(1-cos(alpha)^2-cos(beta)^2-cos(gamma)^2+2*cos(alpha)*cos(beta)*cos(gamma)));

A = [...
    a   b*cos(gamma)    c*cos(beta);...
    0   b*sin(gamma)    -c*F(beta,gamma,alpha)/sin(gamma);...
    0   0               V/(a*b*sin(gamma))];

A(abs(A)<=eps) = 0;

function F=F(alpha,beta,gamma)
F=cos(alpha)*cos(beta)-cos(gamma);