function z = bvnd(par,xy)
% Bivariate normal distribution function
A = par(1);
Cx = par(2);
Cy = par(3);
C = par(4);
rho = par(5);
sigma_x = par(6);
sigma_y = par(7);
x0 = par(8);
y0 = par(9);
%xi = xy(:,:,1);
%yi = xy(:,:,2);
xi = xy(:,1);
yi = xy(:,2);
z = Cx*xi+Cy*yi+C+A*exp(-1/(2*(1-rho)^2)*( (xi-x0).^2/sigma_x^2 + (yi-y0).^2/sigma_y^2 - 2*rho*(xi-x0).*(yi-y0)/(sigma_x*sigma_y)));
