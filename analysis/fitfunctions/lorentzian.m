% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
function y = lorentzian(x,a,b,c)
y = c/pi*b./((x-a).^2+b^2);