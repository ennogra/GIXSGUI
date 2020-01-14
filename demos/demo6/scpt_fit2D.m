load H3__th0.130_sec005_up_001_GapFilled_reshaped.mat

info = H3__th00x2E130_sec005_up_001_GapFilled_reshape;
s = double(s);
skip = (info.CountData ==0) ;      % empty pixels

s(skip) = nan;
x = double(info.XData);
y = double(deg2rad(info.YData));

[x2D,y2D] = meshgrid(x,y);

figure
imagesc(x,y,s);
xlabel(info.XVarName);
ylabel(info.YVarName);

nc = length(y);     % number of 1D curves

xy3D(:,:,1) = x2D;
xy3D(:,:,2) = y2D;


% exclude the none rows
xy2D = [x2D(:),y2D(:)];



%% remove nans
z1D = s(:);
ind = isfinite(z1D);

xy2D_new = xy2D(ind,:);
z1D_new = z1D(ind);



% A = par(1);
% Cx = par(2);
% Cy = par(3);
% C = par(4);
% rho = par(5);
% sigma_x = par(6);
% sigma_y = par(7);
% x0 = par(8);
% y0 = par(9);
%               [A,     Cx,     Cy,     C,      rho,    sigma_x     sigma_y,    x0      y0];
par =           [1,     0,      0,      0,      0,      0.1,        0.1,        1.75,   0];
lb =            [eps,   -Inf,   -Inf,   -Inf,   0,     eps,        eps,        1  ,     0]; 
ub =            [Inf,   Inf,    Inf,    Inf,    0,      1,          50,         2,      0];
fitflag = logical([1    1       1       1       0       1           1           1       0]);
surffit = @(par,xy2D)bvld(par,xy2D);

% parfit = par(fitflag);
% parnotfit = par(~fitflag);
getpar = @(parfit,par,fitflag)eval('par(fitflag)=parfit');
%surffitmin = @(parfit,parnotfit,fitflag)sum(abs((log10(bvnd(par,xy2D_new)) - log10(z1D_new))).^2);
surffitmin = @(par)sum(abs((log10(bvld(par,xy2D_new)) - log10(z1D_new))).^2);

%zcal = surffit(par,xy2D_new);

%return;
%[fittedX1,resnorm,residual,exitflag,output,lambda,jacobian] =  lsqcurvefit(surffit,par,xy2D_new,z1D_new,lb,ub);
                    opts = optimset('fminsearch');
                    opts.Display = 'iter';
                    opts.TolX = 1e-6;
                    opts.TolFun = 1e-6;
                    opts.MaxFunEvals = 8000;
                    opts.MaxIter = 6000;
[fittedX1,fval,exitflag,output]=fminsearchbnd(surffitmin,par,lb,ub,opts);

fittedX1
%%
zcal_new = surffit(fittedX1,xy2D_new);

zcal = nan(size(s(:)));
zcal(ind) = zcal_new;
zcal = reshape(zcal,length(y),length(x));

figure
hold on;
plot3(x2D(:),y2D(:),s(:),'o','markerfacecolor','b','markeredgecolor','w','markersize',5,'displayname','data');
surf(x2D,y2D,zcal,'displayname','fitting');
hidden off;
cameratoolbar;
shading flat;
lighting flat;