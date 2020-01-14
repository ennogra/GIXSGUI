function res = res_bvnd(par,xi,yi,zi)
z = bvnd(par,xi,yi);
diff = ((z)-(zi)).^2; diff = diff(:);
diff(isnan(diff)) = [];
res = sum(diff);
end