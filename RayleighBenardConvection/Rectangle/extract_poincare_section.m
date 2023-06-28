function [x,y,z] = extract_poincare_section(xt,yt,zt,pfun)
% EXTRACT_POINCARE_SECTION This function extracts the intersection points
% of a 3-dimensional trajectory (xt,yt,zt) with a Poincare section
% specified by pfun. Here pfun is a function handle.

% previous and current steps 
xp = xt(1:end-1); xc = xt(2:end);
yp = yt(1:end-1); yc = yt(2:end);
zp = zt(1:end-1); zc = zt(2:end);

% evalute of pfun
pfunp = pfun(xp,yp,zp);
pfunc = pfun(xc,yc,zc);

% check sign change
fun_sign = pfunp.*pfunc;
idx = find(fun_sign<0);

% take midpoints as intersection points
x = 0.5*(xt(idx)+xt(idx+1));
y = 0.5*(yt(idx)+yt(idx+1));
z = 0.5*(zt(idx)+zt(idx+1));

end