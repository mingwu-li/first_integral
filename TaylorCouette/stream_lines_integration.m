function [xt,yt,zt] = stream_lines_integration(x0,y0,z0,tf,vxInterp,vyInterp,vzInterp)
% forward simulation to generate stream lines

ntrajs = numel(x0);
xt = cell([ntrajs,1]);
yt = cell([ntrajs,1]);
zt = cell([ntrajs,1]);
options = odeset('RelTol',1e-5);
for k=1:ntrajs
    init_state = [x0(k);y0(k);z0(k)];
    [~,hist_state] = ode45(@(t,x) odefun(x,vxInterp,vyInterp,vzInterp),...
        [0,tf],init_state,options);
    xt{k} = hist_state(:,1);
    yt{k} = hist_state(:,2);
    zt{k} = hist_state(:,3);
end

end



function dxdt = odefun(x,vxInterp,vyInterp,vzInterp)

[th,r] = cart2pol(x(1),x(2));
if th<0; th = th+2*pi; end

dxdt = zeros(3,1);
dxdt(1) = vxInterp(r,th,x(3));
dxdt(2) = vyInterp(r,th,x(3));
dxdt(3) = vzInterp(r,th,x(3));

end

