function [xt,yt,zt] = stream_lines_integration(flow_handle,x0,y0,z0,tf)
% forward simulation to generate stream lines

ntrajs = numel(x0);
xt = cell([ntrajs,1]);
yt = cell([ntrajs,1]);
zt = cell([ntrajs,1]);
options = odeset('RelTol',1e-5);
for k=1:ntrajs
    init_state = [x0(k);y0(k);z0(k)];
    [~,hist_state] = ode45(@(t,x) odefun(x,flow_handle),[0,tf],init_state,options);
    xt{k} = hist_state(:,1);
    yt{k} = hist_state(:,2);
    zt{k} = hist_state(:,3);
end

end



function dxdt = odefun(x,flow_handle)

[vx,vy,vz] = flow_handle(x(1),x(2),x(3));
dxdt = [vx;vy;vz];

end