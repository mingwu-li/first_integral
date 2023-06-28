function [xt,yt,zt] = stream_lines_integration(x0,y0,z0,t0,tf,velhandles)
% forward simulation to generate stream lines

ntrajs = numel(x0);
xt = cell([ntrajs,1]);
yt = cell([ntrajs,1]);
zt = cell([ntrajs,1]);
unsteady = true;
if isempty(t0); t0=0; unsteady = false; end
options = odeset('RelTol',1e-5);
for k=1:ntrajs
    init_state = [x0(k);y0(k);z0(k)];
    [~,hist_state] = ode45(@(t,x) odefun(t,x,velhandles,unsteady),[t0,tf],init_state,options);
    xt{k} = hist_state(:,1);
    yt{k} = hist_state(:,2);
    zt{k} = hist_state(:,3);
end

end



function dxdt = odefun(t,x,velhandles,unsteady)

if unsteady
    [vx,vy,vz] = RBC_flow(x(1),x(2),x(3),velhandles,t);
else
    [vx,vy,vz] = RBC_flow(x(1),x(2),x(3),velhandles);
end
dxdt = [vx;vy;vz];

end