function [vx,vy,vz] = RBC_flow(xv,yv,zv,velhandles,varargin)

if isempty(varargin)
    vx = velhandles{1}(mod(xv,0.2),yv,mod(zv,0.1));
    vy = velhandles{2}(mod(xv,0.2),yv,mod(zv,0.1));
    vz = velhandles{3}(mod(xv,0.2),yv,mod(zv,0.1));
else
    t  = varargin{1};
    vx = velhandles{1}(mod(xv,0.2),yv,mod(zv,0.1),t);
    vy = velhandles{2}(mod(xv,0.2),yv,mod(zv,0.1),t);
    vz = velhandles{3}(mod(xv,0.2),yv,mod(zv,0.1),t);    
end
end