function [vx,vy,vz] = abc_flow(xv,yv,zv)

vx = sqrt(3)*sin(zv)+cos(yv);
vy = sqrt(2)*sin(xv)+sqrt(3)*cos(zv);
vz = sin(yv)+sqrt(2)*cos(xv);

end
