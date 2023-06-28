function [vx,vy,vz] = euler_flow(x,y,z)

coef = 4*sqrt(2)/(3*sqrt(3));
vx = sin(x-5*pi/6).*cos(y-pi/6).*sin(z)-cos(z-5*pi/6).*sin(x-pi/6).*sin(y);
vy = sin(y-5*pi/6).*cos(z-pi/6).*sin(x)-cos(x-5*pi/6).*sin(y-pi/6).*sin(z);
vz = sin(z-5*pi/6).*cos(x-pi/6).*sin(y)-cos(y-5*pi/6).*sin(z-pi/6).*sin(x);
vx = coef*vx;
vy = coef*vy;
vz = coef*vz;

end
