function y = averaged_relataive_invariance_residual(sol)

[xv,yv,zv] = ndgrid(sol.x,sol.y,sol.z);
[vx,vy,vz] = euler_flow(xv,yv,zv);
v2 = sqrt(vx.^2+vy.^2+vz.^2);
gh = sqrt(sol.GradH(:,:,:,1).^2+sol.GradH(:,:,:,2).^2+sol.GradH(:,:,:,3).^2);

err = sol.GradH(:,:,:,1).*vx+sol.GradH(:,:,:,2).*vy+sol.GradH(:,:,:,3).*vz;

v2 = v2(:); err = err(:); gh = gh(:);
v2(v2==0) = eps;
y = abs(err)./(v2.*gh);
y = mean(y(:));

end
