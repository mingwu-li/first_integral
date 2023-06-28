function y = averaged_relataive_invariance_residual(sol,flow_handle)
% AVERAGED_RELATIVE_INVARIANCE_RESIDUAL This function calculates averaged
% relative invariance residual. Relative invariance residual is defined as
% |gradH*v|/(|v|*|gradH|). Such a relative invariance residual is evaluated
% at a collection of grid points and then this function returns the mean or
% average of these relative invariance residuals.

% extract data and calculate v and gradH
[xv,yv,zv] = ndgrid(sol.x,sol.y,sol.z);
[vx,vy,vz] = flow_handle(xv,yv,zv);
v2 = sqrt(vx.^2+vy.^2+vz.^2);
gh = sqrt(sol.GradH(:,:,:,1).^2+sol.GradH(:,:,:,2).^2+sol.GradH(:,:,:,3).^2);
% extract mean relative residual
err = sol.GradH(:,:,:,1).*vx+sol.GradH(:,:,:,2).*vy+sol.GradH(:,:,:,3).*vz;
y = abs(err)./(v2.*gh);
y = mean(y(:));

end
