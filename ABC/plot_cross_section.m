function plot_cross_section(sol,varargin)
% PLOT_CROSS_SECTION Plot contour plot at cross section x=0, y=0 and z=0

% x=xa (the left end point)
xsamp = sol.x;
ysamp = sol.y;
zsamp = sol.z;
H = sol.H;
nlevels = 30;
[yy,zz] = meshgrid(ysamp,zsamp);
figure()
contourf(yy,zz,squeeze(H(1,:,:))',nlevels,'LineColor','none'); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$y$','FontSize', 20, 'interpreter','latex');
ylabel('$z$','FontSize', 20, 'interpreter','latex');
colorbar

% y=ya (the left end point)
hy = squeeze(H(:,1,:));
[xx,zz] = meshgrid(xsamp,zsamp);
figure()
contourf(xx,zz,hy',nlevels,'LineColor','none'); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$z$','FontSize', 20, 'interpreter','latex');
colorbar

% z=za (the left end point)
hz = squeeze(H(:,:,1));
[xx,yy] = meshgrid(xsamp,ysamp);
figure()
contourf(xx,yy,hz',nlevels,'LineColor','none'); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$y$','FontSize', 20, 'interpreter','latex');
colorbar

end