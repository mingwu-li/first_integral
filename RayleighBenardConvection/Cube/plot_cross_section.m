function plot_cross_section(sol,varargin)

xsamp = sol.x;
ysamp = sol.y;
zsamp = sol.z;
H = sol.H;
nlevels = 60;
% [yy,zz] = meshgrid(ysamp,zsamp);
% figure()
% contour(yy,zz,squeeze(H(1,:,:))',nlevels,'LineWidth',1); axis equal
% % contourf(yy,zz,squeeze(H(1,:,:))',nlevels,'LineColor','none'); axis equal
% set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
% xlabel('$y$','FontSize', 20, 'interpreter','latex');
% ylabel('$z$','FontSize', 20, 'interpreter','latex');
% colorbar
% title('$x=0$ (FEM)', 'interpreter','latex');

hy = squeeze(H(:,1,:));
[xx,zz] = meshgrid(xsamp,zsamp);
figure()
contour(xx,zz,hy',nlevels,'LineWidth',1); axis equal
% contourf(xx,zz,hy',nlevels,'LineColor','none'); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$z$','FontSize', 20, 'interpreter','latex');
colorbar
% 
% hz = squeeze(H(:,:,1));
% [xx,yy] = meshgrid(xsamp,ysamp);
% figure()
% contourf(xx,yy,hz',nlevels,'LineColor','none'); axis equal
% set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
% xlabel('$x$','FontSize', 20, 'interpreter','latex');
% ylabel('$y$','FontSize', 20, 'interpreter','latex');
% colorbar

end