function plot_cross_section(sol,varargin)

xsamp = sol.x;
zsamp = sol.z;
H = sol.H;
nlevels = 60;
hy = squeeze(H(:,1,:));
[xx,zz] = meshgrid(xsamp,zsamp);
figure()
contour(xx,zz,hy',nlevels,'LineWidth',1); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$z$','FontSize', 20, 'interpreter','latex');
colorbar

end