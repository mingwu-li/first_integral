function plot_cross_section(sol,varargin)

rhosamp   = sol.rho;
thetasamp = sol.theta;
psisamp   = sol.psi;
H1   = sol.H;
Href = sol.Href;
nrho = numel(rhosamp);
nth  = numel(thetasamp);
npsi = numel(psisamp);
if numel(varargin)>0
    fitp = varargin{1};
    H1 = polyval(fitp,H1);
end
[RHO,TH] = meshgrid(rhosamp,thetasamp);
yy = RHO.*sin(TH);
zz = RHO.*cos(TH);
[~,idx] = min(abs(psisamp-0.5*pi));
hx = H1(:,:,idx);
hx = reshape(hx,[nrho,nth]);
href = Href(:,:,idx);
href = reshape(href,[nrho,nth]);
hmin = min([hx(:);href(:)]); hmax = max([hx(:);href(:)]);
nlevels = 40;
figure()
% contour(yy,zz,hx',100); axis equal
contourf(yy,zz,hx',nlevels,'LineColor','none'); caxis([hmin,hmax]); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$y$','FontSize', 18, 'interpreter','latex');
ylabel('$z$','FontSize', 18, 'interpreter','latex');
title('$x=0$ (FEM)', 'interpreter','latex');
figure()
contourf(yy,zz,href',nlevels,'LineColor','none'); caxis([hmin,hmax]); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$y$','FontSize', 18, 'interpreter','latex');
ylabel('$z$','FontSize', 18, 'interpreter','latex');
title('$x=0$ (Reference)','FontSize', 18, 'interpreter','latex');


xx = RHO.*sin(TH);
zz = RHO.*cos(TH);
[~,idz] = min(abs(psisamp-0));
hx = H1(:,:,idz);
hx = reshape(hx,[nrho,nth]);
href = Href(:,:,idz);
href = reshape(href,[nrho,nth]);
hmin = min([hx(:);href(:)]); hmax = max([hx(:);href(:)]);
figure()
contourf(xx,zz,hx',nlevels,'LineColor','none'); caxis([hmin,hmax]); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 18, 'interpreter','latex');
ylabel('$z$','FontSize', 18,'FontSize', 18, 'interpreter','latex');
title('$y=0$ (FEM)', 'interpreter','latex');
figure()
contourf(xx,zz,href',nlevels,'LineColor','none'); caxis([hmin,hmax]); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 18, 'interpreter','latex');
ylabel('$z$','FontSize', 18, 'interpreter','latex');
title('$y=0$ (Reference)','FontSize', 18, 'interpreter','latex');

[~,idy] = min(abs(thetasamp-0.5*pi));
hz = H1(:,idy,:);
hz = reshape(hz,[nrho,npsi]);
hzref = Href(:,idy,:);
hzref = reshape(hzref,[nrho,nth]);
hmin = min([hz(:);hzref(:)]); hmax = max([hz(:);hzref(:)]);
[RHO,PSI] = meshgrid(rhosamp,psisamp);
xx = RHO.*cos(PSI);
yy = RHO.*sin(PSI);
figure()
contourf(xx,yy,hz',nlevels,'LineColor','none'); caxis([hmin,hmax]); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 18, 'interpreter','latex');
zlabel('$y$','FontSize', 18, 'interpreter','latex');
title('$z=0$ (FEM)','FontSize', 18, 'interpreter','latex');
figure()
contourf(xx,yy,hzref',nlevels,'LineColor','none'); caxis([hmin,hmax]); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 18, 'interpreter','latex');
ylabel('$y$','FontSize', 18, 'interpreter','latex');
title('$z=0$ (Reference)','FontSize', 18, 'interpreter','latex');

end