function plot_cross_section(sol,varargin)


rhosamp   = sol.rho;
thetasamp = sol.theta;
zsamp   = sol.z;
H1   = sol.H;
Href = sol.Href;
nrho = numel(rhosamp);
nth  = numel(thetasamp);
nz   = numel(zsamp);
if numel(varargin)>0
    fitp = varargin{1};
    H1 = polyval(fitp,H1);
end
Rhosamp = [-flip(rhosamp),rhosamp];
[RHO,Z]   = meshgrid(Rhosamp,zsamp);
yy = RHO;
zz = Z;
[~,idy1] = min(abs(thetasamp-0.5*pi));
[~,idy2] = min(abs(thetasamp-1.5*pi));
hx1 = H1(:,idy1,:);
hx1 = reshape(hx1,[nrho,nz]);
hx2 = H1(:,idy2,:);
hx2 = reshape(hx2,[nrho,nz]);
hx2 = flipud(hx2);
hx  = [hx2',hx1'];
hxref1 = Href(:,idy1,:);
hxref1 = reshape(hxref1,[nrho,nz]);
hxref2 = Href(:,idy2,:);
hxref2 = reshape(hxref2,[nrho,nz]);
hxref2 = flipud(hxref2);
hxref  = [hxref2',hxref1'];
hmin = min([hx(:);hxref(:)]); hmax = max([hx(:);hxref(:)]);
nlevels = 40;
figure()
contourf(yy,zz,hx,nlevels,'LineColor','none'); caxis([hmin,hmax]); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$y$','FontSize', 18, 'interpreter','latex');
ylabel('$z$','FontSize', 18, 'interpreter','latex');
title('$x=0$ (FEniCS)', 'interpreter','latex');
figure()
contourf(yy,zz,hxref,nlevels,'LineColor','none'); caxis([hmin,hmax]); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$y$','FontSize', 18, 'interpreter','latex');
ylabel('$z$','FontSize', 18, 'interpreter','latex');
title('$x=0$ (Reference)','FontSize', 18, 'interpreter','latex');
  


xx = RHO;
zz = Z;
[~,idy3] = min(abs(thetasamp-0));
[~,idy4] = min(abs(thetasamp-pi));
hy1 = H1(:,idy3,:);
hy1 = reshape(hy1,[nrho,nz]);
hy2 = H1(:,idy4,:);
hy2 = reshape(hy2,[nrho,nz]);
hy2 = flipud(hy2);
hy  = [hy2',hy1'];
hyref1 = Href(:,idy3,:);
hyref1 = reshape(hyref1,[nrho,nz]);
hyref2 = Href(:,idy4,:);
hyref2 = reshape(hyref2,[nrho,nz]);
hyref2 = flipud(hyref2);
hyref  = [hyref2',hyref1'];
hmin = min([hy(:);hyref(:)]); hmax = max([hy(:);hyref(:)]);
figure()
contourf(xx,zz,hy,nlevels,'LineColor','none'); caxis([hmin,hmax]); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 18, 'interpreter','latex');
ylabel('$z$','FontSize', 18,'FontSize', 18, 'interpreter','latex');
title('$y=0$ (FEniCS)', 'interpreter','latex');
figure()
contourf(xx,zz,hyref,nlevels,'LineColor','none'); caxis([hmin,hmax]); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 18, 'interpreter','latex');
ylabel('$z$','FontSize', 18,'FontSize', 18, 'interpreter','latex');
title('$y=0$ (Reference)', 'interpreter','latex');

[~,idz] = min(abs(zsamp-0.));
hz = H1(:,:,idz);
hz = reshape(hz,[nrho,nth]);
hzref = Href(:,:,idz);
hzref = reshape(hzref,[nrho,nth]);
hmin  = min([hz(:);hzref(:)]); hmax = max([hz(:);hzref(:)]);
[RHO,TH] = meshgrid(rhosamp,thetasamp);
xx = RHO.*cos(TH);
yy = RHO.*sin(TH);
figure()
contourf(xx,yy,hz',nlevels,'LineColor','none'); caxis([hmin,hmax]); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 18, 'interpreter','latex');
ylabel('$y$','FontSize', 18,'FontSize', 18, 'interpreter','latex');
title('$z=0$ (FEniCS)', 'interpreter','latex');
figure()
contourf(xx,yy,hzref',nlevels,'LineColor','none'); caxis([hmin,hmax]); axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 18, 'interpreter','latex');
ylabel('$y$','FontSize', 18,'FontSize', 18, 'interpreter','latex');
title('$z=0$ (Reference)', 'interpreter','latex');


end