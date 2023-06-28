function plot_cross_section(sol)

rhosamp   = sol.rho;
thetasamp = sol.theta;
zsamp   = sol.z;
H1   = sol.H;
nrho = numel(rhosamp);
nth  = numel(thetasamp);
nz   = numel(zsamp);
[~,idy1] = min(abs(thetasamp-0.5*pi));
[~,idy2] = min(abs(thetasamp-1.5*pi));

Rhosamp = [-flip(rhosamp),-rhosamp(1)+1e-3,0,rhosamp(1)-1e-3,rhosamp];
[RHO,Z] = meshgrid(Rhosamp,zsamp);

hx1 = H1(:,idy1,:);
hx1 = reshape(hx1,[nrho,nz]);
hx2 = H1(:,idy2,:);
hx2 = reshape(hx2,[nrho,nz]); 
hmin = min([hx1(:);hx2(:)]);
hmax = max([hx1(:);hx2(:)]);
hart = NaN;
hx1 = [zeros(2,nz)+hart; hx1];
hx2 = [zeros(1,nz)+hart; hx2];
hx2 = flipud(hx2);
hx  = [hx2',hx1'];
nlevels = 40;
figure()
contourf(RHO,Z,hx,nlevels,'LineColor','none');  axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$y$','FontSize', 18, 'interpreter','latex');
ylabel('$z$','FontSize', 18, 'interpreter','latex');
title('$x=0$ (FEM)', 'interpreter','latex');  



[~,idz] = min(abs(zsamp-0.));
hz = H1(:,:,idz);
hz = reshape(hz,[nrho,nth]);
[RHO,TH] = meshgrid(rhosamp,thetasamp);
xx = RHO.*cos(TH);
yy = RHO.*sin(TH);
figure()
contourf(xx,yy,hz',nlevels,'LineColor','none');  axis equal
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 18, 'interpreter','latex');
ylabel('$y$','FontSize', 18,'FontSize', 18, 'interpreter','latex');
title(['$z=',num2str(zsamp(idz)),'$ (FEM)'], 'interpreter','latex');


end