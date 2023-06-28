clear all

%% plot leading eigenvalue as a function of number of finite elements
ngrids = [10 20 30 40];
c = '0p1';
order = 2;
m = numel(ngrids);
lamd = zeros(m,1);
rsq  = zeros(m,1);
ndof = zeros(m,1);
nele = zeros(m,1);
for k=1:m
    filename = ['SphereDomainNgrids',num2str(ngrids(k)),'Order',num2str(order),'c',c,'.mat'];
    sol = load(filename);
    lamd(k) = sol.lambda(2);
    ndof(k) = sol.numDOF;
    nele(k) = sol.numEle;
    p = polyfit(sol.H(:),sol.Href(:),1);
    Hfit = polyval(p,sol.H(:));
    Hres = sol.Href(:)-Hfit;
    SSresid = sum(Hres.^2);
    SStotal = (length(sol.Href(:))-1)*var(sol.Href(:));
    rsq(k) = 1-SSresid/SStotal; % r square or coefficient of determination
end

figure;
loglog(nele,lamd,'ro-','LineWidth',2,'MarkerSize',10); grid on;
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('Number of Elements','FontSize', 18, 'interpreter','latex');
ylabel('$\lambda_2$','FontSize', 18, 'interpreter','latex');
title('Spherical Vortex','FontSize', 18, 'interpreter','latex');


%% plot isosurfaces and level surfaces
% load data
ngrids = 30; 
coeff = '0p1';
c = 0.1;
order = 2;
filename = ['SphereDomainngrids',num2str(ngrids),'order',num2str(order),'c',coeff,'.mat'];
sol = load(filename);
% linear polynomial fit
p = polyfit(sol.H(:),sol.Href(:),1);
Hfit = polyval(p,sol.H(:));
Hres = sol.Href(:)-Hfit;
SSresid = sum(Hres.^2);
SStotal = (length(sol.Href(:))-1)*var(sol.Href(:));
rsq     = 1-SSresid/SStotal % r square or coefficient of determination
plot_cross_section(sol,p)
% load velocity profile as well
[rho,theta,psi] = meshgrid(sol.rho,sol.theta,sol.psi); % quriy points
H  = reshape(Hfit,size(sol.H));
H  = permute(H,[2,1,3]);
Hmin = min(H(:));
Hmax = max(H(:));
Hsamp = [0.08 0.12];
colors = cool(numel(Hsamp));
tf = zeros(numel(Hsamp),1)+200;
tf(1:2) = [200 100];
for k=1:numel(Hsamp)
    [f,v] = isosurface(rho,theta,psi,H,Hsamp(k));
    figure;
    if numel(v(:,1))>1 && numel(v(1,:))==3
        xv = v(:,1).*sin(v(:,2)).*cos(v(:,3));
        yv = v(:,1).*sin(v(:,2)).*sin(v(:,3));
        zv = v(:,1).*cos(v(:,2));
        h = patch('Faces',f,'Vertices',[xv,yv,zv],'FaceAlpha',0.5);
        set(h,'FaceColor',colors(k,:),'EdgeColor','none');
        % forward simulation of stream lines
        x0 = xv([1,3]);
        y0 = yv([1,3]);
        z0 = zv([1,3]);
        [xt,yt,zt] = stream_lines_integration(x0,y0,z0,tf(k),c);
        hold on
        plot3(xt{1},yt{1},zt{1},'k-','LineWidth',1);
        plot3(xt{2},yt{2},zt{2},'k-','LineWidth',1);
        pause(3)
    end
    view(3); xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
    grid on; camlight; lighting gouraud
    xlabel('$x$','FontSize', 20, 'interpreter','latex');
    ylabel('$y$','FontSize', 20, 'interpreter','latex');
    zlabel('$z$','FontSize', 20, 'interpreter','latex');
    box on
    title(['$\hat{H}_2=',num2str(Hsamp(k)),'$'],'FontSize', 20, 'interpreter','latex');
end

