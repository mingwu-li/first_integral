clear all
%% plot leading eigenvalue as a function of number of finite elements
ngrids = [10 20 30 40];
Om = '1';
order = 2;
m = numel(ngrids);
lamd = zeros(m,1);
rsq  = zeros(m,1);
ndof = zeros(m,1);
nele = zeros(m,1);
for k=1:m
    filename = ['cylinderDomainNgrids',num2str(ngrids(k)),'Order',num2str(order),'Om',Om,'.mat'];
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
title('Cylindrical Vortex','FontSize', 18, 'interpreter','latex');


%% plot isosurfaces and level surfaces
% load data
ngrids = 30; 
Om = '1';
om = 1.;
order = 2;
filename = ['cylinderDomainNgrids',num2str(ngrids),'Order',num2str(order),'Om',Om,'.mat'];
sol = load(filename);

% linear polynomial fit
p = polyfit(sol.H(:),sol.Href(:),1);
Hfit = polyval(p,sol.H(:));
Hres = sol.Href(:)-Hfit;
SSresid = sum(Hres.^2);
SStotal = (length(sol.Href(:))-1)*var(sol.Href(:));
rsq     = 1-SSresid/SStotal % r square or coefficient of determination

plot_cross_section(sol,p)

% load velocity profile
[rho,theta,z] = meshgrid(sol.rho,sol.theta,sol.z); % quriy points

H  = reshape(Hfit,size(sol.H));
H  = permute(H,[2,1,3]);
nlevels = 10;
Hmin = min(H(:));
Hmax = max(H(:));
Hsamp = [0.05 0.13];
colors = cool(numel(Hsamp));
tf = zeros(numel(Hsamp),1)+200;
tf(1:2) = [100 50];
for k=1:numel(Hsamp)
    [f,v] = isosurface(rho,theta,z,H,Hsamp(k));
    if ~isempty(v) && (numel(v(:,1))>1 && numel(v(1,:))==3)
        xv = v(:,1).*cos(v(:,2));
        yv = v(:,1).*sin(v(:,2));
        zv = v(:,3);
        figure;
        h = patch('Faces',f,'Vertices',[xv,yv,zv],'FaceAlpha',0.5);
        set(h,'FaceColor',colors(k,:),'EdgeColor','none');
        % forward simulation of stream lines
        x0 = xv([1,3]);
        y0 = yv([1,3]);
        z0 = zv([1,3]);
        [xt,yt,zt] = stream_lines_integration(x0,y0,z0,tf(k),om);
        hold on
        plot3(xt{1},yt{1},zt{1},'k-','LineWidth',1);
        pause(3)
    end
    view(3); xlim([-1 1]); ylim([-1 1]); zlim([-0.4 0.4]); axis equal;
    grid on; camlight; lighting gouraud
    xlabel('$x$','FontSize', 20, 'interpreter','latex');
    ylabel('$y$','FontSize', 20, 'interpreter','latex');
    zlabel('$z$','FontSize', 20, 'interpreter','latex');
    box on
    title(['${H}_2=',num2str(Hsamp(k)),'$'],'FontSize', 20, 'interpreter','latex');    
end
