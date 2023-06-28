clear all

%% plot leading eigenvalue as a function of number of finite elements
ngrids = [20 30 40 45];
Om = '1';
om = 1.;
order = 2;
m = numel(ngrids);
lamd = zeros(m,1);
rsq  = zeros(m,1);
ndof = zeros(m,1);
nele = zeros(m,1);
for k=1:m
    filename = ['TC_Free_Zseg_Ngrids',num2str(ngrids(k)),'Order',num2str(order),'.mat'];
    sol = load(filename);
    lamd(k) = sol.lambda(2);
    ndof(k) = sol.numDOF;
    nele(k) = sol.numEle;
end

figure;
loglog(nele,lamd,'ro-','LineWidth',2,'MarkerSize',10); grid on;
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('Number of Elements','FontSize', 18, 'interpreter','latex');
ylabel('$\lambda_2$','FontSize', 18, 'interpreter','latex');
title('Cylindrical Vortex','FontSize', 18, 'interpreter','latex');

%% plot isosurfaces and level surfaces
% load data
ngrids = 40; 
order  = 2;
filename = ['TC_Free_Zseg_Ngrids',num2str(ngrids),'Order',num2str(order),'.mat'];
sol = load(filename);
% contour plot at cross sections
plot_cross_section(sol);
[rho,theta,z] = meshgrid(sol.rho,sol.theta,sol.z); % quriy points
H  = sol.H;
H  = abs(H); % take the absolute value
H  = permute(H,[2,1,3]);
nlevels = 20;
Hmin  = min(H(:));
Hmax  = max(H(:));
Hsamp = [0.2 0.5];
% load velocity profile and perform interpolation
solvel   = load('taylorvortex_rphiz_grid.mat');
interp_method = 'linear'; % 'linear'
u_interp = griddedInterpolant({solvel.rgrid,solvel.phigrid,solvel.zgrid},...
    solvel.vx,interp_method,'nearest');
v_interp = griddedInterpolant({solvel.rgrid,solvel.phigrid,solvel.zgrid},...
    solvel.vy,interp_method,'nearest');
w_interp = griddedInterpolant({solvel.rgrid,solvel.phigrid,solvel.zgrid},...
    solvel.vz,interp_method,'nearest');
uinterp = @(x,y,z) u_interp(x,y,mod(z,pi));
vinterp = @(x,y,z) v_interp(x,y,mod(z,pi));
winterp = @(x,y,z) w_interp(x,y,mod(z,pi));
% plot isosurfaces
colors = cool(numel(Hsamp));
tf = zeros(numel(Hsamp),1)+500;
for k=1:numel(Hsamp)
    [f,v] = isosurface(rho,theta,z,H,Hsamp(k));
    if ~isempty(v) && (numel(v(:,1))>1 && numel(v(1,:))==3)
        xv = v(:,1).*cos(v(:,2));
        yv = v(:,1).*sin(v(:,2));
        zv = v(:,3);
        [~,id1] = min(zv); [~,id2] = max(zv);
        figure;
        h = patch('Faces',f,'Vertices',[xv,yv,zv]);
        set(h,'FaceColor',colors(k,:),'EdgeColor','none','FaceAlpha',0.5);
        % forward simulation of stream lines
        x0 = xv([id1,id2]);
        y0 = yv([id1,id2]);
        z0 = zv([id1,id2]);
        [xt,yt,zt] = stream_lines_integration(x0,y0,z0,tf(k),...
            uinterp,vinterp,winterp);
        hold on
        plot3(xt{1},yt{1},zt{1},'k-','LineWidth',1);
        plot3(xt{2},yt{2},zt{2},'k-','LineWidth',1);
        pause(3)
    end
    hold off
    view(3); xlim([-1 1]); ylim([-1 1]); zlim([0 pi]); axis equal;
    grid on; camlight; lighting gouraud
    title(['${H}_2=',num2str(Hsamp(k)),'$'],'FontSize', 20, 'interpreter','latex');
end
