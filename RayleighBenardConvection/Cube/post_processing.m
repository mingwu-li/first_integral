clear all

Ra = '1e5'; 
%% check whether the flow is a quasi-2d flow
data = load(['RB_Re',num2str(Ra),'_square_0093000.mat']);
vel  = data.Velocity;
velx = vel(:,:,:,1);
vely = vel(:,:,:,2);
velz = vel(:,:,:,3);
% show velz is small relative to the other two components
fprintf('norm of vx is %d\n', norm(velx(:)))
fprintf('norm of vy is %d\n', norm(vely(:)))
fprintf('norm of vz is %d\n', norm(velz(:)))

velx_mean = mean(velx,3); % mean along z direction
vely_mean = mean(vely,3); % mean along z direction
[nx,ny,nz] = size(velx);
dvelx = zeros(nx,ny,nz); dvely = zeros(nx,ny,nz); 
for k=1:nz
    dvelx(:,:,k) = abs(velx(:,:,k)-velx_mean);
    dvely(:,:,k) = abs(vely(:,:,k)-vely_mean);
end
dvelx = mean(dvelx,3); dvely = mean(dvely,3); 
% show deviation of velx/y from its mean is tiny
fprintf('norm of |vx-<vx>_z| is %d\n', norm(dvelx(:)))
fprintf('norm of |vy-<vy>_z| is %d\n', norm(dvely(:)))

% contour plot of velocity profiles
[xx,yy,zz] = meshgrid(data.xgrid,data.ygrid,data.zgrid);
nlevels = 30;
velxmin = min(velx(:));
velxmax = max(velx(:));
velxsamp = linspace(velxmin,velxmax,nlevels);
figure; 
for k=1:numel(velxsamp)
    velxk = velxsamp(k);
    isosurface(xx,yy,zz,permute(velx,[2,1,3]),velxk); hold on
end
colorbar
axis equal; box on
axis([0 0.1 0 0.1]);
xlabel('x'); ylabel('y'); zlabel('z'); title('vx');
set(gca,'FontSize',14); set(gca,'LineWidth',1.5);


velymin = min(vely(:));
velymax = max(vely(:));
velysamp = linspace(velymin,velymax,nlevels);
figure; 
for k=1:numel(velysamp)
    velyk = velysamp(k);
    isosurface(xx,yy,zz,permute(vely,[2,1,3]),velyk); hold on
end
colorbar
axis equal; box on
axis([0 0.1 0 0.1]);
xlabel('x'); ylabel('y'); zlabel('z'); title('vy');
set(gca,'FontSize',14); set(gca,'LineWidth',1.5);


velzmin = min(velz(:));
velzmax = max(velz(:));
velzsamp = linspace(velzmin,velzmax,nlevels);
figure; 
for k=1:numel(velzsamp)
    velzk = velzsamp(k);
    isosurface(xx,yy,zz,permute(velz,[2,1,3]),velzk); hold on
end
colorbar
axis equal; box on
axis([0 0.1 0 0.1 0 0.1]);
xlabel('x'); ylabel('y'); zlabel('z'); title('vz');
set(gca,'FontSize',14); set(gca,'LineWidth',1.5);


%% contour and isosurface plots of first integral results
% load first integral results
ny    = 30;
az    = 0.2;
order = 2;
bcs   = 'DPBC'; 
filename = ['RBCnx',num2str(2*ny),'ny',num2str(ny),'nz',num2str(ny),'Order',num2str(order),bcs,Ra,'.mat'];
sol = load(filename);
vel     = permute(vel,[2,1,3,4]);
vel     = double(vel);
[xvel,yvel,zvel] = meshgrid(data.xgrid,data.ygrid,data.zgrid);
[xv,yv,zv] = meshgrid(sol.x,sol.y,sol.z); % quriy points
vx = interp3(xvel,yvel,zvel,vel(:,:,:,1),xv,yv,zv);
vy = interp3(xvel,yvel,zvel,vel(:,:,:,2),xv,yv,zv);
vz = interp3(xvel,yvel,zvel,vel(:,:,:,3),xv,yv,zv);
v2 = sqrt(vx.^2+vy.^2+vz.^2);
% interpolator for GradH and vnorm
GradH1 = permute(sol.GradH1,[2,1,3,4]);
gradH1 = sqrt(GradH1(:,:,:,1).^2+GradH1(:,:,:,2).^2+GradH1(:,:,:,3).^2);
Err1   = GradH1(:,:,:,1).*vx+GradH1(:,:,:,2).*vy+GradH1(:,:,:,3).*vz;
Err1   = abs(Err1);
[Xv,Yv,Zv] = ndgrid(sol.x,sol.y,sol.z); % quriy points
Err1Int    = griddedInterpolant(Xv,Yv,Zv,permute(Err1,[2,1,3]));
GradH1Int  = griddedInterpolant(Xv,Yv,Zv,permute(gradH1,[2,1,3]));
vnormInt   = griddedInterpolant(Xv,Yv,Zv,permute(v2,[2,1,3]));

H1  = sol.H1;
H1  = permute(H1,[2,1,3]);
nlevels = 50;
Hmin = min(H1(:));
Hmax = max(H1(:));
Hsamp = linspace(Hmin,Hmax,nlevels);
threshold = 0.01; % filter threshold
Hg1 = filter_H(Err1Int,GradH1Int,vnormInt,H1,Hsamp,xv,yv,zv,filename,threshold);
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14); grid on;
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$y$','FontSize', 20, 'interpreter','latex');
zlabel('$z$','FontSize', 20, 'interpreter','latex');
title(['RBC $(Ra=1\times10^5, E_A\leq',num2str(threshold),')$'],'FontSize', 20, 'interpreter','latex');

% lauching trajectories 
% frozen time
vxInt = griddedInterpolant(Xv,Yv,Zv,permute(vx,[2,1,3]));
vyInt = griddedInterpolant(Xv,Yv,Zv,permute(vy,[2,1,3]));
vzInt = griddedInterpolant(Xv,Yv,Zv,permute(vz,[2,1,3]));
% time-varying
filevel = ['RB_Re',num2str(Ra),'_square_time_persec.mat'];
solvel  = load(filevel);
ve = solvel.Velocity; ve = double(ve);
ux = ve(:,:,:,:,1); ux = permute(ux,[2,3,4,1]);
uy = ve(:,:,:,:,2); uy = permute(uy,[2,3,4,1]);
uz = ve(:,:,:,:,3); uz = permute(uz,[2,3,4,1]);
ux = ux(2:end-1,2:end-1,2:end-1,:);
uy = uy(2:end-1,2:end-1,2:end-1,:);
uz = uz(2:end-1,2:end-1,2:end-1,:);
xx = solvel.xgrid(2:end-1);
yy = solvel.ygrid(2:end-1);
zz = solvel.zgrid(2:end-1);
[Xvt,Yvt,Zvt,Tt] = ndgrid(xx,yy,zz,solvel.ts); % quriy points
vxtInt = griddedInterpolant(Xvt,Yvt,Zvt,Tt,ux);
vytInt = griddedInterpolant(Xvt,Yvt,Zvt,Tt,uy);
vztInt = griddedInterpolant(Xvt,Yvt,Zvt,Tt,uz);
% plot isosurface
colors = get(0,'defaultaxescolororder');
figure;
[f,v] = isosurface(xv,yv,zv,H1,Hg1(1)); hold on
h = patch('Faces',f,'Vertices',v,'FaceAlpha',0.8); grid on; box on
set(h,'FaceColor',colors(3,:),'EdgeColor','none');
axis equal; view([-37 42]);
axis([0 0.1 0 0.1 0 0.1]);
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$y$','FontSize', 20, 'interpreter','latex');
zlabel('$z$','FontSize', 20, 'interpreter','latex');
camlight; lighting gouraud
      
fk = isosurface(xv,yv,zv,H1,Hg1(1));
vk = fk.vertices;
npts = size(vk,1);
load('Ra1e5_init_h1.mat');
x0 = vk(idxs,1);
y0 = vk(idxs,2);
z0 = vk(idxs,3);
tf = 20;
% steady flow (frozen time)
[xt,yt,zt] = stream_lines_integration(x0,y0,z0,[],tf,{vxInt,vyInt,vzInt});
t0 = solvel.ts(621); % for Ra1e5_0092250
tf = 20;
% unsteady flow (time-varying)
[xut,yut,zut] = stream_lines_integration(x0,y0,z0,t0,t0+tf,{vxtInt,vytInt,vztInt});
hold on
for k=1:10
    ztk = mod(zt{k},0.1);
    dif = ztk(2:end)-ztk(1:end-1);
    dif = abs(dif);
    idk = find(dif>0.05);
    if ~isempty(idk)
        idxs = [0;idk];
        for j=1:numel(idk)
            plot3(xt{k}(idxs(j)+1:idxs(j+1)),yt{k}(idxs(j)+1:idxs(j+1)),...
                ztk(idxs(j)+1:idxs(j+1)),'k-','LineWidth',1);
        end
    else
        plot3(xt{k},yt{k},mod(zt{k},0.1),'k-','LineWidth',1);
    end
end
axis([0.0 0.1 0 0.1 0 0.1]);

hold on
for k=1:10
    ztk = mod(zut{k},0.1);
    dif = ztk(2:end)-ztk(1:end-1);
    dif = abs(dif);
    idk = find(dif>0.05);
    if ~isempty(idk)
        idxs = [0;idk];
        for j=1:numel(idk)
            plot3(xut{k}(idxs(j)+1:idxs(j+1)),yut{k}(idxs(j)+1:idxs(j+1)),...
                ztk(idxs(j)+1:idxs(j+1)),'b-','LineWidth',1);
        end
    else
        plot3(xut{k},yut{k},mod(zut{k},0.1),'b-','LineWidth',1);
    end
end
axis([0.0 0.1 0 0.1 0 0.1]);

