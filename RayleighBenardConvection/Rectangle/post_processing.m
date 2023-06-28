clear all

%% Ra = 5e4 - steady flow
Ra    = '5e4'; ny = 25;
az    = 0.2;
order = 2;
bcs   = 'DPBC'; % PBC
filename = ['RBCnx',num2str(2*ny),'ny',num2str(ny),'nz',num2str(ny),'Order',num2str(order),bcs,Ra,'.mat'];

sol = load(filename);

% load velocity profile as well
filevel = ['RB_Re',num2str(Ra),'.mat'];
% filevel = ['RB_Re',num2str(Ra),'_0092250.mat'];
solvel  = load(filevel);
vel     = solvel.Velocity;
vel     = permute(vel,[2,1,3,4]);
vel     = double(vel);
[xvel,yvel,zvel] = meshgrid(solvel.xgrid,solvel.ygrid,solvel.zgrid);
[xv,yv,zv] = meshgrid(sol.x,sol.y,sol.z); % quriy points
vx = interp3(xvel,yvel,zvel,vel(:,:,:,1),xv,yv,zv);
vy = interp3(xvel,yvel,zvel,vel(:,:,:,2),xv,yv,zv);
vz = interp3(xvel,yvel,zvel,vel(:,:,:,3),xv,yv,zv);
v2 = sqrt(vx.^2+vy.^2+vz.^2);

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
% threshold = 10.05;
% threshold = 0.05;
threshold = 0.005;
Hg1 = filter_H(Err1Int,GradH1Int,vnormInt,H1,Hsamp,xv,yv,zv,filename,threshold);
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14); grid on;
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$y$','FontSize', 20, 'interpreter','latex');
zlabel('$z$','FontSize', 20, 'interpreter','latex');
title('RBC $(H_1, Ra=5\times10^4)$','FontSize', 20, 'interpreter','latex');

% to see wheather the data is a quasi-2d data
data = load(['RB_Re',Ra,'.mat']);
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

% plot of streamlines in cross section (z=constant)
[xx,yy] = meshgrid(data.xgrid,data.ygrid);
figure; hold on
streamslice(xx,yy,transpose(velx(:,:,1)),transpose(vely(:,:,1)),2); grid on
axis equal; box on
axis([0 az 0 0.1]);
xlabel('x'); ylabel('y');
set(gca,'FontSize',14); set(gca,'LineWidth',1.5);
[xx,yy] = meshgrid(sol.x,sol.y);
contour(xx,yy,H1(:,:,1),50);

% we can also look contour plots to check whether the flow is 2D
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
axis([0 az 0 0.1]);
xlabel('x'); ylabel('y'); zlabel('z'); title('ux');
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
axis([0 az 0 0.1]);
xlabel('x'); ylabel('y'); zlabel('z'); title('uy');
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
axis([0 az 0 0.1 0 0.1]);
xlabel('x'); ylabel('y'); zlabel('z'); title('uz');
set(gca,'FontSize',14); set(gca,'LineWidth',1.5);


%% Ra = 1e5 - unsteady flow
Ra    = '1e5'; ny = 30; 
az    = 0.2;
order = 2;
bcs   = 'DPBC'; % PBC
filename  = ['RBCnx',num2str(2*ny),'ny',num2str(ny),'nz',num2str(ny),'Order',num2str(order),bcs,Ra,'.mat'];

sol = load(filename);

% load velocity profile as well
filevel = ['RB_Re',num2str(Ra),'_0092250.mat'];
solvel  = load(filevel);
vel     = solvel.Velocity;
vel     = permute(vel,[2,1,3,4]);
vel     = double(vel);
[xvel,yvel,zvel] = meshgrid(solvel.xgrid,solvel.ygrid,solvel.zgrid);
[xv,yv,zv] = meshgrid(sol.x,sol.y,sol.z); % quriy points
vx = interp3(xvel,yvel,zvel,vel(:,:,:,1),xv,yv,zv);
vy = interp3(xvel,yvel,zvel,vel(:,:,:,2),xv,yv,zv);
vz = interp3(xvel,yvel,zvel,vel(:,:,:,3),xv,yv,zv);
v2 = sqrt(vx.^2+vy.^2+vz.^2);

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
threshold = 0.005;
Hg1 = filter_H(Err1Int,GradH1Int,vnormInt,H1,Hsamp,xv,yv,zv,filename,threshold);
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14); grid on;
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$y$','FontSize', 20, 'interpreter','latex');
zlabel('$z$','FontSize', 20, 'interpreter','latex');
title('RBC $(H_1, Ra=1\times10^5, E_A\leq0.005)$','FontSize', 20, 'interpreter','latex');

if isfield(sol,'H2')
    GradH2 = permute(sol.GradH2,[2,1,3,4]);
    gradH2 = sqrt(GradH2(:,:,:,1).^2+GradH2(:,:,:,2).^2+GradH2(:,:,:,3).^2);
    Err2   = GradH2(:,:,:,1).*vx+GradH2(:,:,:,2).*vy+GradH2(:,:,:,3).*vz;
    Err2   = abs(Err2);
    [Xv,Yv,Zv] = ndgrid(sol.x,sol.y,sol.z); % quriy points
    Err2Int    = griddedInterpolant(Xv,Yv,Zv,permute(Err2,[2,1,3]));
    GradH2Int  = griddedInterpolant(Xv,Yv,Zv,permute(gradH2,[2,1,3]));

    H2  = sol.H2;
    H2  = permute(H2,[2,1,3]);
    nlevels = 50;
    Hmin = min(H2(:));
    Hmax = max(H2(:));
    Hsamp = linspace(Hmin,Hmax,nlevels);
    threshold = 0.005;
    Hg2 = filter_H(Err2Int,GradH2Int,vnormInt,H2,Hsamp,xv,yv,zv,filename,threshold);
    set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
    xlabel('$x$','FontSize', 20, 'interpreter','latex');
    ylabel('$y$','FontSize', 20, 'interpreter','latex');
    zlabel('$z$','FontSize', 20, 'interpreter','latex');
    title(['RBC $(H_2, Ra=1\times10^5, E_A\leq0.005)$'],'FontSize', 20, 'interpreter','latex');
end
