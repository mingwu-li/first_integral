clear all

Ra    = '1e5';
ny    = 30;
az    = 0.2;
order = 2;
% load first integral results
bcs   = 'DPBC'; 
filename = ['LPNormRBCnx',num2str(2*ny),'ny',num2str(ny),'nz',num2str(ny),'order',num2str(order),bcs,Ra,'.mat'];
sol      = load(filename);
% load velocity profile as well
filevel = ['RB_Re',num2str(Ra),'_0092250_laplacian.mat'];
solvel  = load(filevel);
vel     = solvel.LaplacianVelNormalized;
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
nlevels = 40;
Hmin = min(H1(:));
Hmax = max(H1(:));
Hsamp = linspace(Hmin,Hmax,nlevels);
threshold = 0.062; % 0.065 0.07
Hg1 = filter_H(Err1Int,GradH1Int,vnormInt,H1,Hsamp,xv,yv,zv,filename,threshold);
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14); grid on;
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$y$','FontSize', 20, 'interpreter','latex');
zlabel('$z$','FontSize', 20, 'interpreter','latex');
title(['ARBC $(H_1, E_A\leq',num2str(threshold),')$'],'FontSize', 20, 'interpreter','latex');

if isfield(sol,'H4')
    GradH4 = permute(sol.GradH4,[2,1,3,4]);
    gradH4 = sqrt(GradH4(:,:,:,1).^2+GradH4(:,:,:,2).^2+GradH4(:,:,:,3).^2);
    Err4   = GradH4(:,:,:,1).*vx+GradH4(:,:,:,2).*vy+GradH4(:,:,:,3).*vz;
    Err4   = abs(Err4);
    [Xv,Yv,Zv] = ndgrid(sol.x,sol.y,sol.z); % quriy points
    Err4Int    = griddedInterpolant(Xv,Yv,Zv,permute(Err4,[2,1,3]));
    GradH4Int  = griddedInterpolant(Xv,Yv,Zv,permute(gradH4,[2,1,3]));

    H4  = sol.H4;
    H4  = permute(H4,[2,1,3]);
    Hmin = min(H4(:));
    Hmax = max(H4(:));
    Hsamp = linspace(Hmin,Hmax,nlevels);
    Hg4 = filter_H(Err4Int,GradH4Int,vnormInt,H4,Hsamp,xv,yv,zv,filename,threshold);
    set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14); grid on
    xlabel('$x$','FontSize', 20, 'interpreter','latex');
    ylabel('$y$','FontSize', 20, 'interpreter','latex');
    zlabel('$z$','FontSize', 20, 'interpreter','latex');
    title(['ARBC $(H_4, E_A\leq',num2str(threshold),')$'],'FontSize', 20, 'interpreter','latex');
end

%% launch streamlines via numerical integration
vxInt = griddedInterpolant(Xv,Yv,Zv,permute(vx,[2,1,3]));
vyInt = griddedInterpolant(Xv,Yv,Zv,permute(vy,[2,1,3]));
vzInt = griddedInterpolant(Xv,Yv,Zv,permute(vz,[2,1,3]));
colors = get(0,'defaultaxescolororder');
figure;
[f,v] = isosurface(xv,yv,zv,H4,Hg4(2)); hold on
h = patch('Faces',f,'Vertices',v,'FaceAlpha',0.8); grid on; box on
set(h,'FaceColor',colors(3,:),'EdgeColor','none');
axis equal; view(3);
axis([0 0.2 0 0.1 0 0.1]);
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$y$','FontSize', 20, 'interpreter','latex');
zlabel('$z$','FontSize', 20, 'interpreter','latex');
camlight; lighting gouraud
      
fk = isosurface(xv,yv,zv,H4,Hg4(2));
vk = fk.vertices;
% idxs = randi(size(vk,1),5);
% save('LPinit_1e5.mat','idxs');
load('LPinit_1e5.mat','idxs');
x0 = vk(idxs,1);
y0 = vk(idxs,2);
z0 = vk(idxs,3);
tf = 0.1;
[xt,yt,zt] = stream_lines_integration(x0,y0,z0,[],tf,{vxInt,vyInt,vzInt});
hold on;
for k=1:numel(idxs)
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

%%



%% to see wheather the data is a quasi-2d data
data = solvel;
vel  = data.LaplacianVelNormalized;
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

%% contour plots of H1 and H4 at cross section z=0
figure;
subplot(2,1,1)
[xx,yy] = meshgrid(sol.x,sol.y);
contour(xx,yy,H1(:,:,1),50); title('H1 (z=0)')
subplot(2,1,2)
contour(xx,yy,H4(:,:,1),50); title('H4 (z=0)')