clear all

Ra    = '5e4';
ny    = 30;
az    = 0.2;
order = 2;
% load first integral results
bcs      = 'DPBC'; 
filename = ['LPNormRBCnx',num2str(2*ny),'ny',num2str(ny),'nz',num2str(ny),'order',num2str(order),bcs,Ra,'.mat'];
sol      = load(filename);
% load velocity profile as well
filevel = ['RB_Re',num2str(Ra),'_laplacian.mat'];
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
Hmin  = min(H1(:));
Hmax  = max(H1(:));
Hsamp = linspace(Hmin,Hmax,nlevels);
threshold = 0.01; % 0.03 or 0.05
Hg1 = filter_H(Err1Int,GradH1Int,vnormInt,H1,Hsamp,xv,yv,zv,filename,threshold);
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14); grid on;
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$y$','FontSize', 20, 'interpreter','latex');
zlabel('$z$','FontSize', 20, 'interpreter','latex');
title(['ARBC $(H_1, E_A\leq',num2str(threshold),')$'],'FontSize', 20, 'interpreter','latex');

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
    Hmin = min(H2(:));
    Hmax = max(H2(:));
    Hsamp = linspace(Hmin,Hmax,nlevels);
    Hg2 = filter_H(Err2Int,GradH2Int,vnormInt,H2,Hsamp,xv,yv,zv,filename,threshold);
    set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14); grid on
    xlabel('$x$','FontSize', 20, 'interpreter','latex');
    ylabel('$y$','FontSize', 20, 'interpreter','latex');
    zlabel('$z$','FontSize', 20, 'interpreter','latex');
    title(['ARBC $(H_2, E_A\leq',num2str(threshold),')$'],'FontSize', 20, 'interpreter','latex');
end

%% launch streamlines via numerical integration
vxInt = griddedInterpolant(Xv,Yv,Zv,permute(vx,[2,1,3]));
vyInt = griddedInterpolant(Xv,Yv,Zv,permute(vy,[2,1,3]));
vzInt = griddedInterpolant(Xv,Yv,Zv,permute(vz,[2,1,3]));
colors = get(0,'defaultaxescolororder');
figure;
[f,v] = isosurface(xv,yv,zv,H1,Hg1(12)); hold on
h = patch('Faces',f,'Vertices',v,'FaceAlpha',0.8); grid on; box on
set(h,'FaceColor',colors(3,:),'EdgeColor','none');
axis equal; view(3);
axis([0 0.2 0 0.1 0 0.1]);
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$y$','FontSize', 20, 'interpreter','latex');
zlabel('$z$','FontSize', 20, 'interpreter','latex');
camlight; lighting gouraud
      
fk = isosurface(xv,yv,zv,H1,Hg1(12));
vk = fk.vertices;
[~,idx1] = min(abs(vk(:,3)-0.02)+abs(vk(:,1)-0.15));
[~,idx2] = min(abs(vk(:,3)-0.04)+abs(vk(:,1)-0.15));
[~,idx3] = min(abs(vk(:,3)-0.06)+abs(vk(:,1)-0.15));
[~,idx4] = min(abs(vk(:,3)-0.08)+abs(vk(:,1)-0.15));
idxs = [idx1 idx2 idx3 idx4];
x0 = vk(idxs,1);
y0 = vk(idxs,2);
z0 = vk(idxs,3);
% tf = 2;
tf = 50;
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
axis([0.1 0.2 0 0.1 0 0.1]);

% poincare section
pfun = @(x,y,z) y-0.05;
figure; hold on
for k=1:4
    [xp,yp,zp] = extract_poincare_section(xt{k},yt{k},zt{k},pfun);
    plot(xp,zp,'k.','MarkerSize',8);
end
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 20, 'interpreter','latex'); grid on; box on;
ylabel('$z$','FontSize', 20, 'interpreter','latex');


%% repeat the above procedure for the other surface
figure;
[f,v] = isosurface(xv,yv,zv,H2,Hg2(1)); hold on
h = patch('Faces',f,'Vertices',v,'FaceAlpha',0.8); grid on; box on
set(h,'FaceColor',colors(3,:),'EdgeColor','none');
axis equal; view(3);
axis([0 0.2 0 0.1 0 0.1]);
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$y$','FontSize', 20, 'interpreter','latex');
zlabel('$z$','FontSize', 20, 'interpreter','latex');
camlight; lighting gouraud
      
fk = isosurface(xv,yv,zv,H2,Hg2(1));
vk = fk.vertices;
[~,idx1] = min(abs(vk(:,3)-0.02)+abs(vk(:,1)-0.15));
[~,idx2] = min(abs(vk(:,3)-0.04)+abs(vk(:,1)-0.15));
[~,idx3] = min(abs(vk(:,3)-0.06)+abs(vk(:,1)-0.15));
[~,idx4] = min(abs(vk(:,3)-0.08)+abs(vk(:,1)-0.15));
idxs = [idx1 idx2 idx3 idx4];
x0 = vk(idxs,1);
y0 = vk(idxs,2);
z0 = vk(idxs,3);
tf = 50;
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
axis([0.0 0.1 0 0.1 0 0.1]);
title('$H_2$','FontSize', 20, 'interpreter','latex');

% Poincare section
figure; hold on
for k=1:4
    [xp,yp,zp] = extract_poincare_section(xt{k},yt{k},zt{k},pfun);
    plot(xp,zp,'k.','MarkerSize',8);
end
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 20, 'interpreter','latex'); grid on; box on;
ylabel('$z$','FontSize', 20, 'interpreter','latex');

