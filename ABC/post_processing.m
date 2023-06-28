clear all
%% relative invariance error as functions of discretization
relative_invariance_error;


%% plot results for a given discretization
ngrids = 20; % 21;
order  = 2;
ngridk = num2str(ngrids);
filename = ['ABC_nx',ngridk,'ny',ngridk,'nz',ngridk,'Order',num2str(order),'PBC.mat'];
sol = load(filename);
plot_cross_section(sol)

% load velocity profile
[xv,yv,zv] = meshgrid(sol.x,sol.y,sol.z); % quriy points
[vx,vy,vz] = abc_flow(xv,yv,zv);
v2 = sqrt(vx.^2+vy.^2+vz.^2);
H  = sol.H;
H  = permute(H,[2,1,3]);
dx = sol.x(2)-sol.x(1);
dy = sol.y(2)-sol.y(1);
dz = sol.z(2)-sol.z(1);
GradH = permute(sol.GradH,[2,1,3,4]);
gradH = sqrt(GradH(:,:,:,1).^2+GradH(:,:,:,2).^2+GradH(:,:,:,3).^2);
err   = GradH(:,:,:,1).*vx+GradH(:,:,:,2).*vy+GradH(:,:,:,3).*vz;
% reconstruct error interpolant
[npy,npx,npz] = size(vx);
Err   = reshape(abs(err),[npy,npx,npz]);
GradH = reshape(gradH,[npy,npx,npz]);
[Xv,Yv,Zv] = ndgrid(sol.x,sol.y,sol.z); % quriy points
ErrInt     = griddedInterpolant(Xv,Yv,Zv,permute(Err,[2,1,3]));
GradHInt   = griddedInterpolant(Xv,Yv,Zv,permute(GradH,[2,1,3]));
vnormInt   = griddedInterpolant(Xv,Yv,Zv,permute(v2,[2,1,3]));
% filter H is applied
nlevels = 50;
Hmin = min(H(:));
Hmax = max(H(:));
Hsamp = linspace(Hmin,Hmax,nlevels);
threshold = 0.008;
Hg1 = filter_H(ErrInt,GradHInt,vnormInt,H,Hsamp,xv,yv,zv,filename,threshold);
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$y$','FontSize', 20, 'interpreter','latex');
zlabel('$z$','FontSize', 20, 'interpreter','latex');
title(['ABC Flow $(N=25, E_A\leq0.008)$'],'FontSize', 20, 'interpreter','latex');


% plot of sampled isosurfaces along with numerical integration trajectories
colors = get(0,'defaultaxescolororder');
figure;
[f,v] = isosurface(xv,yv,zv,H,Hg1(end-3)); hold on
h = patch('Faces',f,'Vertices',v,'FaceAlpha',0.8); grid on; box on
set(h,'FaceColor',colors(3,:),'EdgeColor','none');

fk = isosurface(xv,yv,zv,H,Hg1(end-3));
vk = fk.vertices;
[~,idx1] = min(abs(vk(:,2)-1));
[~,idx2] = min(abs(vk(:,2)-3));
[~,idx3] = min(abs(vk(:,2)-6));
idxs = [idx1 idx2 idx3];
x0 = vk(idxs,1);
y0 = vk(idxs,2);
z0 = vk(idxs,3);
tf = 70;
[xt,yt,zt] = stream_lines_integration(@abc_flow,x0,y0,z0,tf); % forward simulation
hold on
for k=1:3
    plot3(mod(xt{k},2*pi),mod(yt{k},2*pi),mod(zt{k},2*pi),'k.','LineWidth',1);
end
title(['$H=',num2str(Hg1(end-3)),'$'],'FontSize', 20, 'interpreter','latex');
view(3); grid on; camlight; lighting gouraud
