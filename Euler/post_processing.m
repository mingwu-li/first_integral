clear all

%% plot leading eigenvalue as a function of number of finite elements
% load results
ngrids = [10 20 25 29 35];
m = numel(ngrids);
lamd = zeros(m,1);
rsq  = zeros(m,1);
ndof = zeros(m,1);
nele = zeros(m,1);
ires = zeros(m,1);
order= 2;
for k=1:m
    fprintf('order is %d and k is %d\n',order,k);
    tic
    ngridk = num2str(ngrids(k));
    filename = ['Euler_nx',ngridk,'ny',ngridk,'nz',ngridk,'Order',num2str(order),'PBC.mat'];
    sol = load(filename);
    lamd(k) = sol.lambda(2);
    ndof(k) = sol.numDOF;
    nele(k) = sol.numEle;
    ires(k) = averaged_relataive_invariance_residual(sol);
    toc
end
% plot results
figure;
yyaxis left;
loglog(nele,lamd,'o-','LineWidth',2,'MarkerSize',10); hold on; grid on
yyaxis right
loglog(nele,ires,'d-','LineWidth',2,'MarkerSize',10); 
yyaxis left
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('Number of Elements','FontSize', 18, 'interpreter','latex');
ylabel('$\lambda_2$','FontSize', 18, 'interpreter','latex');
title('ABC Flow','FontSize', 18, 'interpreter','latex');
yyaxis right
ylabel('$E_m$','FontSize', 18, 'interpreter','latex');
legend('$\lambda_2$','$E_m$','interpreter','latex','FontSize',18);


%% plot isosurfaces and level surfaces
ngrids = 35; % 21;
order  = 2;
ngridk = num2str(ngrids);
filename = ['Euler_nx',ngridk,'ny',ngridk,'nz',ngridk,'Order',num2str(order),'PBC.mat'];
sol = load(filename);
plot_cross_section(sol) 
title('$N=25,\mathcal{O}(2)$','FontSize', 20, 'interpreter','latex');

% load velocity profile as well
[xv,yv,zv] = meshgrid(sol.x,sol.y,sol.z); % quriy points
[vx,vy,vz] = euler_flow(xv,yv,zv);
v2 = sqrt(vx.^2+vy.^2+vz.^2);
v2(1,1,1) = eps;
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
% plot isosurfaces
nlevels = 50;
Hmin = min(H(:));
Hmax = max(H(:));
Hsamp = linspace(Hmin,Hmax,nlevels);
threshold = 0.01;Hg1 = filter_H(ErrInt,GradHInt,vnormInt,H,Hsamp,xv,yv,zv,filename,threshold);
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('$x$','FontSize', 20, 'interpreter','latex');
ylabel('$y$','FontSize', 20, 'interpreter','latex');
zlabel('$z$','FontSize', 20, 'interpreter','latex');
title(['Euler Flow $(N=35, E_A\leq0.01)$'],'FontSize', 20, 'interpreter','latex');
colors = get(0,'defaultaxescolororder');
% plot single isosurface along with simulated trajectories
figure;
% hold on
[f,v] = isosurface(xv,yv,zv,H,Hg1(1)); hold on
h = patch('Faces',f,'Vertices',v,'FaceAlpha',0.8); grid on; box on
set(h,'FaceColor',colors(2,:),'EdgeColor','none');

fk = isosurface(xv,yv,zv,H,Hg1(1));
vk = fk.vertices;
npts = size(vk,1);
idxs = randi(npts,100,1);
x0 = vk(idxs,1);
y0 = vk(idxs,2);
z0 = vk(idxs,3);
tf = 100;
[xt,yt,zt] = stream_lines_integration(x0,y0,z0,tf);
hold on
for k=1:5
    plot3(mod(xt{k},2*pi),mod(yt{k},2*pi),mod(zt{k},2*pi),'k.','LineWidth',1);
end
xlim([0 2*pi]); ylim([0 2*pi]); zlim([0. 2*pi]);
view(3); grid on; camlight; lighting gouraud


% figure;
[f,v] = isosurface(xv,yv,zv,H,Hg1(end)); hold on
h = patch('Faces',f,'Vertices',v,'FaceAlpha',0.8); grid on; box on
set(h,'FaceColor',colors(5,:),'EdgeColor','none');
fk = isosurface(xv,yv,zv,H,Hg1(end));
vk = fk.vertices;
[~,idx1] = min(abs(vk(:,3)-1));
[~,idx2] = min(abs(vk(:,3)-3));
[~,idx3] = min(abs(vk(:,3)-6));
idxs = [idx1 idx2 idx3];
x0 = vk(idxs,1);
y0 = vk(idxs,2);
z0 = vk(idxs,3);
tf = 50;
[xt,yt,zt] = stream_lines_integration(x0,y0,z0,tf);
hold on
for k=1:3
    plot3(mod(xt{k},2*pi),mod(yt{k},2*pi),mod(zt{k},2*pi),'k.','LineWidth',1);
end
xlim([0 2*pi]); ylim([0 2*pi]); zlim([0. 2*pi]);
view(3); grid on; camlight; lighting gouraud
