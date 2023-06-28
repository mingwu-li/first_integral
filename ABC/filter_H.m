function Hg1 = filter_H(ErrInt,GradHInt,vnormInt,H,Hsamp,xv,yv,zv,filename,threshold)
% FILTER_H This function apply a filter to H and plot filtered results.
% Here the filter is based on a normalized error metric

% loop over level surfaces
Hg1 = [];
for k=1:numel(Hsamp)
    Hk = Hsamp(k);
    fk = isosurface(xv,yv,zv,H,Hk); 
    if isempty(fk.vertices); continue; end
    errk   = ErrInt(fk.vertices(:,1),fk.vertices(:,2),fk.vertices(:,3));
    gradhk = GradHInt(fk.vertices(:,1),fk.vertices(:,2),fk.vertices(:,3));
    v2k    = vnormInt(fk.vertices(:,1),fk.vertices(:,2),fk.vertices(:,3));
    errk   = errk./(gradhk.*v2k); %./(v2k)
    if sum(errk)/numel(errk)<threshold 
        Hg1 = [Hg1 Hk];
    end
end

% plot level surfaces after filter is applied
Hsamp = Hg1;
figure; 
for k=1:numel(Hsamp)
    Hk = Hsamp(k);
    isosurface(xv,yv,zv,H,Hk); hold on
end
axis equal; box on
axis([min(xv(:)) max(xv(:)) min(yv(:)) max(yv(:)) min(zv(:)) max(zv(:))]);
xlabel('x'); ylabel('y'); 
set(gca,'FontSize',14); set(gca,'LineWidth',1.5);
title([filename,'  filter on H'])

end
