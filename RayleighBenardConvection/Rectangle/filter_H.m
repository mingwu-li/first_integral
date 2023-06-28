function [Hg1,varargout] = filter_H(ErrInt,GradHInt,vnormInt,H,Hsamp,xv,yv,zv,filename,threshold,varargin)


Hg1 = []; Hg2 = [];
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
        if numel(varargin)>0 && strcmp(varargin{1},'wall')
            % compute the fraction of points near the two walls
            nbcpts = sum(abs(fk.vertices(:,2))<typethreshold(1)*0.1)+...
                sum(abs(fk.vertices(:,2))>0.1-typethreshold(1)*0.1);
            % filter the contour surface if the fraction is larger than
            % the typethreshold
            if nbcpts/numel(errk)<=typethreshold(2)
                Hg2 = [Hg2 Hk];
            end
        end
            
    end
end
varargout{1} = Hg2;

Hsamp = Hg1;
figure; 
for k=1:numel(Hsamp)
    Hk = Hsamp(k);
    isosurface(xv,yv,zv,H,Hk); hold on
end
axis equal; box on
axis([0 0.2 0 0.1 0 0.1]);
xlabel('x'); ylabel('y'); 
set(gca,'FontSize',14); set(gca,'LineWidth',1.5);
title([filename,'  filter on H'])

if numel(Hg2)>0
    Hsamp = Hg2;
    figure; 
    for k=1:numel(Hsamp)
        Hk = Hsamp(k);
        isosurface(xv,yv,zv,H,Hk); hold on
    end
    axis equal; box on
    axis([0 0.2 0 0.1 0 0.1]);
    xlabel('x'); ylabel('y'); 
    set(gca,'FontSize',14); set(gca,'LineWidth',1.5);
    title([filename,'  filter on H and v']);
end

end
