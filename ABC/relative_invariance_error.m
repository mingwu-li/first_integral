clear all
ngrids = [5 10 20 25 29];
m = numel(ngrids);
lamd = zeros(m,3);
rsq  = zeros(m,3);
ndof = zeros(m,3);
nele = zeros(m,3);
ires = zeros(m,3);
for order = 2
    for k=1:5
        fprintf('order is %d and k is %d\n',order,k);
        tic
        ngridk = num2str(ngrids(k));
        filename = ['ABC_nx',ngridk,'ny',ngridk,'nz',ngridk,'Order',num2str(order),'PBC.mat'];
        sol = load(filename);
        lamd(k,order) = sol.lambda(2);
        ndof(k,order) = sol.numDOF;
        nele(k,order) = sol.numEle;
        ires(k,order) = averaged_relataive_invariance_residual(sol,@abc_flow);
        toc
    end
end


%% plot results
figure; 
loglog(ndof(:,1),ires(:,1),'bs-','LineWidth',2,'MarkerSize',10); hold on; 
loglog(ndof(:,2),ires(:,2),'ro-','LineWidth',2,'MarkerSize',10); grid on;
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('Number of DOFs','FontSize', 18, 'interpreter','latex');
ylabel('$E_m$','FontSize', 18, 'interpreter','latex');
title('ABC Flow','FontSize', 18, 'interpreter','latex');
figure;
loglog(nele(:,1),ires(:,1),'bs-','LineWidth',2,'MarkerSize',10); hold on
loglog(nele(:,2),ires(:,2),'ro-','LineWidth',2,'MarkerSize',10); grid on;
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('Number of Elements','FontSize', 18, 'interpreter','latex');
ylabel('$E_m$','FontSize', 18, 'interpreter','latex');
title('ABC Flow','FontSize', 18, 'interpreter','latex');


figure;
yyaxis left
loglog(nele(:,2),lamd(:,2),'o-','LineWidth',2,'MarkerSize',10); hold on; grid on
yyaxis right
loglog(nele(:,2),ires(:,2),'d-','LineWidth',2,'MarkerSize',10); 
yyaxis left
set(gca,'LineWidth',1.5); set(gca, 'FontSize', 14);
xlabel('Number of Elements','FontSize', 18, 'interpreter','latex');
ylabel('$\lambda_2$','FontSize', 18, 'interpreter','latex');
title('ABC Flow','FontSize', 18, 'interpreter','latex');
yyaxis right
ylabel('$E_m$','FontSize', 18, 'interpreter','latex');
legend('$\lambda_2$','$E_m$','interpreter','latex','FontSize',18);