
NT = 1e6;

A = importdata('bifurc_kappa=3.667300.dat',' ');
kappa = 3.6673;

nSamples = size(A.data,2);
M = size(A.data,1);

X0 = zeros(M,2);
nUniques = zeros(M,1);
for m=1:M
    X0(m,1) = str2double(A.textdata(m,1));
    X0(m,2) = str2double(A.textdata(m,2));
    
    nUniques(m) = numel(unique(A.data(m,:)));
end

Cycles = unique(nUniques);

% figure
subplot(2,2,4)

for n=1:numel(Cycles)
    hold on
    scatter(X0(nUniques==Cycles(n),1),X0(nUniques==Cycles(n),2),8,'filled')
end

xlabel('$x_{1,0}$','FontSize',12,'Interpreter','latex')
ylabel('$x_{2,0}$','FontSize',12,'Interpreter','latex')
legend({'Period-16','Period-24'},'FontSize',12,'Interpreter','latex')
title('$\kappa = 3.6673$, $p_{\sf s}=0$','FontSize',12,'Interpreter','latex')