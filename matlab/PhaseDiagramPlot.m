%
fname = 'ChaoticWindow_k=3.5703';
openfig([fname,'.fig'])

fig = gcf;
fig.PaperPositionMode = 'auto';
print(gcf,[fname,'.png'],'-dpng','-r300','-painters');




NT = 1e6;
nSamples = 1000;

A = importdata('bifurc_X0=0.310000_0.121100_Ps=0.320000.dat',' ');

kappas = A.textdata(1:end,1);
N = numel(kappas);

for i=1:N
    kappas{i} = str2double(kappas{i});
end

kappas = cell2mat(kappas);

%preparation:
A = A.data;
B = zeros(numel(A),2);

for n=1:N
    B((n-1)*nSamples+1:n*nSamples,1) = kappas(n);
    B((n-1)*nSamples+1:n*nSamples,2) = A(n,:);
end

figure
scatter(B(:,1),B(:,2),4,'filled')

%%
%kc = [1 0 0 0 0 0 0];
kp = [3,3,3,3.0001,3.0007,3.0064,3.0718];%,3.1162];
Ps2 = [1e-5 3.2e-5 1e-4 3.2e-4 1e-3 1e-2 0.1];% 0.15];

%k1 = [0,0,0,3,3,3,0];
k2 = [3.6740,3.6739,3.6737,3.6729,3.6705,3.6629,3.6406,3.5772,3.4354,3.3597];
Ps3 = [1e-5 3.2e-5 1e-4 3.2e-4 1e-3 3.2e-3 1e-2 3.2e-2 0.1 0.15];

kmax = [3.6761,3.6760,3.6757,3.6750,3.6726,3.6651,3.6428,3.5806,3.4474,3.3745,3.168,3.3078];
Ps1 = [1e-5 3.2e-5 1e-4 3.2e-4 1e-3 3.2e-3 1e-2 3.2e-2 0.1 0.15 0.2 0.32];

kJoin = [kp,flip(k2)];
PsJoin = [Ps2,flip(Ps3)];

figure
fill([1e-6 0.33 0.33 1e-6],[3 3 2.95 2.95],'g')
hold on
fill(PsJoin,kJoin,'b')
hold on
fill([Ps1,0.33,0.33,1e-3,1e-2,0.1,flip(Ps3)],[kmax,3.4,3,3,3.0064,3.0718,flip(k2)],'r')

xlabel('$p_{\sf s}$','interpreter','latex','fontsize',12)
ylabel('$\kappa$','interpreter','latex','fontsize',12)

%% Using LyapStats

A = importdata('output_SingleX0=0.310000_0.121100_Ps=0.000000.dat',' ');
NT = 1e6;

kappas = A(:,1);
mles = A(:,3);
cumulants = A(:,4);

% figure
hold on
plot(kappas,mles)
hold on;plot([min(kappas),max(kappas)],[0,0],'k--')
ylabel(strcat('$\lambda(\kappa,p_{\sf s}=0)$'),'FontSize',12,'Interpreter','latex')
xlabel('$\kappa$','FontSize',12,'Interpreter','latex')
kappaPlus = kappas(mles>0);
kappaMax = kappas(find(mles==0,1,'first'))

%title('$X_0=(0.3100,0.1211)$','FontSize',12,'Interpreter','latex')

hold on;plot(kappas,movmean(mles,200))


figure
plot(kappas,cumulants)