%% Scatter plot of ICs:
A = importdata('ScanICs_FixedKappa3.400000.dat',' ');
X1 = A(:,1);
X2 = A(:,2);
T = A(:,3);

figure
subplot(1,3,3)
scatter(X1(T==1e6),X2(T==1e6),8,'filled')
xlim([0,1])
ylim([0,1])
ylabel('$x_{2,0}$','FontSize',12,'Interpreter','latex')
xlabel('$x_{1,0}$','FontSize',12,'Interpreter','latex')
title('$\kappa=3.40$','FontSize',12,'Interpreter','latex')

%%
Ps = 0;
A = importdata('output_SingleX0=0.310000_0.121100.dat',' ');
% kappas = 3.625:1e-6:(3.6761-1e-6);
kappas = 3.6275:1e-8:(3.6277);

t = A.data(:,1);
% tmax = find(t~=NT,1,'first')-1;
tmax = numel(t);

mles = A.data(:,2);
cumulants = A.data(:,3);

figure
yyaxis left;
plot(kappas(1:tmax),mles(1:tmax))
hold on;plot([min(kappas(1:tmax)),max(kappas(1:tmax))],[0,0],'k--')
ylabel(strcat('$\lambda(\kappa,p_{\sf s}=',num2str(Ps),')$'),'FontSize',12,'Interpreter','latex')
xlabel('$\kappa$','FontSize',12,'Interpreter','latex')
title('$X_0=(0.3100,0.1211)$','FontSize',12,'Interpreter','latex')

yyaxis right;
plot(kappas(1:tmax),cumulants(1:tmax))
ylabel('$3(\lambda^{(2)})^2-\lambda^{(4)}$','FontSize',12,'Interpreter','latex')


figure
yyaxis left;
plot(kappas(1:7022),SingleX0.mles(1:7022))
hold on;plot([min(kappas(1:7022)),max(kappas(1:7022))],[0,0],'k--')
ylabel(strcat('$\lambda(\kappa,p_{\sf s}=',num2str(Ps),')$'),'FontSize',12,'Interpreter','latex')
xlabel('$\kappa$','FontSize',12,'Interpreter','latex')

yyaxis right;
plot(kappas(1:7022),SingleX0.cumulants(1:7022))
ylabel('$3(\lambda^{(2)})^2-\lambda^{(4)}$','FontSize',12,'Interpreter','latex')

%%
%initialization:
NT = 1e6;
nSamples = 1e4;

A = importdata('data (1)/output_0.dat',' ',1);

kappas = A.textdata(2:end);
N = numel(kappas);

for i=1:N
    kappas{i} = str2double(kappas{i}(1:end-1));
end

kappas = cell2mat(kappas);

%analysis:
kappa = 3;
m = find(kappas==kappa);

X0 = zeros(nSamples,2);
maxKappa = zeros(nSamples,1);
admisX0 = [];
inadmisX0 = [];
mles = [];
t = [];

%~200s loop:
tic
for n=0:nSamples-1
    filename = ['data (1)/output_',num2str(n),'.dat'];

    temp = importdata(filename);
    A = importdata(filename,' ',1);
    
    maxIters = A.data(:,1);
    maxKappa(n+1) = find(A.data(:,1)~=NT,1,'first');
    X0(n+1,:) = temp(1,:);
    
    if maxIters(m)==NT
        admisX0 = vertcat(admisX0,temp(1,:));
        mles = vertcat(mles,A.data(m,4));
    else
        t = vertcat(t,maxIters(m)/NT);
        inadmisX0 = vertcat(inadmisX0,temp(1,:));
    	continue
    end
end
toc

figure
scatter(X0(:,1),X0(:,2),8,kappas(maxKappa),'filled')
xlabel('$x_{1,0}$','FontSize',12,'Interpreter','latex')
ylabel('$x_{2,0}$','FontSize',12,'Interpreter','latex')
hold on
plot([0,1],[1,0],'k-')
hold on
plot([0,1],[0,1],'k-')
h = colorbar;
set(get(h,'label'),'string','$\kappa^*$','FontSize',12,'Interpreter','latex');

index = find(maxKappa==max(maxKappa))-1;
X0(index,:)
kappas(max(maxKappa)-1)

filename = ['data/output_',num2str(index),'.dat'];
temp = importdata(filename);
A = importdata(filename,' ',1);
    
    

figure
scatter3(kappas(maxKappa),X0(:,1),X0(:,2),8,kappas(maxKappa),'filled')


figure
scatter(admisX0(:,1),admisX0(:,2),8,log10(mles),'filled')

figure
scatter(admisX0(:,1),admisX0(:,2),8,'filled')
hold on
scatter(inadmisX0(:,1),inadmisX0(:,2),8,'filled')

figure
scatter(inadmisX0(:,1),inadmisX0(:,2),8,log10(t),'filled')


n = find(maxKappa==max(maxKappa))-1;
kmax = kappas(max(maxKappa));
filename = ['data/output_',num2str(n),'.dat'];
A = importdata(filename,' ',1);

maxIters = A.data(:,1);
mles = A.data(:,4);
cumulants = A.data(:,5);

figure
plot(kappas(mles~=0),mles(mles~=0))
hold on
plot([kmax,kmax],[1e-25,1],'k--')

figure
semilogy(kappas(mles~=0),abs(cumulants(mles~=0)))

figure
scatter(A.data(:,2),A.data(:,3),8,A.data(:,1)/NT,'filled')

figure
scatter(A.data(:,2),A.data(:,3),8,kappas,'filled')

figure
scatter3(kappas,A.data(:,2),A.data(:,3),8,A.data(:,1)/NT,'filled')

%%
n=5047;
    filename = ['data (1)/output_',num2str(n),'.dat'];

    temp = importdata(filename);
    A = importdata(filename,' ',1);
    
    maxIters = A.data(:,1);
    maxKappa(n+1) = find(A.data(:,1)~=NT,1,'first');
    X0(n+1,:) = temp(1,:);
    
    if maxIters(m)==NT
        admisX0 = vertcat(admisX0,temp(1,:));
        mles = vertcat(mles,A.data(m,4));
    else
        t = vertcat(t,maxIters(m)/NT);
        inadmisX0 = vertcat(inadmisX0,temp(1,:));
    	continue
    end
    