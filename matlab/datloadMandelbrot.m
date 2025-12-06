A = importdata('Mandlebrot_X0=0.310000_0.121100.dat',' ');
X0 = [0.310000,0.121100];
NT = 1e6;

% A = importdata('Mandlebrot_X0=0.336960_0.337970.dat',' ');
% NT = 1e5;

Ps = A(:,1);
kappas = A(:,2);
maxIters = A(:,3);

figure
% subplot(2,1,1)
% scatter(Ps(maxIters~=NT),kappas(maxIters~=NT),8,'filled')
% hold on
% scatter(Ps(maxIters==NT),kappas(maxIters==NT),8,'filled')
% ylim([3,3.7])
% title(strcat('$\vec{x}_0=(',num2str(X0(1)),',',num2str(X0(2)),')$'),'interpreter','latex','fontsize',12)
% ylabel('$\kappa$','interpreter','latex','fontsize',12)
% subplot(2,1,2)
% 
scatter(Ps(maxIters~=NT),kappas(maxIters~=NT),8,log10(maxIters(maxIters~=NT)/NT),'filled')
hold on
scatter(Ps(maxIters==NT),kappas(maxIters==NT),8,'filled')
ylim([3,3.7])
xlabel('$p_{\sf s}$','interpreter','latex','fontsize',12)
ylabel('$\kappa$','interpreter','latex','fontsize',12)
title(strcat('$\vec{x}_0=(',num2str(X0(1)),',',num2str(X0(2)),')$'),'interpreter','latex','fontsize',12)
h = colorbar;

%% Export admissible data, Import Lyapunov data
A = importdata('Mandlebrot_X0=0.310000_0.121100.dat',' ');
X0 = [0.310000,0.121100];
NT = 1e6;

Ps = A(:,1);
kappas = A(:,2);
maxIters = A(:,3);

%Step 1: export admissible data (Ps,kappa) pairs
i = find(maxIters==NT);
ikappa = kappas(i);
iPs = Ps(i);

LyapArgs = horzcat(iPs,ikappa);
save('LyapArgs.dat','LyapArgs','-ascii')

%Step 2: run the bash script MandelLyap.sh to calculate Lyapunovs
%Step 3: import data
filenames = dir('output_SingleX0=0.310000_0.121100_k=*.dat');
filenames = {filenames.name};

nFiles = numel(filenames);
mles = zeros(nFiles,1);
cumulants = zeros(nFiles,1);

for j=1:nFiles
    A = importdata(filenames{j},' ');
    mles(j) = A(3);
    cumulants(j) = A(4);
end

save('output_SingleX0=0.310000_0.121100_RandKappaPs','ikappa','iPs','mles','cumulants','NT')

iMLES = mles;
iMLES(iMLES<0) = -1;
iMLES(iMLES>0) = +1;

[ikappaSort,i] = sort(ikappa); %since dir loads filenames sorted by kappa
iPsSort = iPs(i);

%% 2x2 subplots:
figure
subplot(2,2,1)
scatter(iPsSort(mles<0),ikappaSort(mles<0),8,mles(mles<0),'filled')
xlabel('$p_{\sf s}$','FontSize',12,'Interpreter','latex')
ylabel('$\kappa$','FontSize',12,'Interpreter','latex')
h = colorbar;
ylabel(h,'$\lambda$','FontSize',12,'Interpreter','latex')

subplot(2,2,2)
scatter(iPsSort(mles<0),ikappaSort(mles<0),8,log10(cumulants(mles<0)),'filled')
xlabel('$p_{\sf s}$','FontSize',12,'Interpreter','latex')
ylabel('$\kappa$','FontSize',12,'Interpreter','latex')
h = colorbar;
ylabel(h,'$\log_{10}\left(\lambda^{(4)}-3(\lambda^{(2)})^2\right)$','FontSize',12,'Interpreter','latex')

subplot(2,2,3)
scatter(iPsSort(mles>0),ikappaSort(mles>0),8,mles(mles>0),'filled')
xlabel('$p_{\sf s}$','FontSize',12,'Interpreter','latex')
ylabel('$\kappa$','FontSize',12,'Interpreter','latex')
h = colorbar;
ylabel(h,'$\lambda$','FontSize',12,'Interpreter','latex')

subplot(2,2,4)
scatter(iPsSort(mles>0),ikappaSort(mles>0),8,log10(cumulants(mles>0)),'filled')
xlabel('$p_{\sf s}$','FontSize',12,'Interpreter','latex')
ylabel('$\kappa$','FontSize',12,'Interpreter','latex')
h = colorbar;
ylabel(h,'$\log_{10}\left(\lambda^{(4)}-3(\lambda^{(2)})^2\right)$','FontSize',12,'Interpreter','latex')

%log10 mles:
figure
subplot(2,2,1)
scatter(iPsSort(mles<0),ikappaSort(mles<0),8,log10(-mles(mles<0)),'filled')

subplot(2,2,2)
scatter(iPsSort(mles<0),ikappaSort(mles<0),8,log10(cumulants(mles<0)),'filled')

subplot(2,2,3)
scatter(iPsSort(mles>0),ikappaSort(mles>0),8,log10(mles(mles>0)),'filled')

subplot(2,2,4)
scatter(iPsSort(mles>0),ikappaSort(mles>0),8,log10(cumulants(mles>0)),'filled')






figure
scatter(iPsSort,ikappaSort,8,mles,'filled')
ylabel('$\kappa$','FontSize',12,'Interpreter','latex')
xlabel('$p_{\sf s}$','FontSize',12,'Interpreter','latex')
title('$\vec{x}_0=(0.3100,0.1211)$','FontSize',12,'Interpreter','latex')
h = colorbar;
ylabel(h,'$\lambda$','FontSize',12,'Interpreter','latex')

figure
scatter(iPsSort,ikappaSort,8,iMLES,'filled')
ylabel('$\kappa$','FontSize',12,'Interpreter','latex')
xlabel('$p_{\sf s}$','FontSize',12,'Interpreter','latex')
title('$\vec{x}_0=(0.3100,0.1211)$','FontSize',12,'Interpreter','latex')
h = colorbar('XTickLabel',{'-1','0','+1'},'XTick', -1:1:1);
ylabel(h,'${\rm sgn}(\lambda)$','FontSize',12,'Interpreter','latex')

jMLES = mles;
jMLES(iMLES<0) = 0;
jMLES(iMLES>0) = 1;

figure
MarkerColor = [62 38 168]/256;
scatter(iPsSort(jMLES==0),ikappaSort(jMLES==0),8,MarkerColor,'filled')
hold on
MarkerColor = [249 251 20]/256;
scatter(iPsSort(jMLES>0),ikappaSort(jMLES>0),8,MarkerColor,'filled')
ylabel('$\kappa$','FontSize',12,'Interpreter','latex')
xlabel('$p_{\sf s}$','FontSize',12,'Interpreter','latex')
title('$\vec{x}_0=(0.3100,0.1211)$','FontSize',12,'Interpreter','latex')
legend('$\lambda\leq0$','$\lambda>0$')



figure
scatter(iPsSort,ikappaSort,8,log10(cumulants),'filled')
ylabel('$\kappa$','FontSize',12,'Interpreter','latex')
xlabel('$p_{\sf s}$','FontSize',12,'Interpreter','latex')
title('$\vec{x}_0=(0.3100,0.1211)$','FontSize',12,'Interpreter','latex')
h = colorbar;
ylabel(h,'$\log_{10}\left(\lambda^{(4)}-3(\lambda^{(2)})^2\right)$','FontSize',12,'Interpreter','latex')