
%% %initialization:
NT = 1e6;
nSamples = 100;

% A = importdata('bifurc_multiX0 (1e-4).dat',' ');
% A = importdata('bifurc_X0=0.310000_0.121100 (1).dat',' ');
A = importdata('bifurc_X0=0.178000_0.627700_Ps=0.dat',' ');
% A = importdata('bifurc_X0=0.300000_0.300000.dat',' ');
% A = importdata('bifurc_X0=0.500000_0.100000.dat',' ');

kappas = A.textdata(1:end,1);
N = numel(kappas);

for i=1:N
    kappas{i} = str2double(kappas{i});
%     kappas{i} = str2double(kappas{i}(1:end-1));
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
% hold on
scatter(B(:,1),B(:,2),4,'filled')
xlabel('$\kappa$','FontSize',12,'Interpreter','latex')
ylabel('$\lim_{n\to\infty}x_{1,n}$','FontSize',12,'Interpreter','latex')
% title('$\vec{x}_0=(0.3100,0.1211)$','FontSize',12,'Interpreter','latex')
xlim([2.95,3.6761])


nUniques = zeros(size(kappas));
for i=1:size(A,1)
    nUniques(i) = numel(unique(A(i,:)));
end
hold on
plot(kappas,nUniques)
%% Branching diagram (Baushev et Zhusubaliyev 1992):
delta = 3;
C = horzcat(round(B(:,1),delta),B(:,2));
kappas2 = unique(C(:,1));

M = numel(kappas2);
D = zeros(M,1);

tic
for m=1:M
    D(m) = numel(unique(round(C(C(:,1) == kappas2(m),2),delta)));
end
toc

hold on
plot(kappas2,D)

%% By hand:
% kappaRangeA = [2.95,3,3,3.62762,3.62762,3.66659,3.66659,3.67252,3.67252,3.67367,3.67367,3.67392];
% periodA = [1,1,4,4,8,8,16,16,32,32,64,64];
kappaRangeA = [3,3.62762,3.62762,3.66659,3.66659,3.67252,3.67252,3.67367,3.67367,3.67392];
periodA = [4,4,8,8,16,16,32,32,64,64];

figure
plot(kappaRangeA,periodA)
hold on
plot([min(kappaRangeA),max(kappaRangeA)],[min(periodA),max(periodA)],'.')
xlabel('$\kappa$','FontSize',12,'Interpreter','latex')
ylabel('Cycle period','FontSize',12,'Interpreter','latex')

kappaRangeB = [3.20196,3.21055,3.21055,3.21424,3.21424,3.2149];
periodB = [5,5,10,10,20,20];

hold on
plot(kappaRangeB,periodB)
hold on
plot([min(kappaRangeB),max(kappaRangeB)],[min(periodB),max(periodB)],'.')

kappaRangeC = [3.55624,3.56694,3.56694,3.56953,3.56953,3.57013];
periodC = [12,12,24,24,48,48];

hold on
plot(kappaRangeC,periodC)
hold on
plot([min(kappaRangeC),max(kappaRangeC)],[min(periodC),max(periodC)],'.')

kappaRangeD = [3.51209,3.51385,3.51385,3.51464,3.51464];
periodD = [15,15,30,30,60];

hold on
plot(kappaRangeD,periodD)
hold on
plot([min(kappaRangeD),max(kappaRangeD)],[min(periodD),max(periodD)],'.')

kappaRangeE = [3.66726,3.6678,3.6678,3.66832,3.66832];
periodE = [17,17,34,34,68];
hold on
plot(kappaRangeE,periodE)
hold on
plot([min(kappaRangeE),max(kappaRangeE)],[min(periodE),max(periodE)],'.')















