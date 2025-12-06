A = importdata('output_SingleX0=0.178000_0.627700_Ps=0.dat',' ');

NT = 1e6;

kappas = A(:,1);
mles = A(:,3);
cumulants = A(:,4);

figure
hold on;plot(3.6276*[1,1],[-0.1,0.06],'k-')
hold on;plot(3.6666*[1,1],[-0.1,0.06],'k-')
hold on;plot(3.6725*[1,1],[-0.1,0.06],'k-')
%Zoom lines:
% hold on;plot(3.67252*[1,1],[-0.1,0.06],'k-')
% hold on;plot(3.67367*[1,1],[-0.1,0.06],'k-')
% hold on;plot(3.67392*[1,1],[-0.1,0.06],'k-')

yyaxis left;
plot(kappas,mles)
hold on;plot([min(kappas),max(kappas)],[0,0],'k--')
ylabel(strcat('$\lambda(\kappa,p_{\sf s}=0)$'),'FontSize',12,'Interpreter','latex')
xlabel('$\kappa$','FontSize',12,'Interpreter','latex')
%title('$X_0=(0.3100,0.1211)$','FontSize',12,'Interpreter','latex')
set(gca,'ycolor','k') 

yyaxis right;
plot(kappas,cumulants)
ylabel('$\lambda^{(4)}-3(\lambda^{(2)})^2$','FontSize',12,'Interpreter','latex')
set(gca,'ycolor','r') 



% rectangle('Position', [3.6276, -0.1, 0.0390, 0.16], ...
%                 'FaceColor', [0.23, 0.44, 0.34, 0.4], ...
%                 'EdgeColor', [0, 0, 0, 0]);

% rectangle('Position', [3.6, -0.1, 0.0276, 0.16], ...
%                 'FaceColor', [0.39, 0.47, 0.64, 0.4], ...
%                 'EdgeColor', [0, 0, 0, 0]);      
% rectangle('Position', [3.6737, -0.1, 0.0002, 0.16], ...
%                 'FaceColor', [0.51, 0.38, 0.48, 0.4], ...
%                 'EdgeColor', [0, 0, 0, 0]);

