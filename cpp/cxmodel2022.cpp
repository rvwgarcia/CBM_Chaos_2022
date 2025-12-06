/*
	Produces avalanche statistics.
	MATLAB: 512 nodes, 1e5 timesteps ~=140s

	Updates to come: Output in c-web format.

	\rvwg2022
*/

#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <tuple>
#include <math.h>
#include <fstream>

#include <omp.h>
#include <unistd.h>

using namespace std;
using namespace std::chrono;

const int taur = 2;
double Ps, NT;

double * icGen();
tuple<int, double *, double, double> kin1cbm(double kappa, double X0[2]);

//driver function
int main() {
	return 0;
}


int N = 128;
int kin = 5;
double kappa = 1.12;
double B = 1.4;
int Taur = 1;
double Ps = 1e-4;

NT = 1e5;   %or number of c-webs

WeightMatrix = zeros(N);
WeightKernel = exp(-B*(1:kin));
KernelSum = sum(WeightKernel);

for n=1:N
    Inputs = 1:N;
    Inputs = find(Inputs~=n);
    Inputs = Inputs(randi(N-1,kin,1));
    %WeightMatrix(n,Inputs) = rand(kin,1)/10;
    
    WeightMatrix(n,Inputs) = kappa*WeightKernel/KernelSum;
end

StateVector = zeros(N,1);
StateVectorNext = zeros(N,1);
StateMat = [];

tic
for t=1:NT-1
    Actives = find(StateVector==1);
    
    if numel(Actives)~=0
        StateMat = vertcat(StateMat,[t*ones(size(Actives)),Actives]);
    end
    
    for n=1:N
        temp = WeightMatrix(n,Actives(Actives~=n));
        if numel(temp)==0
            temp=0;
        else
            temp = temp-rand(size(temp));
        end
        
        StateVectorNext(n) = myStep(Ps-rand);
        StateVectorNext(n) = mod(StateVector(n)+1,Taur+1)*...
            myStep(myStep(StateVector(n))+myStep(Ps-rand)+...
            sum(myStep(temp)));
    end
    StateVector = StateVectorNext;
    if mod(t,1e4)==0
        disp([num2str(t),' iterations completed.'])
    end
end
toc

mean_firing_rate = size(StateMat,1)/(N*NT); %average firing rate <x_1>
disp(['Mean firing rate, x_1 = ', num2str(mean_firing_rate),' activations/(neuron*iteration).'])
%% Raster plot:

% figure
% scatter(StateMat(:,1),StateMat(:,2),8,'filled')

%Shows the population density plot:
uTimes = unique(StateMat(:,1));
counts = zeros(size(uTimes));
for t=1:numel(uTimes)
    counts(t) = sum(StateMat(:,1)==uTimes(t));
end

rho1 = counts/N;

figure
plot(uTimes,rho1)
%% Avalanche Stats

AvDurations = [];
AvSizes = [];

s = counts(1);
d = 1;

%slow:
for i=1:numel(uTimes)-1
    if uTimes(i+1)==uTimes(i)+1
        s = s+counts(i+1);
        d = d+1;
    else
        AvSizes = vertcat(AvSizes,s);
        AvDurations = vertcat(AvDurations,d);
        
        s = counts(i+1);
        d = 1;
    end
end

figure
logbin(AvSizes,1.6,'PDF','o');
hold on
logbin(AvDurations,1.6,'PDF','o');
legend('Sizes','Durations')
set(gca,'FontSize',12);
% xlabel('Avalanche Size, s','FontSize',12)
% ylabel('P(s)','FontSize',12)
                
                
num_Avalanches = numel(AvSizes);
disp([num2str(num_Avalanches),' avalanches recorded.'])
