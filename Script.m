% Author : D Yoan L Mekontchou Yomba & Shiaho Song
% Date: 3/3/2018

%% workspace setup
clc; clear all; close all;

% add path
addpath(genpath('.'));
% Queineing process analysis


% Model exponentially distributed arrival service times
% number of arrival processes
N = 10;
Lambda = 5
% generate a vector of exponentially distributed service times
service = simexp(1, N, 20);

% generate n markovian arrival processes
arrival = poissrnd(Lambda,1,N);

%% Numerical Exampls
% Number of arrival processes
N = 4;
% markovian arrival processes (poisson)
arrivals = [0.1, 10, 0.1, 10];
% exponentially distributed mean service times
service_times = [1,2,1,2];

% Say, the transition table is
P = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0];

P_independent = (1/2) .* [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0];
% lack of auto-correlations
P_no_autocorr = (1/4) .* ones(4,4);
a = 0;
b = 3.37;
r = (b-a).*rand(1000,1) + a;
waiting_time = zeros(1,numel(r));
number_cust = zeros(1, numel(r));
mean_sys_time = zeros(1, numel(r));
mean_cust_sys = zeros(1, numel(r));

for i = 1:numel(r)
    traffic_intensity = 0.297 * r(i,1); 
    waiting_time(1,i) = (meanArrivalTime(arrivals,P) * meanSquaredServiceTime(service_times, P))/(2*(1-traffic_intensity));
    waiting_time(2,i) = (meanArrivalTime(arrivals,P_no_autocorr) * meanSquaredServiceTime(service_times, P_no_autocorr))/(2*(1-traffic_intensity));
    waiting_time(3,i) = (meanArrivalTime(arrivals,P_independent) * meanSquaredServiceTime(service_times, P_independent))/(2*(1-traffic_intensity));
    number_cust(1,i) = (meanArrivalTime(arrivals,P)^2 * meanSquaredServiceTime(service_times, P))/(2*(1-traffic_intensity));
    mean_sys_time(1,i) = (meanArrivalTime(arrivals,P) * meanSquaredServiceTime(service_times, P))/(2*(1-traffic_intensity)) + meanServiceTime(service_times, P);
    mean_cust_sys(1, i) = (meanArrivalTime(arrivals,P)^2 * meanSquaredServiceTime(service_times, P))/(2*(1-traffic_intensity)) + traffic_intensity;
end
p = 0.297 .* r;

figure;
scatter(p,waiting_time(1,:),'filled');
hold on;
scatter(p, waiting_time(2,:));
xlabel("p - Traffic Intensity w/ system in steay state")
ylabel("Waiting as a function of p");
title("Negative Correlation : Waiting Time in Queue vs Traffic Intensity - Observance of the effect of the presence of Auto-Correlation");

%%
figure;
scatter(p,waiting_time(1,:),'filled');
hold on;
scatter(p, waiting_time(2,:));
xlabel("p - Traffic Intensity w/ system in steay state")
ylabel("Waiting as a function of p");
title("Waiting Time in Queue vs Traffic Intensity - Observance of the effect of the presence of Auto-Correlation");
%%
figure;
stem(p,waiting_time(1,:));
xlabel("p - Traffic Intensity w/ system in steay state")
ylabel("Waiting as a function of p");
title("Waiting Time in Queue vs Traffic Intensity");

figure;
stem(p, number_cust(1, :));
xlabel("p - Traffic Intensity w/ system in steay state")
ylabel("Mean Number of Waiting customers in the buffer as a function of traffic intensity");
title("Number of Waiting Customers in the Queue vs Traffic Intensity");

figure;
stem(p, mean_sys_time(1, :));
xlabel("p - Traffic Intensity w/ system in steay state")
ylabel("Mean time spent in system as a function of traffic intensity");
title("Mean time spent in system vs the traffic intensity");

figure;
stem(p, mean_cust_sys(1, :));
xlabel("p - Traffic Intensity w/ system in steay state")
ylabel("Mean Number of Customers in system as a function of traffic intensity");
title("Mean Number of Customers in the system vs Traffic Intensity");
%% Observe Autorrelation behavior within arrival times as n approached inf
lambda = 0.9;
exp_service = exprnd(lambda, 1, 1500000);
figure;
autocorr(exp_service);

% sample auto-correlation
poisson_arrivals = poissrnd(lambda, 1, 1500000);
figure;title("Autocorrelation In Poisson distributed arrival ");
autocorr(poisson_arrivals);

