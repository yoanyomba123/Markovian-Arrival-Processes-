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
service_times = [2,1,2,1];

% Say, the transition table is
P = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0];

% Below is used to calculate cross-correlation between the inter-arrival
% and service time
[V D] = eig( P.' );

% st represents the State Transition Table
st = V(:,1).';

st = abs(st)./sum(abs(st));

lamb_neg_1 = sum(st.*(arrivals.^-1));

mean_service_time = sum(st.*service_times);

numerator = abs(sum(st.*(arrivals-lamb_neg_1).*(service_times-mean_service_time)));

denominator = abs((sum(st.*((arrivals.^-1-lamb_neg_1).^2))*sum(st.*((service_times-mean_service_time).^2)))^0.5);

corres = numerator/denominator


% compute the stationary distribution of the above example
syms a  a0 b b0 c c0 d d0 positive;
 
% assumption that all transition probabilities are less than one
assumeAlso([a0,b0,c0,d0] == 1/2);

% define the transition matrix bounding the transition probabilities
p = sym(zeros(4,4));
p(1,1:4) = [0,a0,0,1-a0];
p(2,1:4) = [1-b0,0,b0,0];
p(3,1:4) = [0,1-c0,0,c0];
p(4,1:4) = [d0,0,d0,0];
disp(p);

%% compute all possible analytical stationary distributions of the MC
% extract eigenvectors
[V, D] = eig(p);
% extract eigenvalues
ix = find(isAlways(diag(D) == 1,'unknown', 'error'));
diag(D(ix,ix))
%% extract stationary probabilities
Probability = V(:,ix);
%% 
% extact the stationary distributions

% compute the autocorrelation of th service times
