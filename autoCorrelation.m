function autocorrel = autoCorrelation(s1, s2, s1_serv, s2_serv, P, n)
% Computes the autocorrelation of two samples
% inputs : s1 = sample 1
% .        s2 = sample 2
%          s1_serv = service time distribution of sample 1
%          s2_serv = service time distribution of sample 2
%          P = transition probability matrix
% .        n = nth time step

% acquire the eigen values and vector for respective transition probability
[V, D] = eig(P);

% compute the nstate transition probability table
% still faulty
P_nstep = inv(V')*(P - D)*V';

% Only obtain real parts of complex variable
P_nstep = real(P_nstep);


% acquire the steady state distribution for both samples
st_dist_s1 = V(:, 1);
st_dist_s2 = V(:, 1);
st_dis_s1 = abs(st_dist_s1) ./ sum(abs(st_dist_s1));
st_dis_s2 = abs(st_dist_s2) ./ sum(abs(st_dist_s2));

% obtain second moment of sample 1
second_moment = var(s1);

numerator = 0;
denominator = 0;
sum_service_time = 0;
sum_second_moment = 0;
for i = 1 : numel(s1)
   for j = 1 : nume(s2)
       numerator = numerator + st_dis_s1(i) * (P_nstep(i, j)^n - st_dis_s2(j)) * s1_serv(i) * s2_serv(i); 
       sum_service_time = sum_service_time + (st_dis_s1(i) * s1_serv(i));
       sum_second_moment = sum_second_moment + (st_dis_s1(i) * second_moment(i));
   end
end

denominator = (sum_service_time * sum_service_time) + sum_second_moment;
autocorrel = numerator/denominator;
end

