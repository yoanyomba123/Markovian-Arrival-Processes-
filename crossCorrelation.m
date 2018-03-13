function crosscorrel = crossCorrelation(arrivals, service_times, P)
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

crosscorrel = numerator/denominator
end

