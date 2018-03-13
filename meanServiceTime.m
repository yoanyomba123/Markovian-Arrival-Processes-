function mean_service_time = meanServiceTime(service_times,P)
[V D] = eig( P.' );

% st represents the State Transition Table
st = V(:,1).';

st = abs(st)./sum(abs(st));


mean_service_time = sum(st.*service_times);
end

