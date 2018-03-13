function mean_arrival_time = meanArrivalTime(arrival_times,P)
[V D] = eig( P.' );

% st represents the State Transition Table
st = V(:,1).';

st = abs(st)./sum(abs(st));


mean_arrival_time = 1/sum(st.*(1./arrival_times));
end

