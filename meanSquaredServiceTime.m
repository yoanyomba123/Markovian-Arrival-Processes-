function mean_squared_service_time = meanSquaredServiceTime(service_times,P)
[V D] = eig( P.' );

% st represents the State Transition Table
st = V(:,1).';

st = abs(st)./sum(abs(st));


mean_squared_service_time = sum(st.*(service_times.^2));
end


