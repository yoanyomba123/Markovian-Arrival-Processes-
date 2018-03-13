function st_dist = steady_state_dist(sample, P)
% D Yoan L Mekontchou Yomba & Shihao Song
% Date 3/8/2018
% computes the steady state distribution of a same sample 

[V D] = eig( P.' );
% st represents the State Transition Table
st_dist = V(:,1).';

end

