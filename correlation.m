% Number of matter sites
%L = 12;
n_l = L/2; % Center site
times = 0:0.1:300; % Time evolution steps

% Operators in eigenbasis
occupationOps = cell(1, L);
for m = 1:L
    occupationOps{m} = M' * matterField{m} * M;
end

% Initialize storage for correlators
correlator_abs = zeros(L, length(times));

% Compute connected density-density correlator
for tInd = 1:length(times)
    t = times(tInd);
    psi_t = diag(exp(-1i * E * t)) * psi0Rotated;
    
    % Expectation values
    n_t = zeros(1, L);
    for m = 1:L
        n_t(m) = real(psi_t' * (occupationOps{m} * psi_t));
    end
    
    % Compute connected correlator
    for j = 1:L
        corr = real(psi_t' * (occupationOps{n_l} * occupationOps{j} * psi_t));
        correlator_abs(j, tInd) = abs(corr - n_t(n_l) * n_t(j));
    end
end

% Transpose for correct orientation in imagesc
correlator_abs = correlator_abs';  

% Create the density plot
figure;
imagesc(1:L, times, correlator_abs);
colorbar;
xlabel('Lattice Site');
ylabel('Time');
title('Connected Density-Density Correlator |C_j^c (t)|');
set(gca, 'YDir', 'normal'); % Ensures time runs from bottom to top
colormap(jet); % Adjust colormap for better visualization
